/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.hh"

// BAT
#include "BAT/BCMath.h"
#include "BAT/BCAux.h"
#include "BAT/BCTF1Prior.h"
#include "BAT/BCTH1Prior.h"

// ROOT
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TParameter.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TParameter.h"

GerdaFitter::GerdaFitter(json outconfig) : config(outconfig) {
    // open log file
    auto outdir = config["output-dir"].get<std::string>();
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/gerda-fitter-" + config["id"].get<std::string>() + "-";

    BCLog::SetLogLevelFile(BCLog::debug);
    BCLog::OpenLog(prefix + "output.log");
    BCLog::OutSummary("Saving results in " + outdir);
    BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

    this->SetName(config["id"]);

    /*
     * create fit parameters
     */

    for (auto& el : config["fit"]["parameters"].items()) {
        // do we already have a parameter with the same name?
        bool already_exists = false;
        for (unsigned int idx = 0; idx < this->GetNParameters(); ++idx) {
            if (this->GetParameters().At(idx).GetName() == el.key()) {
                already_exists = true;
                break;
            }
        }

        if (already_exists) {
            BCLog::OutDetail("model parameter '" + el.key() + "' already exists, skipping");
        }
        else {
            BCLog::OutDetail("adding model parameter '" + el.key() + "' (\""
                + el.value().value("long-name", "") + "\" [" + el.value().value("units", "")
                + "]) in range = [" + std::to_string(el.value()["range"][0].get<double>()) + ","
                + std::to_string(el.value()["range"][1].get<double>()) + "]");

            this->GetParameters().Add(
                el.key(),
                el.value()["range"][0].get<double>(),
                el.value()["range"][1].get<double>(),
                el.value().value("long-name", ""),
                el.value().value("units", "")
            );
        }

        // assign prior
        if (el.value().contains("prior")) {
            auto& prior_cfg = el.value()["prior"];
            TH1* prior_hist = nullptr;
            BCPrior* prior = nullptr;
            std::string expr;

            if (prior_cfg.contains("TFormula") and prior_cfg.contains("histogram")) {
                throw std::runtime_error("please choose either \"TFormula\" or \"histogram\" for parameter " + el.key());
            }
            // a ROOT histogram is given
            if (prior_cfg.contains("histogram")) {
                expr = prior_cfg["histogram"].get<std::string>();
                if (expr.find(':') == std::string::npos) {
                    throw std::runtime_error("invalid \"histogram\" format for parameter " + el.key());
                }
                auto filename = expr.substr(0, expr.find_first_of(':'));
                auto objname = expr.substr(expr.find_first_of(':')+1, std::string::npos);
                TFile _tff(filename.c_str());
                if (!_tff.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
                auto obj = _tff.Get(objname.c_str());
                if (!obj) throw std::runtime_error("could not find object '" + objname + "' in file " + filename);
                if (obj->InheritsFrom(TH1::Class())) {
                    prior_hist = dynamic_cast<TH1*>(obj);
                    // the following is a workaround for BAT's bad implementation of the BCTH1Prior constructor
                    prior_hist->SetDirectory(nullptr);
                    _tff.Close();
                }
                else throw std::runtime_error("object '" + objname + "' in file " + filename + " is not an histogram");
                BCLog::OutDebug("found histogram prior '" + expr + "' for parameter '" + el.key() + "'");
            }
            // is a tformula TODO: i cannot make it work, it segfaults
            else if (prior_cfg.contains("TFormula")) {
                TF1* _tformula = nullptr;
                expr = prior_cfg["TFormula"].get<std::string>();
                // if there's a list of parameters after
                if (expr.find(':') != std::string::npos) {
                    auto formula = expr.substr(0, expr.find_first_of(':'));
                    auto parlist = expr.substr(expr.find_first_of(':')+1, std::string::npos);
                    std::vector<double> parlist_vec;
                    _tformula = new TF1(
                        "f1_prior", formula.c_str(),
                        el.value()["range"][0].get<double>(),
                        el.value()["range"][1].get<double>()
                    );

                    if (!_tformula->IsValid()) throw std::runtime_error("invalid TFormula given");

                    // eventually set parameters, if any
                    while (!parlist.empty()) {
                        auto val = std::stod(parlist.substr(0, parlist.find_first_of(',')));
                        parlist_vec.push_back(val);

                        if (parlist.find(',') == std::string::npos) parlist = "";
                        else parlist.erase(0, parlist.find_first_of(',')+1);
                    }
                    if ((int)parlist_vec.size() != _tformula->GetNpar()) {
                        throw std::runtime_error("number of values specified in '"
                                + el.key() + "' prior TFormula does not match number of TFormula parameters");
                    }

                    for (size_t j = 0; j < parlist_vec.size(); ++j) _tformula->SetParameter(j, parlist_vec[j]);
                }
                // it's just the tformula
                else {
                    _tformula = new TF1(
                        "f1_prior", expr.c_str(),
                        el.value()["range"][0].get<double>(),
                        el.value()["range"][1].get<double>()
                    );
                    if (_tformula->GetNpar() > 0) {
                        throw std::runtime_error("TFormula specified for parameter '"
                                + el.key() + "' is parametric but no parameters were specified");
                    }
                }
                BCLog::OutDebug("found TFormula prior '" + expr + "' for parameter '" + el.key() + "'");

                // I did not manage to use BCTF1Prior from BAT, it segfaults for misterious reasons
                // so I just make a temporary histogram out of the TFormula
                prior_hist = new TH1D("h1_prior", "h1_prior", 1000,
                    el.value()["range"][0].get<double>(),
                    el.value()["range"][1].get<double>()
                );
                for (int j = 1; j <= prior_hist->GetNbinsX(); j++) {
                    prior_hist->SetBinContent(j, _tformula->Eval(prior_hist->GetBinCenter(j)));
                }
                prior_hist->SetDirectory(nullptr);
                BCLog::OutDebug("produced temporary prior 1000-bin histogram for parameter '" + el.key() + "'");
            }
            else throw std::runtime_error("supported prior types: 'histogram', 'TFormula'");

            prior = new BCTH1Prior(prior_hist);
            delete prior_hist;
            BCLog::OutDetail("assigned prior '" + expr + "' to parameter '" + el.key() + "'");

            if (prior != nullptr) {
                if (!prior->IsValid()) {
                    throw std::runtime_error("invalid prior set for parameter " + el.key());
                }
                this->GetParameters().Back().SetPrior(prior);
            }
            else {
                this->GetParameters().Back().SetPriorConstant();
                BCLog::OutDetail("assigned uniform prior to parameter '" + el.key() + "'");
            }
        }
        else {
            this->GetParameters().Back().SetPriorConstant();
            BCLog::OutDetail("assigned uniform prior to parameter '" + el.key() + "'");
        }
    }

    /*
     * read in histograms
     */

    // loop over files with data histograms
    for (auto& el : config["fit"]["theoretical-expectations"].items()) {
        BCLog::OutDebug("opening data file " + el.key());
        TFile _tf(el.key().c_str());
        if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + el.key());

        // loop over requested histograms in file
        for (auto& elh : el.value().items()) {
            BCLog::OutDebug("getting data histogram '" + elh.key() + "' from " + el.key());
            auto th = dynamic_cast<TH1*>(_tf.Get(elh.key().c_str()));
            if (!th) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + el.key());
            th->SetDirectory(nullptr);
            auto basename = el.key().substr(
                    el.key().find_last_of('/')+1,
                    el.key().find_last_of('.') - el.key().find_last_of('/') - 1);
            basename += "_" + std::string(th->GetName());
            th->SetName(this->SafeROOTName(basename).c_str());

            dataset _current_ds;
            _current_ds.data = th;

            // get rebin factor
            auto rebin = elh.value().value("rebin-factor", 1);
            BCLog::OutDetail("using rebin factor = " + std::to_string(rebin) + " for dataset '" + elh.key() + "'");

            // eventually get a global value for the gerda-pdfs path
            auto gerda_pdfs_path = elh.value().value("gerda-pdfs", ".");

            BCLog::OutDebug("getting the requested PDFs");
            // loop over requested components
            for (auto& it : elh.value()["components"]) {

                prefix = it.value("prefix", gerda_pdfs_path);

                /* START INTERMEZZO */
                // utility to sum over the requested parts (with weight) given isotope
                auto sum_parts = [&](std::string i) {
                    std::string true_iso = i;
                    if (i.find('-') != std::string::npos) true_iso = i.substr(0, i.find('-'));

                    TH1* sum = nullptr;
                    if (it["part"].is_object()) {
                        // compute sum of weights
                        double sumw = 0;
                        for (auto& p : it["part"].items()) sumw += p.value().get<double>();

                        for (auto& p : it["part"].items()) {
                            // get volume name
                            auto path_to_part = prefix + "/" + p.key();
                            if (path_to_part.back() == '/') path_to_part.pop_back();
                            auto part = path_to_part.substr(path_to_part.find_last_of('/')+1);
                            path_to_part.erase(path_to_part.find_last_of('/'));
                            auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                            auto filename = prefix + "/" + p.key() + "/" + true_iso + "/" + "pdf-"
                                + volume + "-" + part + "-" + i + ".root";
                            BCLog::OutDebug("opening file " + filename);
                            BCLog::OutDebug("summing object '" + elh.key() + " with weight "
                                    + std::to_string(p.value().get<double>()/sumw));
                            // get histogram
                            auto thh = this->GetFitComponent(filename, elh.key(), _current_ds.data);
                            // add it with weight
                            if (!sum) {
                                sum = thh;
                                sum->SetDirectory(nullptr); // please do not delete it when the TFile goes out of scope
                                sum->Scale(p.value().get<double>()/sumw);
                            }
                            else sum->Add(thh, p.value().get<double>()/sumw);
                        }
                        return sum;
                    }
                    else if (it["part"].is_string()) {
                        // get volume name
                        auto path_to_part = prefix + "/" + it["part"].get<std::string>();
                        if (path_to_part.back() == '/') path_to_part.pop_back();
                        auto part = path_to_part.substr(path_to_part.find_last_of('/')+1);
                        path_to_part.erase(path_to_part.find_last_of('/'));
                        auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                        auto filename = prefix + "/" + it["part"].get<std::string>() + "/" + true_iso + "/" + "pdf-"
                            + volume + "-" + part + "-" + i + ".root";
                        BCLog::OutDebug("getting object '" + elh.key() + "' in file " + filename);
                        // get histogram
                        auto thh = this->GetFitComponent(filename, elh.key(), _current_ds.data);
                        return thh;
                    }
                    else throw std::runtime_error("unexpected 'part' value found in [\"fit\"][\""
                            + el.key() + "\"][\"" + elh.key() + "\"]");
                };
                /* END INTERMEZZO */

                // loop over requested isotopes on the relative part
                for (auto& iso : it["components"].items()) {
                    BCLog::OutDetail("adding a PDF for the '" + elh.key() + "' dataset with '"
                        + iso.key() + "' as scaling parameter (debug mode for details)");
                    // throw exception if component name is not in list
                    bool exists = false;
                    int comp_idx = 0;
                    for (unsigned int idx = 0; idx < this->GetNParameters(); ++idx) {
                        if (this->GetParameters().At(idx).GetName() == iso.key()) {
                            comp_idx = idx;
                            exists = true;
                            break;
                        }
                    }
                    if (!exists) throw std::runtime_error(
                        "fit parameter '" + iso.key() + "' not found, is it defined in \"parameters\"?"
                    );

                    // it's a user defined file
                    if (it.contains("root-file")) {
                        auto obj_name = iso.value().value("hist-name", elh.key());
                        auto thh = this->GetFitComponent(it["root-file"].get<std::string>(), obj_name, _current_ds.data);
                        _current_ds.comp.insert({comp_idx, thh});

                    }
                    else if (iso.value().contains("TFormula")) {
                        if (!iso.value()["TFormula"].is_string()) {
                            throw std::runtime_error("The \"TFormula\" key must be a string");
                        }
                        TF1 _tfunc(
                            "func", iso.value()["TFormula"].get<std::string>().c_str(),
                            _current_ds.data->GetXaxis()->GetXmin(),
                            _current_ds.data->GetXaxis()->GetXmax()
                        );
                        if (_tfunc.GetNpar() > 0) throw std::runtime_error("TFormula parameters are not supported yet");
                        auto thh = new TH1D(
                            iso.key().c_str(), iso.key().c_str(),
                            _current_ds.data->GetNbinsX(),
                            _current_ds.data->GetXaxis()->GetXmin(),
                            _current_ds.data->GetXaxis()->GetXmax()
                        );
                        for (int b = 1; b < thh->GetNbinsX(); ++b) {
                            thh->SetBinContent(b, _tfunc.Eval(thh->GetBinCenter(b)));
                        }
                        thh->SetDirectory(nullptr);
                        thh->SetName(this->SafeROOTName(iso.key()).c_str());
                        _current_ds.comp.insert({comp_idx, thh});
                    }
                    else { // look into gerda-pdfs database
                        if (iso.value()["isotope"].is_string()) {
                            auto comp = sum_parts(iso.value()["isotope"]);
                            comp->SetName(this->SafeROOTName(iso.key() + "_" + std::string(comp->GetName())).c_str());
                            _current_ds.comp.insert({comp_idx, comp});
                        }
                        else if (iso.value()["isotope"].is_object()) {
                            double sumwi = 0;
                            for (auto& i : iso.value()["isotope"].items()) sumwi += i.value().get<double>();

                            TH1* comp = nullptr;
                            for (auto& i : iso.value()["isotope"].items()) {
                                BCLog::OutDebug("scaling pdf for " + i.key() + " by a factor "
                                        + std::to_string(i.value().get<double>()/sumwi));
                                if (!comp) {
                                    comp = sum_parts(i.key());
                                    comp->Scale(i.value().get<double>()/sumwi);
                                }
                                else comp->Add(sum_parts(i.key()), i.value().get<double>()/sumwi);

                            }
                            comp->SetName(this->SafeROOTName(iso.key() + "_" + std::string(comp->GetName())).c_str());
                            _current_ds.comp.insert({comp_idx, comp});
                        }
                        else throw std::runtime_error("unexpected entry " + iso.value()["isotope"].dump() + "found in [\"fit\"][\""
                                + el.key() + "\"][\"" + elh.key() + "\"][\"components\"][\"" + iso.key() + "\"][\"isotope\"]");

                        if (!_current_ds.comp[comp_idx]) {
                            throw std::runtime_error("invalid pointer found in component list at position " + std::to_string(comp_idx));
                        }
                    }

                    // eventually rebin
                    _current_ds.comp[comp_idx]->Rebin(rebin);
                }
            }

            // eventually rebin data
            _current_ds.data->Rebin(rebin);

            // set fit range
            if (elh.value().contains("fit-range")) {
                if (elh.value()["fit-range"].is_array()) {
                    if (elh.value()["fit-range"][0].is_array()) {
                        BCLog::OutDebug("\"fit-range\" is an array of arrays");
                        for (auto& r : elh.value()["fit-range"]) {
                            if (!r.is_array()) throw std::runtime_error("non-array member detected in \"fit-range\"");
                            if (!r[0].is_number() or !r[1].is_number()) {
                                throw std::runtime_error("not-a-number value detected in \"fit-range\"");
                            }
                            _current_ds.brange.push_back({_current_ds.data->FindBin(r[0]), _current_ds.data->FindBin(r[1])});
                            BCLog::OutDetail("Adding fit range [" +
                                std::to_string(_current_ds.brange.back().first) + ", " +
                                std::to_string(_current_ds.brange.back().second) + "] (bins)"
                            );
                        }
                    }
                    if (elh.value()["fit-range"][0].is_number()) {
                        BCLog::OutDebug("\"fit-range\" is an array of numbers");
                        _current_ds.brange.push_back({
                            _current_ds.data->FindBin(elh.value()["fit-range"][0]),
                            _current_ds.data->FindBin(elh.value()["fit-range"][1])
                        });
                        BCLog::OutDetail("Setting fit range to [" +
                            std::to_string(_current_ds.brange.back().first) + ", " +
                            std::to_string(_current_ds.brange.back().second) + "] (bins)"
                        );
                    }
                }
                else throw std::runtime_error("wrong \"fit-range\" format, array expected");
            }
            else {
                _current_ds.brange.push_back({1, _current_ds.data->GetNbinsX()});
            }

            // now sanity checks on the range
            double _last = - std::numeric_limits<double>::infinity();
            for (auto& r : _current_ds.brange) {
                if (r.first > _last) _last = r.first;
                else throw std::runtime_error("illegal ranges in \"fit-range\" detected, did you specify them in ascending order?");
                if (r.second > _last) _last = r.second;
                else throw std::runtime_error("illegal ranges in \"fit-range\" detected, did you specify them in ascending order?");
            }

            data.push_back(_current_ds);

            BCLog::OutDebug("data and pdf components got so far:");
            this->DumpData();
        }
    }
}

GerdaFitter::~GerdaFitter() {
    for (auto& h : data) {
        delete h.data;
        for (auto& hh : h.comp) delete hh.second;
    }
}

double GerdaFitter::LogLikelihood(const std::vector<double>& parameters) {
    double logprob = 0;
    // loop over datasets
    for (auto& it : data) {
        for (auto& r : it.brange) {
            for (int b = r.first; b < r.second; ++b) {
                // compute theoretical prediction for bin 'b'
                double pred = 0;
                for (auto& h : it.comp) pred += parameters[h.first]*h.second->GetBinContent(b);
                logprob += BCMath::LogPoisson(it.data->GetBinContent(b), pred);
            }
        }
    }
    return logprob;
}

TH1* GerdaFitter::GetFitComponent(std::string filename, std::string objectname, TH1* dataformat) {
    TFile _tf(filename.c_str());
    if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
    auto obj = _tf.Get(objectname.c_str());
    if (!obj) throw std::runtime_error("could not find object '" + objectname + "' in file " + filename);
    // please do not delete it when the TFile goes out of scope
    if (obj->InheritsFrom(TH1::Class())) {
        auto _th = dynamic_cast<TH1*>(obj);
        if (_th->GetDimension() > 1) throw std::runtime_error("TH2/TH3 are not supported yet");
        if (_th->GetNbinsX() != dataformat->GetNbinsX() or
            _th->GetDimension() != dataformat->GetDimension() or
            _th->GetXaxis()->GetXmin() != dataformat->GetXaxis()->GetXmin() or
            _th->GetXaxis()->GetXmax() != dataformat->GetXaxis()->GetXmax()) {
            throw std::runtime_error("histogram '" + objectname + "' in file " + filename
                + " and corresponding data histogram do not have the same number of bins and/or same ranges");
        }
        TParameter<Long64_t>* _nprim = nullptr;
        if (objectname.substr(0, 3) == "M1_") {
            _nprim = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesEdep"));
        }
        else if (objectname.substr(0, 3) == "M2_") {
            _nprim = dynamic_cast<TParameter<Long64_t>*>(_tf.Get("NumberOfPrimariesCoin"));
        }
        long long int nprim = (_nprim) ? _nprim->GetVal() : 1;
        _th->Scale(1./nprim);

        _th->SetDirectory(nullptr);
        return _th;
    }
    else if (obj->InheritsFrom(TF1::Class())) {
        auto _th = new TH1D(obj->GetName(), obj->GetTitle(), dataformat->GetNbinsX(),
            dataformat->GetXaxis()->GetXmin(), dataformat->GetXaxis()->GetXmax());
        for (int b = 1; b < _th->GetNbinsX(); ++b) {
            _th->SetBinContent(b, dynamic_cast<TF1*>(obj)->Eval(_th->GetBinCenter(b)));
        }
        _th->SetDirectory(nullptr);
        return _th;
    }
    else {
        throw std::runtime_error("object '" + objectname + "' in file " + filename + " isn't of type TH1 or TF1");
    }
}

void GerdaFitter::SetIntegrationProperties(json j) {
    if (j.is_null()) return;

    this->SetIntegrationMethod(j.value("method", BCIntegrate::BCIntegrationMethod::kIntDefault));
    this->SetCubaIntegrationMethod(j.value("cuba-method", BCIntegrate::BCCubaMethod::kCubaDefault));

    auto method = j.value("method", "kIntDefault");
    BCLog::OutDebug("Integration method set: " + method);
    auto cubamethod = j.value("cuba-method", "kCubaDefault");
    BCLog::OutDebug("Cuba integration method set: " + cubamethod);

    if (j["integrator-settings"].is_object()) {
        auto& jj = j["integrator-settings"];
        if (jj[method].is_object()) {
            if (this->GetIntegrationMethod() != BCIntegrate::BCIntegrationMethod::kIntCuba) {
                if (jj[method]["niter-min"].is_number()) {
                    this->SetNIterationsMin(jj[method]["niter-min"].get<long long>());
                }
                if (jj[method]["niter-max"].is_number()) {
                    this->SetNIterationsMax(jj[method]["niter-max"].get<long long>());
                }
                if (jj[method]["rel-precision"].is_number()) {
                    this->SetRelativePrecision(jj[method]["rel-precision"].get<double>());
                }
                if (jj[method]["abs-precision"].is_number()) {
                    this->SetAbsolutePrecision(jj[method]["abs-precision"].get<double>());
                }
            }
            else {
                if (jj[method][cubamethod].is_object()) {
                    if (jj[method][cubamethod]["niter-min"].is_number()) {
                        this->SetNIterationsMin(jj[method][cubamethod]["niter-min"].get<long long>());
                    }
                    if (jj[method][cubamethod]["niter-max"].is_number()) {
                        this->SetNIterationsMax(jj[method][cubamethod]["niter-max"].get<long long>());
                    }
                    if (jj[method][cubamethod]["rel-precision"].is_number()) {
                        this->SetRelativePrecision(jj[method][cubamethod]["rel-precision"].get<double>());
                    }
                    if (jj[method][cubamethod]["abs-precision"].is_number()) {
                        this->SetAbsolutePrecision(jj[method][cubamethod]["abs-precision"].get<double>());
                    }
                }
                auto& jjj = jj[method][cubamethod];

                auto set_base_props = [&jjj](BCCubaOptions::General& m) {
                    if (jjj["ncomp"].is_number()) m.ncomp = jjj["ncomp"].get<int>();
                    if (jjj["flags"].is_number()) m.flags = jjj["flags"].get<int>();
                };

                switch (this->GetCubaIntegrationMethod()) {

                    case BCIntegrate::BCCubaMethod::kCubaVegas : {
                        auto o = this->GetCubaVegasOptions();
                        set_base_props(o);
                        if (jjj["nstart"]   .is_number()) o.nstart    = jjj["nstart"]   .get<int>();
                        if (jjj["nincrease"].is_number()) o.nincrease = jjj["nincrease"].get<int>();
                        if (jjj["nbatch"]   .is_number()) o.nbatch    = jjj["nbatch"]   .get<int>();
                        if (jjj["gridno"]   .is_number()) o.gridno    = jjj["gridno"]   .get<int>();
                        this->SetCubaOptions(o);
                    }
                    case BCIntegrate::BCCubaMethod::kCubaSuave : {
                        auto o = this->GetCubaSuaveOptions();
                        set_base_props(o);
                        if (jjj["nmin"]    .is_number()) o.nmin     = jjj["nmin"]    .get<int>();
                        if (jjj["nnew"]    .is_number()) o.nnew     = jjj["nnev"]    .get<int>();
                        if (jjj["flatness"].is_number()) o.flatness = jjj["flatness"].get<double>();
                        this->SetCubaOptions(o);
                    }
                    case BCIntegrate::BCCubaMethod::kCubaDivonne : {
                        auto o = this->GetCubaDivonneOptions();
                        set_base_props(o);
                        if (jjj["key1"]        .is_number()) o.key1         = jjj["key3"]        .get<int>();
                        if (jjj["key2"]        .is_number()) o.key2         = jjj["key2"]        .get<int>();
                        if (jjj["key3"]        .is_number()) o.key3         = jjj["key1"]        .get<int>();
                        if (jjj["maxpass"]     .is_number()) o.maxpass      = jjj["maxpass"]     .get<int>();
                        if (jjj["border"]      .is_number()) o.border       = jjj["border"]      .get<double>();
                        if (jjj["maxchisq"]    .is_number()) o.maxchisq     = jjj["maxchisq"]    .get<double>();
                        if (jjj["mindeviation"].is_number()) o.mindeviation = jjj["mindeviation"].get<double>();
                        this->SetCubaOptions(o);
                    }
                    case BCIntegrate::BCCubaMethod::kCubaCuhre : {
                        auto o = this->GetCubaCuhreOptions();
                        set_base_props(o);
                        if (jjj["key"].is_number()) o.key = jjj["key"].get<int>();
                        this->SetCubaOptions(o);
                    }
                    case BCIntegrate::BCCubaMethod::kCubaDefault : break;
                    case BCIntegrate::BCCubaMethod::NCubaMethods : break;
                }
            }
        }
    }
}

void GerdaFitter::SaveHistograms(std::string filename) {
    TFile tf(filename.c_str(), "recreate");
    for (auto& it : data) {
        tf.mkdir(it.data->GetName());
        tf.cd(it.data->GetName());
        it.data->Write("fitted_data");

        TH1* sum = nullptr;
        for (auto& h : it.comp) {
            h.second->Scale(this->GetBestFitParameters()[h.first]);
            // compute total model
            if (!sum) sum = dynamic_cast<TH1*>(h.second->Clone());
            else      sum->Add(h.second);
            h.second->Write();
        }
        sum->SetTitle("total_model");
        sum->Write("total_model");

        // write fit range
        TParameter<double> range_low("fit_range_lower", it.data->GetBinLowEdge(it.brange[0].first));
        TParameter<double> range_upp("fit_range_upper", it.data->GetBinLowEdge(it.brange.back().second)
                + it.data->GetBinWidth(it.brange.back().second));
        range_low.Write();
        range_upp.Write();
        tf.cd();
    }
}

void GerdaFitter::DumpData() {
    for (auto& it : data) {
        BCLog::OutDebug(it.data->GetName());
        size_t i = 0;
        for (auto& comp : it.comp) {
            auto msg = (i == it.comp.size()-1 ? " └─ [" : " ├─ [")
                + std::to_string(comp.first) + "] " + std::string(comp.second->GetName());
            BCLog::OutDebug(msg);
            i++;
        }
    }
}

double GerdaFitter::GetFastPValue(const std::vector<double>& parameters, long niter) {

    BCLog::OutSummary("Computing fast (corrected) p-value with " + std::to_string(niter) + " iterations");

    if (parameters.size() != this->GetParameters().Size()) {
        throw std::runtime_error("GetFastPValue: input number of parameters does not match the number of defined model parameters");
    }

    std::vector<unsigned> observed;
    std::vector<double>   expected;
    int nbins = 0;

    for (auto& it : data) {
        // compute total model
        TH1* sum = nullptr;
        for (auto& h : it.comp) {
            h.second->Scale(parameters[h.first]);
            // compute total model
            if (!sum) sum = dynamic_cast<TH1*>(h.second->Clone());
            else      sum->Add(h.second);
        }
        for (auto& r : it.brange) {
            for (int b = r.first; b < r.second; ++b) {
                observed.push_back(it.data->GetBinContent(b));
                expected.push_back(sum->GetBinContent(b));
                nbins++;
            }
        }
    }

    TRandom3 rnd(0);
    double pvalue = BCMath::CorrectPValue(
        BCMath::FastPValue(observed, expected, niter, rnd.GetSeed()),
        this->GetNFreeParameters(), nbins
    );

    return pvalue;
}

void GerdaFitter::PrintOptimizationSummary() {
    // ┌─┬┐
    // │ ││
    // ├─┼┤
    // └─┴┘
    BCLog::OutSummary(Form("Optimization summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Best fit parameters (global mode):");
    auto logl_best = this->LogLikelihood(this->GetBestFitParameters());
    BCLog::OutSummary("  LogLikelihood value at global mode: " + std::to_string(logl_best));

    auto best = this->GetBestFitParameters();

    std::string line;
    int maxnamelength = this->GetMaximumParameterNameLength(false);
    line = "    ┌"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┬───────────────────┐";
    BCLog::OutSummary(line);
    BCLog::OutSummary(Form("    │ %-*s │ global mode       │", maxnamelength+5, "parameters"));
    line = "    ├"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┼───────────────────┤";
    BCLog::OutSummary(line);
    for (size_t i = 0; i < best.size(); ++i) {
        auto val = Form("%.*g", 2, best[i]);
        auto err = Form("%.*g", 2, this->GetBestFitParameterErrors()[i]);
        line = Form("    │ [%2i] %-*s │ %7s ± %-7s │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                val, err
        );
        if (this->GetParameter(i).Fixed()) line += " (fixed)";
        BCLog::OutSummary(line);
    }
    line = "    └"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┴───────────────────┘";
    BCLog::OutSummary(line);
}

void GerdaFitter::PrintShortMarginalizationSummary() {
    // ┌─┬┐
    // │ ││
    // ├─┼┤
    // └─┴┘
    BCLog::OutSummary(Form("Marginalization summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Parameters values (marginalized mode):");

    std::string line;
    int maxnamelength = this->GetMaximumParameterNameLength(false);

    line = "    ┌"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┬─────────────────────────────┐";
    BCLog::OutSummary(line);
    BCLog::OutSummary(Form("    │ %-*s │ marg. mode ± 1σ             │", maxnamelength+5, "parameters"));
    line = "    ├"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┼─────────────────────────────┤";
    BCLog::OutSummary(line);
    for (size_t i = 0; i < this->GetNParameters(); ++i) {
        auto val  = Form("%.*g", 2, this->GetMarginalized(i).GetLocalMode());
        auto err1 = Form("%.*g", 2, this->GetMarginalized(i).GetQuantile(0.16));
        auto err2 = Form("%.*g", 2, this->GetMarginalized(i).GetQuantile(1-0.16));
        line = Form("    │ [%2i] %-*s │ %7s + %-7s - %-7s │",
                int(i),
                maxnamelength,
                this->GetVariable(i).GetName().data(),
                val, err2, err1
        );
        if (this->GetParameter(i).Fixed()) line += " (fixed)";
        BCLog::OutSummary(line);
    }
    line = "    └"; for (int i = 0; i < maxnamelength+7; ++i) line += "─"; line += "┴─────────────────────────────┘";
    BCLog::OutSummary(line);
}

std::string GerdaFitter::SafeROOTName(const std::string orig) {
    TString torig(orig);

    for (auto& c : {'.', '-', '/', ':', '|', '+'}) torig.ReplaceAll(c, '_');
    for (auto& c : {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}) {
        if (torig[0] == c) torig[0] = 'N';
    }

    return std::string(torig.Data());
}
