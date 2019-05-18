/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.h"

// BAT
#include "BAT/BCMath.h"
#include "BAT/BCAux.h"

// ROOT
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TParameter.h"

GerdaFitter::GerdaFitter(json config) : config(config) {
    // open log file
    auto outdir = config["output-dir"].get<std::string>();
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/gerda-fitter-" + config["id"].get<std::string>() + "-";

    BCLog::SetLogLevelFile(BCLog::debug);
    BCLog::OpenLog(prefix + "output.log");
    BCLog::OutSummary("Saving results in " + outdir);
    BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

    this->SetName(config["id"]);

    int comp_idx = 0;
    // loop over files with data histograms
    for (auto& el : config["fit"].items()) {
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
            th->SetName(BCAux::SafeName(basename).c_str());

            dataset _current_ds;
            _current_ds.data = th;

            // get rebin factor
            auto rebin = elh.value().value("rebin-factor", 1);
            BCLog::OutDetail("using rebin factor = " + std::to_string(rebin));

            // eventually get a global value for the gerda-pdfs path
            auto gerda_pdfs_path = elh.value().value("gerda-pdfs", ".");

            BCLog::OutDebug("getting the requested pdfs");
            // loop over requested components
            for (auto& it : elh.value()["components"]) {

                auto prefix = it.value("prefix", gerda_pdfs_path);

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
                            path_to_part.erase(path_to_part.find_last_of('/'));
                            auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                            auto filename = prefix + "/" + p.key() + "/" + true_iso + "/" + "pdf-"
                                + volume + "-" + p.key() + "-" + i + ".root";
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
                    BCLog::OutDebug("building pdf for entry " + iso.key());

                    // it's a user defined file
                    if (it.contains("root-file")) {
                        auto obj_name = iso.value().value("hist-name", elh.key());
                        auto thh = this->GetFitComponent(it["root-file"].get<std::string>(), obj_name, _current_ds.data);
                        _current_ds.comp.insert({comp_idx, thh});

                    }
                    else { // look into gerda-pdfs database
                        if (iso.value()["isotope"].is_string()) {
                            auto comp = sum_parts(iso.value()["isotope"]);
                            comp->SetName((iso.key() + "_" + std::string(comp->GetName())).c_str());
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
                            comp->SetName((iso.key() + "_" + std::string(comp->GetName())).c_str());
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

                    // create a corresponding fit parameter
                    BCLog::OutDebug("adding model parameter '" + iso.key() + "' (\""
                        + iso.value().value("long-name", "") + "\" [" + iso.value().value("units", "")
                        + "]) in range = [" + std::to_string(iso.value()["parameter-range"][0].get<double>()) + ","
                        + std::to_string(iso.value()["parameter-range"][1].get<double>()) + "]");

                    this->AddParameter(
                        iso.key(),
                        iso.value()["parameter-range"][0].get<double>(),
                        iso.value()["parameter-range"][1].get<double>(),
                        iso.value().value("long-name", ""),
                        iso.value().value("units", "")
                    );
                    this->GetParameters().Back().SetPriorConstant();
                    comp_idx++;
                }
            }

            // eventually rebin data
            _current_ds.data->Rebin(rebin);

            // set fit range
            _current_ds.brange = {1, _current_ds.data->GetNbinsX()};
            if (elh.value().contains("fit-range")) {
                if (elh.value()["fit-range"].is_array()) {
                    if (!elh.value()["fit-range"][0].is_null()) {
                        _current_ds.brange.first = _current_ds.data->FindBin(elh.value()["fit-range"][0]);
                        BCLog::OutDebug("Setting lower fit bin range to " + std::to_string(_current_ds.brange.first));
                    }
                    if (!elh.value()["fit-range"][1].is_null()) {
                        _current_ds.brange.second = _current_ds.data->FindBin(elh.value()["fit-range"][1]);
                        BCLog::OutDebug("Setting upper fit bin range to " + std::to_string(_current_ds.brange.second));
                    }
                }
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
        for (int b = it.brange.first; b < it.brange.second; ++b) {
            // compute theoretical prediction for bin 'b'
            double pred = 0;
            for (auto& h : it.comp) pred += parameters[h.first]*h.second->GetBinContent(b);
            logprob += BCMath::LogPoisson(it.data->GetBinContent(b), pred);
        }
    }
    return logprob;
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
        tf.cd();
    }
}

TH1* GerdaFitter::GetFitComponent(std::string filename, std::string objectname, TH1* data) {
    TFile _tf(filename.c_str());
    if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
    auto obj = _tf.Get(objectname.c_str());
    if (!obj) throw std::runtime_error("could not find object '" + objectname + "' in file " + filename);
    // please do not delete it when the TFile goes out of scope
    if (obj->InheritsFrom(TH1::Class())) {
        auto _th = dynamic_cast<TH1*>(obj);
        if (_th->GetDimension() > 1) throw std::runtime_error("TH2/TH3 are not supported yet");
        if (_th->GetNbinsX() != data->GetNbinsX() or
            _th->GetDimension() != data->GetDimension() or
            _th->GetXaxis()->GetXmin() != data->GetXaxis()->GetXmin() or
            _th->GetXaxis()->GetXmax() != data->GetXaxis()->GetXmax()) {
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
        auto _th = new TH1D(obj->GetName(), obj->GetTitle(), data->GetNbinsX(),
            data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
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

void GerdaFitter::PrintExtendedFitSummary()
{
    BCLog::OutSummary("---------------------------------------------------");
    BCLog::OutSummary(Form("Fit summary for model \'%s\':", GetName().data()));
    BCLog::OutSummary(Form("  Number of parameters:  Npar = %i", GetNParameters()));
    BCLog::OutSummary("  Best fit parameters (global):");

    auto best = this->GetBestFitParameters();
    for (size_t i = 0; i < best.size(); ++i) {
        std::string line = Form("    %-*s : %.*g ± %.*g",
                this->GetMaximumParameterNameLength(best.size() > this->GetNParameters()),
                this->GetVariable(i).GetName().data(),
                2, best[i],
                2, this->GetBestFitParameterErrors()[i]
        );
        if (this->GetParameter(i).Fixed()) line += " (fixed)";
        BCLog::OutSummary(line);
    }

    BCLog::OutSummary("---------------------------------------------------");
}
