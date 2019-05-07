/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.h"

// BAT
#include "BAT/BCLog.h"
#include "BAT/BCMath.h"

// ROOT
#include "TFile.h"
#include "TH1.h"

GerdaFitter::GerdaFitter(json config) {
    BCLog::SetLogLevelFile(BCLog::detail);
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
            th->SetName((std::string(th->GetName()) + "_" + basename).c_str());
            std::pair<TH1*, std::map<int, TH1*>> _pair = {th, {}};

            BCLog::OutDebug("getting the requested pdfs");
            // loop over requested components
            for (auto& it : elh.value()) {
                auto prefix = it.value("prefix", ".");

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
                            // open pdf file
                            TFile _tff(filename.c_str());
                            if (!_tff.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
                            BCLog::OutDebug("summing object '" + elh.key() + " with weight "
                                    + std::to_string(p.value().get<double>()/sumw));
                            // get histogram
                            auto thh = dynamic_cast<TH1*>(_tff.Get(elh.key().c_str()));
                            if (!thh) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + filename);
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
                        path_to_part.erase(path_to_part.find_last_of('/'));
                        auto volume = path_to_part.substr(path_to_part.find_last_of('/')+1);
                        auto filename = prefix + "/" + it["part"].get<std::string>() + "/" + i + "/" + "pdf-"
                            + volume + "-" + it["part"].get<std::string>() + "-" + i + ".root";
                        BCLog::OutDebug("getting object '" + elh.key() + "' in file " + filename);
                        // open pdf file
                        TFile _tff(filename.c_str());
                        if (!_tff.IsOpen()) throw std::runtime_error("invalid ROOT file: " + filename);
                        // get histogram
                        auto thh = dynamic_cast<TH1*>(_tff.Get(elh.key().c_str()));
                        if (!thh) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + filename);
                        thh->SetDirectory(nullptr); // please do not delete it when the TFile goes out of scope
                        return thh;
                    }
                    else throw std::runtime_error("unexpected 'part' value found in [\"fit\"][\""
                            + el.key() + "\"][\"" + elh.key() + "\"]");
                };

                // loop over requested isotopes on the relative part
                for (auto& iso : it["components"].items()) {
                    BCLog::OutDebug("building pdf for entry " + iso.key());

                    if (iso.value()["isotope"].is_string()) {
                        auto comp = sum_parts(iso.value()["isotope"]);
                        comp->SetName((std::string(comp->GetName()) + "_" + iso.key()).c_str());
                        _pair.second.insert({comp_idx, comp});
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
                        comp->SetName((std::string(comp->GetName()) + "_" + iso.key()).c_str());
                        _pair.second.insert({comp_idx, comp});
                    }
                    else throw std::runtime_error("unexpected entry " + iso.value()["isotope"].dump() + "found in [\"fit\"][\""
                            + el.key() + "\"][\"" + elh.key() + "\"][\"compoents\"][\"" + iso.key() + "\"][\"isotope\"]");

                    // create a corresponding fit parameter
                    BCLog::OutDebug("adding model parameter '" + iso.key() + "' (\""
                        + iso.value().value("long-name", "") + "\" [" + iso.value().value("units", "")
                        + "]) in range = [" + std::to_string(iso.value()["range"][0].get<double>()) + ","
                        + std::to_string(iso.value()["range"][1].get<double>()) + "]");

                    this->AddParameter(
                        iso.key(),
                        iso.value()["range"][0].get<double>(),
                        iso.value()["range"][1].get<double>(),
                        iso.value().value("long-name", ""),
                        iso.value().value("units", "")
                    );
                    this->GetParameters().Back().SetPriorConstant();
                    comp_idx++;
                }
            }
            data.insert(_pair);

            BCLog::OutDebug("data and pdf components got so far:");
            this->DumpData();
        }
    }
}

GerdaFitter::~GerdaFitter() {
    for (auto& h : data) {
        delete h.first;
        for (auto& hh : h.second) delete hh.second;
    }
}

double GerdaFitter::LogLikelihood(const std::vector<double>& parameters) {
    double logprob = 0;
    double f = 0;
    for (auto& it : data) {
        for (int b = 1; b < it.first->GetNbinsX(); ++b) {
            for (auto& h : it.second) {
                f += parameters[h.first]*h.second->GetBinContent(b);
            }
            logprob += BCMath::LogPoisson(it.first->GetBinContent(b), f);
        }
    }
    return logprob;
}

void GerdaFitter::DumpData() {
    for (auto& it : data) {
        BCLog::OutDebug(it.first->GetName());
        size_t i = 0;
        for (auto& comp : it.second) {
            auto msg = (i == it.second.size()-1 ? " └─ [" : " ├─ [")
                + std::to_string(comp.first) + "] " + std::string(comp.second->GetName());
            BCLog::OutDebug(msg);
            i++;
        }
    }
}
