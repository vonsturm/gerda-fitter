/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.h"

// BAT
#include "BAT/BCLog.h"

// ROOT
#include "TFile.h"
#include "TH1.h"

GerdaFitter::GerdaFitter(json config) {
    this->SetName(config["id"]);

    // loop over files with data histograms
    for (auto& el : config["fit"].items()) {
        TFile _tf(el.key().c_str());
        if (!_tf.IsOpen()) throw std::runtime_error("invalid ROOT file: " + el.key());

        // loop over requested histograms in file
        for (auto& elh : el.value().items()) {
            auto th = dynamic_cast<TH1*>(_tf.Get(elh.key().c_str()));
            if (!th) throw std::runtime_error("could not find object '" + elh.key() + "' in file " + el.key());
            th->SetDirectory(0);
            data.insert({th, {}});
            BCLog::OutDebug("Adding data histogram '" + elh.key() + "' from " + el.key());

            // loop over requested components
            for (auto& it : elh.value()) {
                auto prefix = it.value("prefix", ".");

                // loop over requested isotopes on the relative part
                for (auto& iso : it["isotopes"]) {
                    if (iso.is_string()) {
                        TFile _tff((prefix + "").c_str());
                    }
                    else if (iso.is_object()) {
                    }
                    else throw std::runtime_error("unexpected entry " + iso.dump() + "found in [\"fit\"][\"" + el.key() + "\"][\"" + elh.key() + "\"][\"isotopes\"]");
                }
            }
        }
    }
}

GerdaFitter::~GerdaFitter() {
    for (auto& h : data) {
        delete h.first;
        for (auto& hh : h.second) delete hh;
    }
}

double GerdaFitter::LogLikelihood(const std::vector<double>& parameters) {
    return 0;
}
