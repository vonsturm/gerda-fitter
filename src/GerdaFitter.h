/* GerdaFitter.h
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#ifndef _GERDA_FITTER_H
#define _GERDA_FITTER_H

// STL
#include <vector>
#include <map>
#include <string>

// ROOT
#include "TH1.h"

// BAT
#include "BAT/BCModel.h"

#include "json.hpp"
using json = nlohmann::json;

class GerdaFitter : public BCModel {

    public:

    // delete dangerous constructors
    GerdaFitter           ()                   = delete;
    GerdaFitter           (GerdaFitter const&) = delete;
    GerdaFitter& operator=(GerdaFitter const&) = delete;

    // custom constructor
    GerdaFitter(json metadata);
    ~GerdaFitter();

    // methods from BCModel to be overloaded
    double LogLikelihood(const std::vector<double>& parameters);

    // map that associates data histogram with list of components
    std::map<TH1*, std::vector<TH1*>> data;
};

#endif
