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
#include <string>

// ROOT
#include "TH1.h"

// BAT
#include "BAT/BCModel.h"

class GerdaFitter : public BCModel {

    public:

    // delete dangerous constructors
    GerdaFitter           ()                   = delete;
    GerdaFitter           (GerdaFitter const&) = delete;
    GerdaFitter& operator=(GerdaFitter const&) = delete;
    // use default destructor
    ~GerdaFitter          ()                   = default;

    // custom constructor
    GerdaFitter(std::string name);

    // methods from BCModel to be overloaded
    double LogLikelihood(const std::vector<double>& parameters);

    std::vector<TH1*> data;
    std::vector<TH1*> pred;
};
