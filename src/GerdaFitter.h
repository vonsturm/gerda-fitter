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

NLOHMANN_JSON_SERIALIZE_ENUM(BCEngineMCMC::Precision, {
    {BCEngineMCMC::kQuick,    "kQuick"},
    {BCEngineMCMC::kLow,      "kLow"},
    {BCEngineMCMC::kMedium,   "kMedium"},
    {BCEngineMCMC::kHigh,     "kHigh"},
    {BCEngineMCMC::kVeryHigh, "kVeryHigh"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCLog::LogLevel, {
    {BCLog::debug,   "debug"},
    {BCLog::detail,  "detail"},
    {BCLog::summary, "summary"},
    {BCLog::warning, "warning"},
    {BCLog::error,   "error"},
    {BCLog::nothing, "nothing"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCIntegrationMethod, {
    {BCEngineMCMC::kIntMonteCarlo, "kIntMonteCarlo"},
    {BCEngineMCMC::kIntGrid,       "kIntGrid"},
    {BCEngineMCMC::kIntLaplace,    "kIntLaplace"},
    {BCEngineMCMC::kIntCuba,       "kIntCuba"},
    {BCEngineMCMC::kIntDefault,    "kIntDefault"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCCubaMethod, {
    {BCEngineMCMC::kCubaDivonne, "kCubaDivonne"},
    {BCEngineMCMC::kCubaVegas,   "kCubaVegas"},
    {BCEngineMCMC::kCubaSuave,   "kCubaSuave"},
    {BCEngineMCMC::kCubaCuhre,   "kCubaCuhre"},
    {BCEngineMCMC::kCubaDefault, "kCubaDefault"},
})

struct dataset {
    TH1* data;
    std::pair<int,int> brange;
    std::map<int, TH1*> comp;
};

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

    void SetIntegrationProperties(json config);
    void PrintExtendedFitSummary();
    void SaveHistograms(std::string filename);

    std::vector<dataset> data;
    json config;

    private:

    void DumpData();
    TH1* GetFitComponent(std::string filename, std::string objectname, TH1* data);
};

#endif
