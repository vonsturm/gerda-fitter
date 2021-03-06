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
#include "TFormula.h"

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
    {BCIntegrate::kIntMonteCarlo, "kIntMonteCarlo"},
    {BCIntegrate::kIntGrid,       "kIntGrid"},
    {BCIntegrate::kIntLaplace,    "kIntLaplace"},
    {BCIntegrate::kIntCuba,       "kIntCuba"},
    {BCIntegrate::kIntDefault,    "kIntDefault"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCCubaMethod, {
    {BCIntegrate::kCubaDivonne, "kCubaDivonne"},
    {BCIntegrate::kCubaVegas,   "kCubaVegas"},
    {BCIntegrate::kCubaSuave,   "kCubaSuave"},
    {BCIntegrate::kCubaCuhre,   "kCubaCuhre"},
    {BCIntegrate::kCubaDefault, "kCubaDefault"},
})

NLOHMANN_JSON_SERIALIZE_ENUM(BCIntegrate::BCOptimizationMethod, {
    {BCIntegrate::kOptEmpty,      "kOptEmpty"},
    {BCIntegrate::kOptSimAnn,     "kOptSimAnn"},
    {BCIntegrate::kOptMetropolis, "kOptMetropolis"},
    {BCIntegrate::kOptMinuit,     "kOptMinuit"},
    {BCIntegrate::kOptDefault,    "kOptDefault"},
})

struct dataset {
    TH1* data;                               // histogram holding data
    std::vector<std::pair<int,int>> brange;  // histogram range (bin index!)
    std::map<int, TH1*> comp;                // catalog of fit components
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
    void CalculateObservables(const std::vector<double> & parameters);

    void SetIntegrationProperties(json config);
    void PrintOptimizationSummary();
    void PrintShortMarginalizationSummary();
    void SaveHistograms(std::string filename);
    void WriteResultsTree(std::string filename);
    double GetFastPValue(const std::vector<double>& parameters, long niter);

    std::vector<dataset> data;
    json config;

    private:

    std::map<std::string,TFormula> obs_tformulas;

    std::string SafeROOTName(const std::string original);
    void DumpData();
    TH1* GetFitComponent(std::string filename, std::string objectname, TH1* data);
};

#endif
