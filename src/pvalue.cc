/* pvalue.cc
 *
 * simple MCMC to compute p-value
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 31/03/2017
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <random>
#include <chrono>
#include <omp.h>

#include <TFile.h>
#include <TH1D.h>

#include "BAT/BCMath.h"
#include "Fit2nbbLV.h"

double GetPValue(Fit2nbbLV& model, BCEngineMCMC::Precision level, bool save) {

    int nth = 4;
    omp_set_num_threads(nth);

    long int Niter;
    if ( level == BCEngineMCMC::kLow ) Niter = 1E05;
    else                               Niter = 1E06;

    std::cout << "Summary : Calculate p-value with 4 threads (" << Niter << " iterations).\n";
    // get loglikelihood after marginalization
    auto bestpar = model.GetBestFitParameters();
    double logP0 = model.LogLikelihood(bestpar); // ~ -4276
    std::cout << "Summary : logP0 = " << logP0 << std::endl;
    std::cout << "Summary : " << std::flush;

    // get the best fitted function
    auto meanBEGe = model.GetFittedFncBEGe(bestpar);
    auto meanCOAX = model.GetFittedFncCOAX(bestpar);

    //auto dbins = model.GetBinning();
    int nbins = model.GetNbins();
    int downBin = model.GetDownBin();
    int upBin = model.GetUpBin();

    // define a starting hist for the markov chain
    std::vector<int> lCountsBEGe(nbins);
    std::vector<int> lCountsCOAX(nbins);
    for ( int i = 0; i < nbins; ++i ) {
        lCountsBEGe[i] = (int)meanBEGe[i];
        lCountsCOAX[i] = (int)meanCOAX[i];
    }

    // counter for p-value
    long int pv = 0;
    double logP = 0;

    // random numbers generator
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> distr(0,1);

    // write logPs on file (for each thread)
    std::vector<std::unique_ptr<std::ofstream>> f;
    if (save) {
        for ( int i = 0; i < 4; i++ ) {
            std::string str = std::string(std::getenv("GERDACPTDIR")) + "/out/thread" + std::to_string(i) + ".out";
            f.emplace_back( new std::ofstream(str.c_str()));
        }
    }

    auto start = std::chrono::system_clock::now();
    // markov chain, Niter = number of datasets
#pragma omp parallel for reduction(+:pv) firstprivate(meanBEGe,lCountsBEGe,meanCOAX,lCountsCOAX) private(logP)
    for ( int i = 0; i < Niter; ++i ) {
        double r;
        logP = 0;

        // update bin content and loglikelihood
        for ( int j = downBin; j <= upBin; ++j ) {

            // metropolis test BEGe
            if ( distr(eng) < 0.5 ) { // increment
                r = meanBEGe[j]/(lCountsBEGe[j]+1);
                if      ( r >= 1         ) lCountsBEGe[j]++;
                else if ( r > distr(eng) ) lCountsBEGe[j]++;
            }

            else { // decrement
                if ( lCountsBEGe[j] == 1 ) r = 1/meanBEGe[j];
                else r = (lCountsBEGe[j]-1)/meanBEGe[j];
                if      ( r >= 1         ) lCountsBEGe[j]--;
                else if ( r > distr(eng) ) lCountsBEGe[j]--;
            }
            logP += BCMath::LogPoisson(lCountsBEGe[j], meanBEGe[j]);

            // metropolis test COAX
            if ( distr(eng) < 0.5 ) { // increment
                r = meanCOAX[j]/(lCountsCOAX[j]+1);
                if      ( r >= 1         ) lCountsCOAX[j]++;
                else if ( r > distr(eng) ) lCountsCOAX[j]++;
            }

            else { // decrement
                if ( lCountsCOAX[j] == 1 ) r = lCountsCOAX[j]/meanCOAX[j];
                else r = (lCountsCOAX[j]-1)/meanCOAX[j];
                if      ( r >= 1         ) lCountsCOAX[j]--;
                else if ( r > distr(eng) ) lCountsCOAX[j]--;
            }
            logP += BCMath::LogPoisson(lCountsCOAX[j], meanCOAX[j]);
        }

        if (save) *f[omp_get_thread_num()] << logP << std::endl;
        // test
        if ( logP < logP0 ) pv++;
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);
    std::cout << " [" << elapsed.count() << "s]\n";

    return pv*1./Niter;
}
