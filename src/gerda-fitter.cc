/* gerda-fitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: 24/03/2017
 *
 */

#include "GerdaFitter.h"

// STL
#include <iostream>
#include <getopt.h>
#include <chrono>

// BAT
#include "BAT/BCLog.h"

int main(int argc, char** argv) {

    /*
     * get command line args
     */

    std::string progname(argv[0]);

    auto usage = [&]() {
        std::cerr << "USAGE: " << progname << " [-h|--help] json-config\n";
    };

    const char* const short_opts = ":h";
    const option long_opts[] = {
        { "help",  no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0   }
    };

    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                usage();
                return 1;
        }
    }

    // extra arguments
    std::vector<std::string> args;
    for(; optind < argc; optind++){
        args.emplace_back(argv[optind]);
    }

    if (args.empty() or args.size() > 1) {usage(); return 1;}

    std::ifstream fconfig(args[0]);
    if (!fconfig.is_open()) {
        BCLog::OutError("config file " + args[0] + " does not exist");
        return 1;
    }
    json config;
    fconfig >> config;

    /*
     * main routine
     */

    GerdaFitter* model;
    try {
        model = new GerdaFitter(config);
    }
    catch(std::exception& e) {
        BCLog::OutError(e.what());
        BCLog::OutError("caught exception while initializing model, aborting...");
        return 1;
    }

    BCLog::SetLogLevelScreen(config.value("logging", BCLog::summary));

    // set precision (number of samples in Markov chain)
    model->SetPrecision(config.value("precision", BCEngineMCMC::kMedium));

    // run MCMC and marginalize posterior w/r/t all parameters and all
    // combinations of two parameters
    // model->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    auto start = std::chrono::system_clock::now();
    model->MarginalizeAll();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);
    BCLog::OutSummary("Time spent: " + std::to_string(elapsed.count()) + "s");

    // run mode finding, by default using Minuit
    model->FindMode(BCIntegrate::kOptMinuit, model->GetBestFitParameters());
    model->PrintExtendedFitSummary();

    // OUTPUT
    auto outdir = config["output-dir"].get<std::string>();
    std::system(("mkdir -p " + outdir).c_str());
    auto prefix = outdir + "/gerda-fitter-" + config["id"].get<std::string>() + "-";

    // draw parameter plot
    model->PrintParameterPlot(prefix + "parameters.pdf");
    model->PrintParameterLatex(prefix + "parameters.tex");
    model->PrintCorrelationPlot(prefix + "par-correlation.pdf");

    // draw/save all marginalized distributions
    model->WriteMarginalizedDistributions(prefix + "marginalized.root", "recreate");
    model->SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    model->PrintKnowledgeUpdatePlots(prefix + "know-update.pdf");
    model->SaveHistograms(prefix + "histograms.root");

    BCLog::OutSummary("Exiting");
    // close log file
    BCLog::CloseLog();

    delete model;

    return 0;
}
