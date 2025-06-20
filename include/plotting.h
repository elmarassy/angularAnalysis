//
// Created by Mero Elmarassy on 6/18/25.
//

#ifndef CERN_PLOTTING_H
#define CERN_PLOTTING_H

#include <string>
#include <map>

void plotPulls(const std::map<std::string, std::pair<std::vector<double>, std::vector<double>>>& results, TFile* outFile);


void makePullProjections(std::map<std::string, std::pair<double, double>> currentPullResults, RooFitResult* fitResult, RooSimultaneous& pdf, RooDataSet* fitData, TDirectory* pullDir, int toyNumber);




#endif //CERN_PLOTTING_H
