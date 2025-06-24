//
// Created by Mero Elmarassy on 6/18/25.
//

#ifndef CERN_RUN_H
#define CERN_RUN_H

void resetPdf(RooSimultaneous* pdf, const RooArgSet& observables);

void runFit(const char* directory, double B0proportion, int nEvents, int nToys, bool masslessApproximation);

void runPulls(const char* saveFile, RooSimultaneous* pdf, RooCategory type, const RooArgSet& observables,
              const std::map<std::string, double>& trueValues, bool useAsymmetricError, double B0proportion, int nToys, int nEvents);

#endif //CERN_RUN_H
