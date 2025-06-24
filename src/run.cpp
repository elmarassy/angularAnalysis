//
// Created by Mero Elmarassy on 6/18/25.
//

#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooFitResult.h>
#include <TH1F.h>
#include <TFile.h>
#include <TStyle.h>

#include "../include/run.h"
#include "../include/individualTimeDependent.h"
#include "../include/parameterValues.h"
#include "../include/plotting.h"

using namespace RooFit;

void runFit(const char* directory, double B0proportion, int nEvents, int nToys, bool masslessApproximation) {

    std::unique_ptr<RooAbsReal> K1s;
    std::unique_ptr<RooAbsReal> K2s;
    std::unique_ptr<RooAbsReal> K2c;
    std::unique_ptr<RooAbsReal> W2s;
    std::unique_ptr<RooAbsReal> W2c;
    std::unique_ptr<RooAbsReal> H2s;
    std::unique_ptr<RooAbsReal> H2c;
    std::unique_ptr<RooAbsReal> Z1s;
    std::unique_ptr<RooAbsReal> Z2s;
    std::unique_ptr<RooAbsReal> Z2c;

    auto makeParam = [=](const char *name, bool floating=true) {
        if (!floating) {
            auto temp = std::make_unique<RooRealVar>(name, name, testValues.at(name));
            temp->setError(0);
            return temp;
        }
        auto temp = std::make_unique<RooRealVar>(name, name, testValues.at(name), -1 -4*abs(testValues.at(name)), 1 + 4*abs(testValues.at(name)));
        temp->setError(0);
        return temp;
    };

    auto x = makeParam("x", false);
    auto y = makeParam("y", false);

    auto K1c = makeParam("K1c");
    auto K3 = makeParam("K3");
    auto K4 = makeParam("K4");
    auto K5 = makeParam("K5");
    auto K6s = makeParam("K6s");
    auto K7 = makeParam("K7");
    auto K8 = makeParam("K8");
    auto K9 = makeParam("K9");

    auto W1s = makeParam("W1s");
    auto W1c = makeParam("W1c");
    auto W3 = makeParam("W3");
    auto W4 = makeParam("W4");
    auto W5 = makeParam("W5");
    auto W6s = makeParam("W6s");
    auto W7 = makeParam("W7");
    auto W8 = makeParam("W8");
    auto W9 = makeParam("W9");

    auto H1s = makeParam("H1s");
    auto H1c = makeParam("H1c");
    auto H3 = makeParam("H3");
    auto H4 = makeParam("H4");
    auto H5 = makeParam("H5");
    auto H6s = makeParam("H6s");
    auto H7 = makeParam("H7");
    auto H8 = makeParam("H8");
    auto H9 = makeParam("H9");

    auto Z1c = makeParam("Z1c");
    auto Z3 = makeParam("Z3");
    auto Z4 = makeParam("Z4");
    auto Z5 = makeParam("Z5");
    auto Z6s = makeParam("Z6s");
    auto Z7 = makeParam("Z7");
    auto Z8 = makeParam("Z8");
    auto Z9 = makeParam("Z9");

    if (masslessApproximation) {
        K1s = std::make_unique<RooFormulaVar>("K1s", "(3.0/4.0) * (1 - y*y - K1c + y*H1c) + y*H1s", RooArgSet(*K1c, *H1s, *H1c, *y));
        K2s = std::make_unique<RooFormulaVar>("K2s", "(1.0/4.0) * (1 - y*y - K1c + y*H1c) + y*H1s / 3.0", RooArgSet(*K1c, *H1s, *H1c, *y));
        K2c = std::make_unique<RooFormulaVar>("K2c", "-K1c", RooArgSet(*K1c));

        W2s = std::make_unique<RooFormulaVar>("W2s", "W1s / 3.0", RooArgSet(*W1s));
        W2c = std::make_unique<RooFormulaVar>("W2c", "-W1c", RooArgSet(*W1c));

        H2s = std::make_unique<RooFormulaVar>("H2s", "H1s/3", RooArgSet(*H1s));
        H2c = std::make_unique<RooFormulaVar>("H2c", "-H1c", RooArgSet(*H1c));

//        if (normalizeBothPdfs) {
//            Z1s = std::make_unique<RooFormulaVar>("Z1s", "(3.0/16.0) * ((16.0*W1s/3.0 + 4*W1c)/x - 4*Z1c)", RooArgSet(*W1s, *W1c, *Z1c, *x));
//            Z2s = std::make_unique<RooFormulaVar>("Z2s", "(1.0/16.0) * ((16.0*W1s/3.0 + 4*W1c)/x - 4*Z1c)", RooArgSet(*W1s, *W1c, *Z1c, *x));
//        } else {
        Z1s = makeParam("Z1s");
        Z2s = std::make_unique<RooFormulaVar>("Z2s", "Z1s/3", RooArgSet(*Z1s));
//        }
        Z2c = std::make_unique<RooFormulaVar>("Z2c", "-Z1c", RooArgSet(*Z1c));
    } else {
        K2c = makeParam("K2c");
        K2s = makeParam("K2s");
        W2s = makeParam("W2s");
        W2c = makeParam("W2c");
        H2s = makeParam("H2s");
        H2c = makeParam("H2c");
        Z1s = makeParam("Z1s");
        Z2s = makeParam("Z2s");
        Z2c = makeParam("Z2c");

        K1s = std::make_unique<RooFormulaVar>("K1s", ""
                                                     "(4*(1-y*y) - 3*K1c + K2c + 2*K2s + y*(3*H1c + 6*H1s - H2c - 2*H2s))/6",
                                              RooArgList(*K2c, *K1c, *K2s, *H1s, *H1c, *H2s, *H2c, *y));
    }

    RooRealVar cosThetaK("cosThetaK", "cosThetaK", -1, 1);
    RooRealVar cosThetaL("cosThetaL", "cosThetaL", -1, 1);
    RooRealVar phi("phi", "phi", -TMath::Pi(), TMath::Pi());
    RooRealVar t("t", "t", 0, 10);

    timeDependent b("b", "b", cosThetaL, cosThetaK, phi, t, 1, *x, *y,
                    *K1s, *K1c, *K2s, *K2c, *K3, *K4, *K5, *K6s, *K7, *K8, *K9,
                    *W1s, *W1c, *W2s, *W2c, *W3, *W4, *W5, *W6s, *W7, *W8, *W9,
                    *H1s, *H1c, *H2s, *H2c, *H3, *H4, *H5, *H6s, *H7, *H8, *H9,
                    *Z1s, *Z1c, *Z2s, *Z2c, *Z3, *Z4, *Z5, *Z6s, *Z7, *Z8, *Z9);

    timeDependent bBar("bBar", "bBar", cosThetaL, cosThetaK, phi, t, -1, *x, *y,
                       *K1s, *K1c, *K2s, *K2c, *K3, *K4, *K5, *K6s, *K7, *K8, *K9,
                       *W1s, *W1c, *W2s, *W2c, *W3, *W4, *W5, *W6s, *W7, *W8, *W9,
                       *H1s, *H1c, *H2s, *H2c, *H3, *H4, *H5, *H6s, *H7, *H8, *H9,
                       *Z1s, *Z1c, *Z2s, *Z2c, *Z3, *Z4, *Z5, *Z6s, *Z7, *Z8, *Z9);

    RooCategory type("type", "type");
    type.defineType("B0");
    type.defineType("B0Bar");

    auto observables = RooArgSet(cosThetaL, cosThetaK, phi, t);

    RooSimultaneous simPdf("simPdf", "Simultaneous PDF", {{"B0", &b}, {"B0Bar", &bBar}}, type);

//    RooFitResult* result = simPdf.fitTo(combData, Save(), Strategy(2), Hesse(true), Save(true), PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(12), Precision(1e-5));
//
//    auto bFitData = static_cast<RooDataSet*>(bData.get());
//    auto bBarFitData = static_cast<RooDataSet*>(bBarData.get());
//
//    plotInitial(cosThetaL, bFitData, *b.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/plots/B0/cosThetaL.png");
//    plotInitial(cosThetaK, bFitData, *b.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/plots/B0/cosThetaK.png");
//    plotInitial(phi, bFitData, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/plots/B0/phi.png");
//    plotInitial(t, bFitData, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/plots/B0/t.png");
//    plotInitial(cosThetaL, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/plots/B0bar/cosThetaL.png");
//    plotInitial(cosThetaK, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/plots/B0bar/cosThetaK.png");
//    plotInitial(phi, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/plots/B0bar/phi.png");
//    plotInitial(t, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/plots/B0bar/t.png");
//    plot(cosThetaL, bFitData, result, *b.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/fits/B0/cosThetaL.png");
//    plot(cosThetaK, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/fits/B0/cosThetaK.png");
//    plot(phi, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/fits/B0/phi.png");
//    plot(t, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/fits/B0/t.png");
//    plot(cosThetaL, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/fits/B0bar/cosThetaL.png");
//    plot(cosThetaK, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/fits/B0bar/cosThetaK.png");
//    plot(phi, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/fits/B0bar/phi.png");
//    plot(t, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/fits/B0bar/t.png");
    const char* title = Form("%s/run_%d_%.3f_%d_%d.root", directory, masslessApproximation, B0proportion, nEvents, nToys);
    std::cout << "Started running; saving to file " << title << std::endl;
    runPulls(title, &simPdf, type, observables, testValues, false, B0proportion, nToys, nEvents);
//    plotPulls("testing.root", pulls, masslessApproximation);
}


void runPulls(const char* saveFile, RooSimultaneous* pdf, RooCategory type, const RooArgSet& observables,
              const std::map<std::string, double>& trueValues, bool useAsymmetricError, double B0proportion=0.5, int nToys = 100, int nEvents = 500) {

    auto outFile = TFile::Open(saveFile, "RECREATE");
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.09);
    gStyle->SetLabelSize(0.07, "XY");

    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Could not open file '" << saveFile << "'" << std::endl;
        return;
    }
    TDirectory* projectionDir = outFile->mkdir("projections");
    std::map<std::string, std::pair<std::vector<double>, std::vector<double>>> results;

    for (auto& p: trueValues) {
        std::pair<std::vector<double>, std::vector<double>> result;
        results[p.first] = result;
    }

    auto fullObservables = dynamic_cast<RooArgSet*>(observables.Clone());
    fullObservables->add(type);
//    double proportion = pdf->getPdf("B0")->createIntegral(observables)->getVal();
    std::cout << "Generating data with B0 proportion = " << B0proportion << std::endl;
    int bEvents = int(nEvents * B0proportion);
    for (int i = 0; i < nToys; i++) {

        std::unique_ptr<RooDataSet> bData{pdf->getPdf("B0")->generate(observables, bEvents)};
        std::unique_ptr<RooDataSet> bBarData{pdf->getPdf("B0Bar")->generate(observables, nEvents - bEvents)};

        RooDataSet toy("combData", "combined data", observables, Index(type),
                       Import({{"B0", bData.get()}, {"B0Bar", bBarData.get()}}));

        std::unique_ptr<RooFitResult> fitResult(pdf->fitTo(toy, Save(true), Minos(useAsymmetricError), Hesse(true), Strategy(2), PrintLevel(-1), EvalBackend("legacy"), NumCPU(10)));

        if (fitResult->status() != 0) {
            i --;
            std::cout << "Invalid minimum (status " << fitResult->status() << "), retrying..." << std::endl;
            continue;
        }
        std::map<std::string, std::pair<double, double>> currentPullResults;
        const double saveValue = 3;
        bool saveProjections = false;

        for (auto& p : trueValues) {
            auto* par = dynamic_cast<RooRealVar*>(fitResult->floatParsFinal().find(p.first.c_str()));
            double val;
            if (useAsymmetricError) {
                if (par && par->getAsymErrorLo() < 0 && par->getAsymErrorHi() > 0) {
                    val = par->getVal();
                    results[p.first].first.push_back(val);
                    currentPullResults[p.first].first = val;
                    const auto error = (par->getVal() - p.second > 0)? par->getAsymErrorLo(): par->getAsymErrorHi();
                    results[p.first].second.push_back(error);
                    currentPullResults[p.first].second = error;
                    if (std::abs((val - p.second)) / error > saveValue) saveProjections = true;
                }
            } else {
                if (par && par->getError() > 0) {
                    val = par->getVal();
                    results[p.first].first.push_back(val);
                    currentPullResults[p.first].first = val;
                    const auto error = par->getError();
                    results[p.first].second.push_back(error);
                    currentPullResults[p.first].second = error;
                    if (std::abs((val - p.second)) / error > saveValue) saveProjections = true;
                }
            }
        }
        if (saveProjections) {
            std::cout << "Large pull value in toy #" << i << ", saving projections" << std::endl;
            auto dir = projectionDir->mkdir(Form("Toy #%d", i));
            makePullProjections(currentPullResults, &*fitResult, *pdf, &toy, dir, i);
        }
        resetPdf(pdf, *fullObservables);
        std::cout << "Completed toy #" << i << std::endl;
    }

    plotPulls(results, outFile);
    outFile->Close();
    std::cout << "Closed file" << std::endl;
}


void resetPdf(RooSimultaneous* pdf, const RooArgSet& observables) {
    for (auto p: *pdf->getParameters(observables)) {
        auto param = dynamic_cast<RooRealVar *>(p);
        param->setVal(testValues.at(param->GetName()));
    }
}