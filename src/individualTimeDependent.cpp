//
// Created by Mero Elmarassy on 6/10/25.
//
#include <TF1.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooMCStudy.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooFFTConvPdf.h>
#include <TFile.h>

#include "../include/individualTimeDependent.h"
//#include "../include/parameterValues.h"

using namespace RooFit;

timeDependent::timeDependent(const char *name, const char *title,
                             RooAbsReal& cosThetaL, RooAbsReal& cosThetaK, RooAbsReal& phi, RooAbsReal& t,
                             int sign, RooAbsReal& x, RooAbsReal& y,
                             RooAbsReal& K1s, RooAbsReal& K1c, RooAbsReal& K2s, RooAbsReal& K2c, RooAbsReal& K3,
                             RooAbsReal& K4, RooAbsReal& K5, RooAbsReal& K6s, RooAbsReal& K7, RooAbsReal& K8, RooAbsReal& K9,
                             RooAbsReal& W1s, RooAbsReal& W1c, RooAbsReal& W2s, RooAbsReal& W2c, RooAbsReal& W3,
                             RooAbsReal& W4, RooAbsReal& W5, RooAbsReal& W6s, RooAbsReal& W7, RooAbsReal& W8, RooAbsReal& W9,
                             RooAbsReal& H1s, RooAbsReal& H1c, RooAbsReal& H2s, RooAbsReal& H2c, RooAbsReal& H3,
                             RooAbsReal& H4, RooAbsReal& H5, RooAbsReal& H6s, RooAbsReal& H7, RooAbsReal& H8, RooAbsReal& H9,
                             RooAbsReal& Z1s, RooAbsReal& Z1c, RooAbsReal& Z2s, RooAbsReal& Z2c, RooAbsReal& Z3,
                             RooAbsReal& Z4, RooAbsReal& Z5, RooAbsReal& Z6s, RooAbsReal& Z7, RooAbsReal& Z8, RooAbsReal& Z9):
        RooAbsPdf(name, title),
        sign_(sign),
        cosThetaL_("costhetal", "costhetal", this, cosThetaL),
        cosThetaK_("costhetak", "costhetak", this, cosThetaK),
        phi_("phi", "phi", this, phi),
        t_("t", "t", this, t),
        x_("x", "x", this, x),
        y_("y", "y", this, y),
        K1s_("K1s", "K1s", this, K1s),
        K1c_("K1c", "K1c", this, K1c),
        K2s_("K2s", "K2s", this, K2s),
        K2c_("K2c", "K2c", this, K2c),
        K3_("K3", "K3", this, K3),
        K4_("K4", "K4", this, K4),
        K5_("K5", "K5", this, K5),
        K6s_("K6s", "K6s", this, K6s),
        K7_("K7", "K7", this, K7),
        K8_("K8", "K8", this, K8),
        K9_("K9", "K9", this, K9),
        W1s_("W1s", "W1s", this, W1s),
        W1c_("W1c", "W1c", this, W1c),
        W2s_("W2s", "W2s", this, W2s),
        W2c_("W2c", "W2c", this, W2c),
        W3_("W3", "W3", this, W3),
        W4_("W4", "W4", this, W4),
        W5_("W5", "W5", this, W5),
        W6s_("W6s", "W6s", this, W6s),
        W7_("W7", "W7", this, W7),
        W8_("W8", "W8", this, W8),
        W9_("W9", "W9", this, W9),
        H1s_("H1s", "H1s", this, H1s),
        H1c_("H1c", "H1c", this, H1c),
        H2s_("H2s", "H2s", this, H2s),
        H2c_("H2c", "H2c", this, H2c),
        H3_("H3", "H3", this, H3),
        H4_("H4", "H4", this, H4),
        H5_("H5", "H5", this, H5),
        H6s_("H6s", "H6s", this, H6s),
        H7_("H7", "H7", this, H7),
        H8_("H8", "H8", this, H8),
        H9_("H9", "H9", this, H9),
        Z1s_("Z1s", "Z1s", this, Z1s),
        Z1c_("Z1c", "Z1c", this, Z1c),
        Z2s_("Z2s", "Z2s", this, Z2s),
        Z2c_("Z2c", "Z2c", this, Z2c),
        Z3_("Z3", "Z3", this, Z3),
        Z4_("Z4", "Z4", this, Z4),
        Z5_("Z5", "Z5", this, Z5),
        Z6s_("Z6s", "Z6s", this, Z6s),
        Z7_("Z7", "Z7", this, Z7),
        Z8_("Z8", "Z8", this, Z8),
        Z9_("Z9", "Z9", this, Z9) {}

timeDependent::timeDependent(timeDependent const &other, const char *name):
        RooAbsPdf(other, name),
        sign_(other.sign_),
        cosThetaL_("costhetal", this, other.cosThetaL_),
        cosThetaK_("costhetak", this, other.cosThetaK_),
        phi_("phi", this, other.phi_),
        t_("t", this, other.t_),
        x_("x", this, other.x_),
        y_("y", this, other.y_),
        K1s_("K1s", this, other.K1s_),
        K1c_("K1c", this, other.K1c_),
        K2s_("K2s", this, other.K2s_),
        K2c_("K2c", this, other.K2c_),
        K3_("K3", this, other.K3_),
        K4_("K4", this, other.K4_),
        K5_("K5", this, other.K5_),
        K6s_("K6s", this, other.K6s_),
        K7_("K7", this, other.K7_),
        K8_("K8", this, other.K8_),
        K9_("K9", this, other.K9_),
        W1s_("W1s", this, other.W1s_),
        W1c_("W1c", this, other.W1c_),
        W2s_("W2s", this, other.W2s_),
        W2c_("W2c", this, other.W2c_),
        W3_("W3", this, other.W3_),
        W4_("W4", this, other.W4_),
        W5_("W5", this, other.W5_),
        W6s_("W6s", this, other.W6s_),
        W7_("W7", this, other.W7_),
        W8_("W8", this, other.W8_),
        W9_("W9", this, other.W9_),
        H1s_("H1s", this, other.H1s_),
        H1c_("H1c", this, other.H1c_),
        H2s_("H2s", this, other.H2s_),
        H2c_("H2c", this, other.H2c_),
        H3_("H3", this, other.H3_),
        H4_("H4", this, other.H4_),
        H5_("H5", this, other.H5_),
        H6s_("H6s", this, other.H6s_),
        H7_("H7", this, other.H7_),
        H8_("H8", this, other.H8_),
        H9_("H9", this, other.H9_),
        Z1s_("Z1s", this, other.Z1s_),
        Z1c_("Z1c", this, other.Z1c_),
        Z2s_("Z2s", this, other.Z2s_),
        Z2c_("Z2c", this, other.Z2c_),
        Z3_("Z3", this, other.Z3_),
        Z4_("Z4", this, other.Z4_),
        Z5_("Z5", this, other.Z5_),
        Z6s_("Z6s", this, other.Z6s_),
        Z7_("Z7", this, other.Z7_),
        Z8_("Z8", this, other.Z8_),
        Z9_("Z9", this, other.Z9_) {}


void plot(RooRealVar& variable, RooDataSet* fitData, RooFitResult* result, RooAbsPdf& pdf, char* location) {
    RooPlot* fitFrame = variable.frame(25);
    auto* fitResult = const_cast<RooFitResult *>(result);
    pdf.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow), Normalization(fitData->sumEntries(), RooAbsReal::NumEvent));
    pdf.plotOn(fitFrame, Normalization(fitData->sumEntries(), RooAbsReal::NumEvent));
    fitData->plotOn(fitFrame);
    auto fitCanvas = new TCanvas("Fit", "Variable", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs(location);
    delete fitCanvas;
}

void plotInitial(RooRealVar& variable, RooDataSet* fitData, RooAbsPdf& pdf, char* location) {
    RooPlot* fitFrame = variable.frame(25);
    pdf.plotOn(fitFrame, Normalization(fitData->sumEntries(), RooAbsReal::NumEvent));
    fitData->plotOn(fitFrame);
    auto fitCanvas = new TCanvas("Fit", "Variable", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs(location);
    delete fitCanvas;
}


//
//
//void runFit(int nEvents, int nToys, bool masslessApproximation, bool normalizeBothPdfs) {
//
//    std::unique_ptr<RooAbsReal> K1s;
//    std::unique_ptr<RooAbsReal> K2s;
//    std::unique_ptr<RooAbsReal> K2c;
//    std::unique_ptr<RooAbsReal> W2s;
//    std::unique_ptr<RooAbsReal> W2c;
//    std::unique_ptr<RooAbsReal> H2s;
//    std::unique_ptr<RooAbsReal> H2c;
//    std::unique_ptr<RooAbsReal> Z1s;
//    std::unique_ptr<RooAbsReal> Z2s;
//    std::unique_ptr<RooAbsReal> Z2c;
//
//    auto makeParam = [=](const char *name, bool floating=true) {
//        if (!floating) {
//            auto temp = std::make_unique<RooRealVar>(name, name, testValues.at(name));
//            temp->setError(0);
//            return temp;
//        }
//        auto temp = std::make_unique<RooRealVar>(name, name, testValues.at(name), -1 -4*abs(testValues.at(name)), 1 + 4*abs(testValues.at(name)));
//        temp->setError(0);
//        return temp;
//    };
//
//    auto x = makeParam("x", false);
//    auto y = makeParam("y", false);
//
//    auto K1c = makeParam("K1c");
//    auto K3 = makeParam("K3");
//    auto K4 = makeParam("K4");
//    auto K5 = makeParam("K5");
//    auto K6s = makeParam("K6s");
//    auto K7 = makeParam("K7");
//    auto K8 = makeParam("K8");
//    auto K9 = makeParam("K9");
//
//    auto W1s = makeParam("W1s");
//    auto W1c = makeParam("W1c");
//    auto W3 = makeParam("W3");
//    auto W4 = makeParam("W4");
//    auto W5 = makeParam("W5");
//    auto W6s = makeParam("W6s");
//    auto W7 = makeParam("W7");
//    auto W8 = makeParam("W8");
//    auto W9 = makeParam("W9");
//
//    auto H1s = makeParam("H1s");
//    auto H1c = makeParam("H1c");
//    auto H3 = makeParam("H3");
//    auto H4 = makeParam("H4");
//    auto H5 = makeParam("H5");
//    auto H6s = makeParam("H6s");
//    auto H7 = makeParam("H7");
//    auto H8 = makeParam("H8");
//    auto H9 = makeParam("H9");
//
////    auto Z1s = makeParam("Z1s");
//    auto Z1c = makeParam("Z1c");
//    auto Z3 = makeParam("Z3");
//    auto Z4 = makeParam("Z4");
//    auto Z5 = makeParam("Z5");
//    auto Z6s = makeParam("Z6s");
//    auto Z7 = makeParam("Z7");
//    auto Z8 = makeParam("Z8");
//    auto Z9 = makeParam("Z9");
//
//    if (masslessApproximation) {
//        K1s = std::make_unique<RooFormulaVar>("K1s", "(3.0/4.0) * (1 - y*y - K1c + y*H1c) + y*H1s", RooArgSet(*K1c, *H1s, *H1c, *y));
//        K2s = std::make_unique<RooFormulaVar>("K2s", "(1.0/4.0) * (1 - y*y - K1c + y*H1c) + y*H1s / 3.0", RooArgSet(*K1c, *H1s, *H1c, *y));
//        K2c = std::make_unique<RooFormulaVar>("K2c", "-K1c", RooArgSet(*K1c));
//
//        W2s = std::make_unique<RooFormulaVar>("W2s", "W1s / 3.0", RooArgSet(*W1s));
//        W2c = std::make_unique<RooFormulaVar>("W2c", "-W1c", RooArgSet(*W1c));
//
//        H2s = std::make_unique<RooFormulaVar>("H2s", "H1s/3", RooArgSet(*H1s));
//        H2c = std::make_unique<RooFormulaVar>("H2c", "-H1c", RooArgSet(*H1c));
//
//        if (normalizeBothPdfs) {
//            Z1s = std::make_unique<RooFormulaVar>("Z1s", "(3.0/16.0) * ((16.0*W1s/3.0 + 4*W1c)/x - 4*Z1c)", RooArgSet(*W1s, *W1c, *Z1c, *x));
//            Z2s = std::make_unique<RooFormulaVar>("Z2s", "(1.0/16.0) * ((16.0*W1s/3.0 + 4*W1c)/x - 4*Z1c)", RooArgSet(*W1s, *W1c, *Z1c, *x));
//        } else {
//            Z1s = makeParam("Z1s");
//            Z2s = std::make_unique<RooFormulaVar>("Z2s", "Z1s/3", RooArgSet(*Z1s));
//        }
//        Z2c = std::make_unique<RooFormulaVar>("Z2c", "-Z1c", RooArgSet(*Z1c));
//    } else {
//        K2c = makeParam("K2c");
//        K2s = makeParam("K2s");
//        W2s = makeParam("W2s");
//        W2c = makeParam("W2c");
//        H2s = makeParam("H2s");
//        H2c = makeParam("H2c");
//        Z1s = makeParam("Z1s");
//        Z2s = makeParam("Z2s");
//        Z2c = makeParam("Z2c");
//
//        K1s = std::make_unique<RooFormulaVar>("K1s", ""
//                                                     "(4*(1-y*y) - 3*K1c + K2c + 2*K2s + y*(3*H1c + 6*H1s - H2c - 2*H2s))/6",
//                                              RooArgList(*K2c, *K1c, *K2s, *H1s, *H1c, *H2s, *H2c, *y));
//    }
//
//    RooRealVar cosThetaK("cosThetaK", "cosThetaK", -1, 1);
//    RooRealVar cosThetaL("cosThetaL", "cosThetaL", -1, 1);
//    RooRealVar phi("phi", "phi", -TMath::Pi(), TMath::Pi());
//    RooRealVar t("t", "t", 0, 10);
//
//    timeDependent b("b", "b", cosThetaL, cosThetaK, phi, t, 1, *x, *y,
//                    *K1s, *K1c, *K2s, *K2c, *K3, *K4, *K5, *K6s, *K7, *K8, *K9,
//                    *W1s, *W1c, *W2s, *W2c, *W3, *W4, *W5, *W6s, *W7, *W8, *W9,
//                    *H1s, *H1c, *H2s, *H2c, *H3, *H4, *H5, *H6s, *H7, *H8, *H9,
//                    *Z1s, *Z1c, *Z2s, *Z2c, *Z3, *Z4, *Z5, *Z6s, *Z7, *Z8, *Z9);
//
//    timeDependent bBar("bBar", "bBar", cosThetaL, cosThetaK, phi, t, -1, *x, *y,
//                       *K1s, *K1c, *K2s, *K2c, *K3, *K4, *K5, *K6s, *K7, *K8, *K9,
//                       *W1s, *W1c, *W2s, *W2c, *W3, *W4, *W5, *W6s, *W7, *W8, *W9,
//                       *H1s, *H1c, *H2s, *H2c, *H3, *H4, *H5, *H6s, *H7, *H8, *H9,
//                       *Z1s, *Z1c, *Z2s, *Z2c, *Z3, *Z4, *Z5, *Z6s, *Z7, *Z8, *Z9);
//
//    RooCategory type("type", "type");
//    type.defineType("B0");
//    type.defineType("B0Bar");
//
//    auto observables = RooArgSet(cosThetaL, cosThetaK, phi, t);
//
////    std::unique_ptr<RooDataSet> bData{b.generate(observables, nEvents)};
////
////    std::unique_ptr<RooDataSet> bBarData{bBar.generate(observables, nEvents)};
////
////    RooDataSet combData("combData", "combined data", RooArgSet(cosThetaL, cosThetaK, phi, t), Index(type),
////                        Import({{"B0", bData.get()}, {"B0Bar", bBarData.get()}}));
//
//    RooSimultaneous simPdf("simPdf", "Simultaneous PDF", {{"B0", &b}, {"B0Bar", &bBar}}, type);
//
////    RooFitResult* result = simPdf.fitTo(combData, Save(), Strategy(2), Hesse(true), Save(true), PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(12), Precision(1e-5));
////
////    auto bFitData = static_cast<RooDataSet*>(bData.get());
////    auto bBarFitData = static_cast<RooDataSet*>(bBarData.get());
////
////    plotInitial(cosThetaL, bFitData, *b.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/plots/B0/cosThetaL.png");
////    plotInitial(cosThetaK, bFitData, *b.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/plots/B0/cosThetaK.png");
////    plotInitial(phi, bFitData, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/plots/B0/phi.png");
////    plotInitial(t, bFitData, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/plots/B0/t.png");
////    plotInitial(cosThetaL, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/plots/B0bar/cosThetaL.png");
////    plotInitial(cosThetaK, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/plots/B0bar/cosThetaK.png");
////    plotInitial(phi, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/plots/B0bar/phi.png");
////    plotInitial(t, bBarFitData, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/plots/B0bar/t.png");
////    plot(cosThetaL, bFitData, result, *b.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/fits/B0/cosThetaL.png");
////    plot(cosThetaK, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/fits/B0/cosThetaK.png");
////    plot(phi, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/fits/B0/phi.png");
////    plot(t, bFitData, result, *b.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/fits/B0/t.png");
////    plot(cosThetaL, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaK, phi, t)), "timeDependent/fits/B0bar/cosThetaL.png");
////    plot(cosThetaK, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, phi, t)), "timeDependent/fits/B0bar/cosThetaK.png");
////    plot(phi, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, t)), "timeDependent/fits/B0bar/phi.png");
////    plot(t, bBarFitData, result, *bBar.createProjection(RooArgSet(cosThetaL, cosThetaK, phi)), "timeDependent/fits/B0bar/t.png");
//
//    std::map<std::string, TH1F*> pulls = RunPullStudy(&simPdf, type, observables, testValues, false, nToys, nEvents);
//    plotPulls("testing.root", pulls, masslessApproximation);
//}