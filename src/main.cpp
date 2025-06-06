#include <iostream>
#include <cmath>
#include <chrono>
#include <algorithm>

#include "TRandom3.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TError.h"

#include "../include/fit.h"

using namespace RooFit;


void runTimeIndependent() {
    RooRealVar costhetal("costhetal", "costhetal", -1.0, 1.0);
    RooRealVar costhetak("costhetak", "costhetak", -1.0, 1.0);
    RooRealVar phi("phi", "phi", -TMath::Pi(), +TMath::Pi());

    RooRealVar J1s("J1s", "J1s", 0.3, 0, 10);
    RooRealVar J1c("J1c", "J1c", 0.8, -10, 10);
    RooRealVar J2s("J2s", "J2s", 0, -10, 10);
    RooFormulaVar J2c("J2c", "3*J1c + 6*J1s - 2*J2s - 4", RooArgList(J1c, J1s, J2s));

    RooRealVar J3("J3", "J3", 0, -10, 10);
    RooRealVar J4("J4", "J4", 0, -10, 10);
    RooRealVar J5("J5", "J5", 0, -10, 10);
    RooRealVar J6("J6", "J6", 0, -10, 10);
    RooRealVar J7("J7", "J7", 0, -10, 10);
    RooRealVar J8("J8", "J8", 0, -10, 10);
    RooRealVar J9("J9", "J9", 0, -10, 10);

    J1s.setError(0.01);
    J1c.setError(0.01);
    J2s.setError(0.01);
    J3.setError(0.01);
    J4.setError(0.01);
    J5.setError(0.01);
    J6.setError(0.01);
    J7.setError(0.01);
    J8.setError(0.01);
    J9.setError(0.01);

    timeIndependent model("model", "model", costhetal, costhetak, phi, J1s, J1c, J2s, J2c, J3, J4, J5, J6, J7, J8, J9);


    J1s.setVal(0.3);
    J1c.setVal(0.8);
    J2s.setVal(0.0);
    J3.setVal(0.0);
    J4.setVal(0.0);
    J5.setVal(0.0);
    J6.setVal(0.0);
    J7.setVal(0.0);
    J8.setVal(0.0);
    J9.setVal(0.0);

    J1s.setError(0.01);
    J1c.setError(0.01);
    J2s.setError(0.01);
    J3.setError(0.01);
    J4.setError(0.01);
    J5.setError(0.01);
    J6.setError(0.01);
    J7.setError(0.01);
    J8.setError(0.01);
    J9.setError(0.01);

    auto t_before_toystudy = std::chrono::high_resolution_clock::now();

    unsigned int ngen = 10000;
    unsigned int nruns = 10;
    RooMCStudy *mcstudy =
            new RooMCStudy(model, RooArgSet(costhetal, costhetak, phi),Binned(false), FitOptions(Strategy(2),Hesse(true), Save(true),PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(10)));

    mcstudy->generateAndFit(nruns, ngen, true);

    auto t_after_toystudy = std::chrono::high_resolution_clock::now();

    const RooArgSet* fitParams = mcstudy->fitParams(0);

    TIterator* iter = fitParams->createIterator();
    RooRealVar* var = nullptr;

    while ((var = (RooRealVar*)iter->Next())) {
        const char* name = var->GetName();
        try {
            RooPlot* pullFrame = mcstudy->plotPull(*var, -4, 4);
            if (!pullFrame || pullFrame->numItems() == 0) {
                continue;
            }
            TCanvas* pullCanvas = new TCanvas(Form("c_%s", name), Form("Pull of %s", name), 1600, 1200);
            pullFrame->Draw();
            pullCanvas->SaveAs(Form("pulls2/%s.png", name));
            delete pullCanvas;
        } catch (const std::exception& e) {
            continue;
        }
    }
    delete iter;


    RooPlot* fitFrame = costhetal.frame(100);
    RooDataSet* toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    RooFitResult* fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    TCanvas* fitCanvas = new TCanvas("Fit", "CosThetaL", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/cosThetaL.png");
    delete fitCanvas;

    fitFrame = costhetak.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "CosThetaK", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/cosThetaK.png");
    delete fitCanvas;

    fitFrame = phi.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "Phi", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/phi.png");
    delete fitCanvas;

    std::cout << "Processed " << nruns << " toys with " << ngen << " events in " << std::chrono::duration<double, std::milli>(t_after_toystudy-t_before_toystudy).count()/1000 << "s" << std::endl;

}

void runTimeDependentUntagged() {
    RooRealVar costhetal("costhetal", "costhetal", -1.0, 1.0);
    RooRealVar costhetak("costhetak", "costhetak", -1.0, 1.0);
    RooRealVar phi("phi", "phi", -TMath::Pi(), +TMath::Pi());
    RooRealVar t("t", "t", 0, 10);

    RooRealVar gamma("gamma", "gamma", 1);
    RooRealVar y("y", "y", 0.124);

    RooRealVar SJ1s("SJ1s", "SJ1s", 0.3, 0, 10);
    RooRealVar SJ2s("SJ2s", "SJ2s", 0, -10, 10);
    RooRealVar SJ2c("SJ2c", "SJ2c", 0, -10, 10);
    RooRealVar SJ3("SJ3", "SJ3", 0, -10, 10);
    RooRealVar SJ4("SJ4", "SJ4", 0, -10, 10);
    RooRealVar DJ5("DJ5", "DJ5", 0, -10, 10);
    RooRealVar DJ6s("DJ6s", "DJ6s", 0, -10, 10);
    RooRealVar SJ7("SJ7", "SJ7", 0, -10, 10);
    RooRealVar DJ8("DJ8", "DJ8", 0, -10, 10);
    RooRealVar DJ9("DJ9", "DJ9", 0, -10, 10);
    RooRealVar h1s("h1s", "h1s", 0, -10, 10);
    RooRealVar h1c("h1c", "h1c", 0, -10, 10);
    RooRealVar h2s("h2s", "h2s", 0, -10, 10);
    RooRealVar h2c("h2c", "h2c", 0, -10, 10);
    RooRealVar h3("h3", "h3", 0, -10, 10);
    RooRealVar h4("h4", "h4", 0, -10, 10);
    RooRealVar h5("h5", "h5", 0, -10, 10);
    RooRealVar h6s("h6s", "h6s", 0, -10, 10);
    RooRealVar h7("h7", "h7", 0, -10, 10);
    RooRealVar h8("h8", "h8", 0, -10, 10);
    RooRealVar h9("h9", "h9", 0, -10, 10);

    RooFormulaVar SJ1c("SJ1c", "(4*(1 - y**2) - (2*h2s*y) - (h2c*y) + (6*h1s*y) + (3*h1c*y) + (2*SJ2s) + (SJ2c) - (6*SJ1s)) / 3", RooArgList(y, h1s, h1c, h2s, h2c, SJ1s, SJ2s, SJ2c));

    SJ1s.setError(0.01);
    SJ2s.setError(0.01);
    SJ2c.setError(0.01);
    SJ3.setError(0.01);
    SJ4.setError(0.01);
    DJ5.setError(0.01);
    DJ6s.setError(0.01);
    SJ7.setError(0.01);
    DJ8.setError(0.01);
    DJ9.setError(0.01);
    h1s.setError(0.01);
    h2s.setError(0.01);
    h2c.setError(0.01);
    h3.setError(0.01);
    h4.setError(0.01);
    h5.setError(0.01);
    h6s.setError(0.01);
    h7.setError(0.01);
    h8.setError(0.01);
    h9.setError(0.01);

    timeDependentUntagged model("untagged", "untagged", costhetal, costhetak, phi, t, gamma, y, SJ1s, SJ1c, SJ2s, SJ2c, SJ3, SJ4, DJ5, DJ6s, SJ7, DJ8, DJ9,
                                h1s, h1c, h2s, h2c, h3, h4, h5, h6s, h7, h8, h9);

    auto t_before_toystudy = std::chrono::high_resolution_clock::now();

    unsigned int ngen = 10000;
    unsigned int nruns = 1000;
    RooMCStudy *mcstudy =
            new RooMCStudy(model, RooArgSet(costhetal, costhetak, phi, t),Binned(false), FitOptions(Strategy(2),Hesse(true), Save(true),PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(10)));

    mcstudy->generateAndFit(nruns, ngen, true);

    auto t_after_toystudy = std::chrono::high_resolution_clock::now();

    const RooArgSet* fitParams = mcstudy->fitParams(0);

    TIterator* iter = fitParams->createIterator();
    RooRealVar* var = nullptr;

    while ((var = (RooRealVar*)iter->Next())) {
        const char* name = var->GetName();
        try {
            RooPlot* pullFrame = mcstudy->plotPull(*var, -4, 4);
            if (!pullFrame || pullFrame->numItems() == 0) {
                continue;
            }
            TCanvas* pullCanvas = new TCanvas(Form("c_%s", name), Form("Pull of %s", name), 1600, 1200);
            pullFrame->Draw();
            pullCanvas->SaveAs(Form("pulls2/%s.png", name));
            delete pullCanvas;
        } catch (const std::exception& e) {
            continue;
        }
    }
    delete iter;


    RooPlot* fitFrame = costhetal.frame(100);
    RooDataSet* toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    RooFitResult* fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    TCanvas* fitCanvas = new TCanvas("Fit", "CosThetaL", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/cosThetaL.png");
    delete fitCanvas;

    fitFrame = costhetak.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "CosThetaK", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/cosThetaK.png");
    delete fitCanvas;

    fitFrame = phi.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "Phi", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/phi.png");
    delete fitCanvas;

    fitFrame = t.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "t", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits/t.png");
    delete fitCanvas;

    std::cout << "Processed " << nruns << " toys with " << ngen << " events in " << std::chrono::duration<double, std::milli>(t_after_toystudy-t_before_toystudy).count()/1000 << "s" << std::endl;

}

void runTimeDependentTagged() {
    RooRealVar costhetal("costhetal", "costhetal", -1.0, 1.0);
    RooRealVar costhetak("costhetak", "costhetak", -1.0, 1.0);
    RooRealVar phi("phi", "phi", -TMath::Pi(), +TMath::Pi());
    RooRealVar t("t", "t", 0, 10);

    RooRealVar gamma("gamma", "gamma", 1);
    RooRealVar x("x", "x", 26.93);

    RooRealVar DJ1s("DJ1s", "DJ1s", 0.3, 0, 10);
    RooRealVar DJ2s("DJ2s", "DJ2s", 0, -10, 10);
    RooRealVar DJ2c("DJ2c", "DJ2c", 0, -10, 10);
    RooRealVar DJ3("DJ3", "DJ3", 0, -10, 10);
    RooRealVar DJ4("DJ4", "DJ4", 0, -10, 10);
    RooRealVar SJ5("SJ5", "SJ5", 0, -10, 10);
    RooRealVar SJ6s("SJ6s", "SJ6s", 0, -10, 10);
    RooRealVar DJ7("DJ7", "DJ7", 0, -10, 10);
    RooRealVar SJ8("SJ8", "SJ8", 0, -10, 10);
    RooRealVar SJ9("SJ9", "SJ9", 0, -10, 10);
    RooRealVar s1s("s1s", "s1s", 0, -10, 10);
    RooRealVar s1c("s1c", "s1c", 0, -10, 10);
    RooRealVar s2s("s2s", "s2s", 0, -10, 10);
    RooRealVar s2c("s2c", "s2c", 0, -10, 10);
    RooRealVar s3("s3", "s3", 0, -10, 10);
    RooRealVar s4("s4", "s4", 0, -10, 10);
    RooRealVar s5("s5", "s5", 0, -10, 10);
    RooRealVar s6s("s6s", "s6s", 0, -10, 10);
    RooRealVar s7("s7", "s7", 0, -10, 10);
    RooRealVar s8("s8", "s8", 0, -10, 10);
    RooRealVar s9("s9", "s9", 0, -10, 10);

    RooFormulaVar DJ1c("DJ1c", "(4 + 4*x**2 - 2*s2s*x - s2c*x + 6*s1s*x + 3*s1c*x + 2*DJ2s + DJ2c - 6*DJ1s ) / 3", RooArgList(x, s1s, s1c, s2s, s2c, DJ1s, DJ2s, DJ2c));

    DJ1s.setError(0.01);
    DJ2s.setError(0.01);
    DJ2c.setError(0.01);
    DJ3.setError(0.01);
    DJ4.setError(0.01);
    SJ5.setError(0.01);
    SJ6s.setError(0.01);
    DJ7.setError(0.01);
    SJ8.setError(0.01);
    SJ9.setError(0.01);
    s1s.setError(0.01);
    s2s.setError(0.01);
    s2c.setError(0.01);
    s3.setError(0.01);
    s4.setError(0.01);
    s5.setError(0.01);
    s6s.setError(0.01);
    s7.setError(0.01);
    s8.setError(0.01);
    s9.setError(0.01);

    timeDependentTagged model("tagged", "tagged", costhetal, costhetak, phi, t, gamma, x, DJ1s, DJ1c, DJ2s, DJ2c, DJ3, DJ4, SJ5, SJ6s, DJ7, SJ8, SJ9,
                                s1s, s1c, s2s, s2c, s3, s4, s5, s6s, s7, s8, s9);

    auto t_before_toystudy = std::chrono::high_resolution_clock::now();

    unsigned int ngen = 10000;
    unsigned int nruns = 1000;
    RooMCStudy *mcstudy =
            new RooMCStudy(model, RooArgSet(costhetal, costhetak, phi, t),Binned(false), FitOptions(Strategy(2),Hesse(true), Save(true),PrintEvalErrors(0), EvalBackend("legacy"), NumCPU(10)));

    mcstudy->generateAndFit(nruns, ngen, true);

    auto t_after_toystudy = std::chrono::high_resolution_clock::now();

    const RooArgSet* fitParams = mcstudy->fitParams(0);

    TIterator* iter = fitParams->createIterator();
    RooRealVar* var = nullptr;

    while ((var = (RooRealVar*)iter->Next())) {
        const char* name = var->GetName();
        try {
            RooPlot* pullFrame = mcstudy->plotPull(*var, -4, 4);
            if (!pullFrame || pullFrame->numItems() == 0) {
                continue;
            }
            TCanvas* pullCanvas = new TCanvas(Form("c_%s", name), Form("Pull of %s", name), 1600, 1200);
            pullFrame->Draw();
            pullCanvas->SaveAs(Form("pulls3/%s.png", name));
            delete pullCanvas;
        } catch (const std::exception& e) {
            continue;
        }
    }
    delete iter;


    RooPlot* fitFrame = costhetal.frame(100);
    RooDataSet* toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    RooFitResult* fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    TCanvas* fitCanvas = new TCanvas("Fit", "CosThetaL", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits3/cosThetaL.png");
    delete fitCanvas;

    fitFrame = costhetak.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "CosThetaK", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits3/cosThetaK.png");
    delete fitCanvas;

    fitFrame = phi.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "Phi", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits3/phi.png");
    delete fitCanvas;

    fitFrame = t.frame(100);
    toyData = static_cast<RooDataSet*>(mcstudy->genData(0));
    fitResult = const_cast<RooFitResult *>(mcstudy->fitResult(0));
    model.plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow));
    model.plotOn(fitFrame);
    toyData->plotOn(fitFrame);
    fitCanvas = new TCanvas("Fit", "t", 1600, 1200);
    fitFrame->Draw();
    fitCanvas->SaveAs("fits3/t.png");
    delete fitCanvas;

    std::cout << "Processed " << nruns << " toys with " << ngen << " events in " << std::chrono::duration<double, std::milli>(t_after_toystudy-t_before_toystudy).count()/1000 << "s" << std::endl;

}

int main() {
//    gErrorIgnoreLevel = kFatal + 1;
//    RooMsgService::instance().setStreamStatus(0, false);
//    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL );

    runTimeDependentTagged();

    return 0;
}