//
// Created by Mero Elmarassy on 6/18/25.
//

#include <TF1.h>
#include <TGraphErrors.h>
#include <TLine.h>
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
#include <RooRealVar.h>
#include <RooGaussian.h>

#include "../include/plotting.h"
#include "../include/parameterValues.h"


std::pair<double, double> makeHist(const char* title, const std::vector<double>& values, double trueValue, bool fitGauss) {
    double fitMean;
    double fitSigma;

    double upper = *std::max_element(values.begin(), values.end());
    double lower = *std::min_element(values.begin(), values.end());
    double mid = (upper + lower) / 2;
    double spread = (upper - lower) / 2;
    auto hist = new TH1F(title, title, 100, mid - 1.5*spread, mid + 1.5*spread);
    for (auto value: values) {
        hist->Fill(value);
    }
    hist->SetLineColor(kBlue);
    hist->SetFillColorAlpha(kBlue, 0.3);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(kBlue);
    hist->SetTitle(hist->GetName());
    hist->GetXaxis()->SetTitle(title);
    hist->Draw("hist");

    double height = 0;
    if (fitGauss) {
        RooRealVar x("x", "x", mid - 1.5*spread, mid + 1.5*spread);
        RooRealVar mean("mean", "mean", hist->GetMean(), lower, upper);
        RooRealVar sigma("sigma", "sigma", hist->GetRMS(), 0.001, 2 * hist->GetRMS());
        RooGaussian gauss("gauss", "gauss", x, mean, sigma);
        RooDataSet data("data", "data", RooArgSet(x));
        for (double val : values) {
            x.setVal(val);
            data.add(RooArgSet(x));
        }
        gauss.fitTo(data, RooFit::PrintLevel(-1));
        RooPlot* frame = x.frame();
        gauss.plotOn(frame, RooFit::Normalization(data.sumEntries(), RooAbsReal::NumEvent));
        frame->Draw("same");
        gPad->Update();
        height = frame->GetMaximum();
        TPaveText *pave = new TPaveText(0.65, 0.7, 0.9, 0.9, "NDC");
        pave->AddText(Form("mean = %.3f +/- %.3f", mean.getVal(), mean.getError()));
        pave->AddText(Form("#sigma = %.3f +/- %.3f", sigma.getVal(), sigma.getError()));
        pave->Draw("same");
        fitMean = mean.getVal();
        fitSigma = sigma.getVal();
    }
    if (trueValue > mid - 1.5*spread && trueValue < mid + 1.5*spread) {
        TLine *line = new TLine(trueValue, 0, trueValue, height);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }
    return std::make_pair(fitMean, fitSigma);
}


void plotSummary(const char* title, const std::vector<std::pair<double, double>>& results,
                         const std::vector<std::string>& labels, bool showTrueValues=false) {
    int n = results.size();

    TH1D* h = new TH1D(title, title, n, 0, n);

    for (int i = 0; i < n; ++i) {
        h->SetBinContent(i + 1, results[i].first);
        h->SetBinError(i + 1, results[i].second);
        h->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
    }

    h->SetStats(0);
    gStyle->SetEndErrorSize(8);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kGray);
    h->SetMarkerSize(0.8);
    h->Draw("E1P");

    if (showTrueValues) {
        TH1D* trueVals = new TH1D("trueVals", "", n, 0, n);

        for (int i = 0; i < n; ++i) {
            trueVals->SetBinContent(i + 1, testValues.at(labels[i]));
            trueVals->SetBinError(i + 1, 0.001);
            h->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
        }

        trueVals->SetStats(0);
        trueVals->SetMarkerStyle(20);
        trueVals->SetMarkerColor(kRed);
        trueVals->SetLineColor(kRed);
        trueVals->SetMarkerSize(0);
        trueVals->Draw("same");
    }

    TLine* line = new TLine(0, 0, n, 0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw("same");
}


void plotPulls(const std::map<std::string, std::pair<std::vector<double>, std::vector<double>>>& results, TFile* outFile) {
    std::map<std::string, std::vector<std::string>> groups;
    outFile->cd();
    TDirectory* pullDir = outFile->mkdir("pulls");
    TDirectory* valueDir = outFile->mkdir("values");
    TCanvas* summary = new TCanvas("Summary",
                                   "Pull summary",
                                   4800, 3200);
    summary->Divide(2, 2);

    for (const auto& [param, hist] : results) {
        std::string prefix(1, param[0]);
        if (std::isupper(prefix[0]) && !hist.first.empty()) {
            groups[prefix].push_back(param);
        }
    }


    int summaryIndex = 1;
    for (const auto& [type, params]: groups) {
        std::vector<std::pair<double, double>> pullResults;
        std::vector<std::string> names;
        TCanvas* pulls = new TCanvas(Form("pulls_%s", type.c_str()),
                                 Form("Pulls for %s parameters", type.c_str()),
                                 4800, 3200);
        TCanvas* values = new TCanvas(Form("values_%s", type.c_str()),
                                 Form("Values for %s parameters", type.c_str()),
                                 4800, 3200);
        int n = params.size();
        int nCols = 3;
        int nRows = (n + nCols - 1) / nCols;
        pulls->Divide(nCols, nRows);
        values->Divide(nCols, nRows);

        int padPulls = 1;
        int padVals = 1;
        for (const auto& param : params) {
            const auto& hist = results.at(param);
            pulls->cd(padPulls++);
            std::vector<double> pullValues;
            pullValues.reserve(hist.first.size());
            for (int i = 0; i < hist.first.size(); i++) {
                pullValues.push_back((hist.first[i] - testValues.at(param))/hist.second[i]);
            }
            pullResults.push_back(makeHist(Form("Pull for %s", param.c_str()), pullValues, 0, true));
            names.push_back(param);
            values->cd(padVals++);
            makeHist(Form("Values for %s", param.c_str()), hist.first, testValues.at(param), false);
        }
        summary->cd(summaryIndex++);
        plotSummary(Form("%s", type.c_str()), pullResults, names);
        pulls->Update();
        values->Update();
        pullDir->cd();
        pulls->Write();
        valueDir->cd();
        values->Write();
    }
    summary->Update();
    outFile->cd();
    summary->Write();
}

void plot(const std::vector<RooRealVar*>& observables, RooDataSet* fitData, RooFitResult* result, RooAbsPdf& pdf, char* title) {
    using namespace RooFit;
    auto helper = [&](RooRealVar* obs) {
        RooPlot* fitFrame = obs->frame(strcmp(obs->GetName(), "t") == 0? 100: 50);
        auto* fitResult = const_cast<RooFitResult *>(result);
        RooArgSet projSet;
        for (auto otherObs: observables) {
            if (otherObs != obs) projSet.add(*otherObs);
        }
        auto proj = pdf.createProjection(projSet);
        proj->plotOn(fitFrame, VisualizeError(*fitResult, 1), FillColor(kYellow), Normalization(fitData->sumEntries(), RooAbsReal::NumEvent));
        proj->plotOn(fitFrame, Normalization(fitData->sumEntries(), RooAbsReal::NumEvent));
        fitFrame->SetMarkerSize(0.2);
        fitData->plotOn(fitFrame);
        fitFrame->SetTitle(obs->GetName());
        fitFrame->Draw();
    };
    TCanvas* c = new TCanvas(title, title, 4800, 3200);
    c->Divide(2, 2);
    int index = 1;
    for (auto obs: observables) {
        c->cd(index++);
        helper(obs);
    }
    c->Update();
    c->Write();
}


void makePullProjections(std::map<std::string, std::pair<double, double>> currentPullResults, RooFitResult* fitResult, RooSimultaneous& pdf, RooDataSet* fitData, TDirectory* pullDir, int toyNumber) {
    std::map<std::string, std::vector<std::string>> groups;

    pullDir->cd();

    TCanvas* valueSummary = new TCanvas("Value summary",
                                   "Value summary",
                                   4800, 3200);

    valueSummary->Divide(2, 2);
    TCanvas* pullSummary = new TCanvas("Pull summary",
                                        "Pull summary",
                                        4800, 3200);

    pullSummary->Divide(2, 2);
    int summaryIndex = 1;
    for (const auto& pair: testValues) {
        if (currentPullResults.find(pair.first) != currentPullResults.end()) {
            std::string prefix(1, pair.first[0]);
            groups[prefix].push_back(pair.first);
        }
    }

    for (const auto& group: groups) {
        std::vector<std::pair<double, double>> groupValueSummary;
        std::vector<std::pair<double, double>> groupPullSummary;

        for (const auto& param: group.second) {
            groupValueSummary.push_back(currentPullResults.at(param));
            auto temp = std::pair<double, double>((currentPullResults.at(param).first - testValues.at(param) ) / currentPullResults.at(param).second, 0);
            groupPullSummary.push_back(temp);
        }
        pullSummary->cd(summaryIndex);
        plotSummary(Form("%s", group.first.c_str()), groupPullSummary, group.second, false);
        valueSummary->cd(summaryIndex++);
        plotSummary(Form("%s", group.first.c_str()), groupValueSummary, group.second, true);
    }
    pullSummary->Write();
    valueSummary->Write();

    auto* obsList = pdf.getPdf("B0")->getObservables(fitData);
    std::vector<RooRealVar*> observables;
    observables.push_back(dynamic_cast<RooRealVar*>(obsList->find("cosThetaL")));
    observables.push_back(dynamic_cast<RooRealVar*>(obsList->find("cosThetaK")));
    observables.push_back(dynamic_cast<RooRealVar*>(obsList->find("phi")));
    observables.push_back(dynamic_cast<RooRealVar*>(obsList->find("t")));
    RooDataSet* bData = (RooDataSet*) fitData->reduce("type == type::B0");
    RooDataSet* bBarData = (RooDataSet*) fitData->reduce("type == type::B0Bar");

    plot(observables, bData, fitResult, *pdf.getPdf("B0"), "B0 projections");
    plot(observables, bBarData, fitResult, *pdf.getPdf("B0Bar"), "B0Bar projections");
}

