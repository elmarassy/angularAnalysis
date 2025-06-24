//
// Created by Mero Elmarassy on 6/18/25.
//

#include <TF1.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooMCStudy.h>
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2.h"
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
#include <TLatex.h>
#include <map>

#include "../include/plotting.h"
#include "../include/parameterValues.h"


std::pair<std::pair<double, double>, std::pair<double, double>> makeHist(const char* title, const std::vector<double>& values, double trueValue, bool fitGauss) {
    double fitMean;
    double fitSigma;
    double errMean;
    double errSigma;

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
        errMean = mean.getError();
        fitSigma = sigma.getVal();
        errSigma = sigma.getError();
    }
    if (trueValue > mid - 1.5*spread && trueValue < mid + 1.5*spread) {
        TLine *line = new TLine(trueValue, 0, trueValue, height);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw("same");
    }
    return std::make_pair(std::make_pair(fitMean, errMean), std::make_pair(fitSigma, errSigma));
}

void plotSummary(const char* title, const std::vector<std::pair<double, double>>& results,
                      const std::vector<std::string>& labels, bool showTrueValues=false) {
    int n = results.size();

    TH1D *h = new TH1D(title, title, n, 0, n);

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
        TH1D *trueVals = new TH1D("trueVals", "", n, 0, n);

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

void plotValueSummary(const char* title, const std::vector<std::pair<double, double>>& results,
                         const std::vector<std::string>& labels) {
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


    TLine* line = new TLine(0, 0, n, 0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw("same");
}



void plotPullSummary(const char* title, const std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>>& results,
                 const std::vector<std::string>& labels) {
    int n = results.size();
    gStyle->SetEndErrorSize(8);
    TH1D* h = new TH1D("mean", title, 2*n, 0, n);

    for (int i = 0; i < n; ++i) {
        h->SetBinContent(2*i + 1, results[i].first.first);
        h->SetBinError(2*i + 1, results[i].first.second);
        h->GetXaxis()->SetBinLabel(2*i + 1, labels[i].c_str());
        h->GetXaxis()->SetBinLabel(2*i, "");
    }

    h->SetStats(false);
    h->SetMarkerStyle(20);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0.8);
    h->Draw("E1P");

    TH1D* sigma = new TH1D("sigma", "sigma", 2*n, 0, n);

    for (int i = 0; i < n; ++i) {
        sigma->SetBinContent(2*i + 2, results[i].second.first - 1);
        sigma->SetBinError(2*i + 2, results[i].second.second);
        h->GetXaxis()->SetBinLabel(2*i + 1, labels[i].c_str());
        h->GetXaxis()->SetBinLabel(2*i, "");
    }

    sigma->SetStats(false);
    gStyle->SetEndErrorSize(8);
    sigma->SetMarkerStyle(20);
    sigma->SetMarkerColor(kRed);
    sigma->SetLineColor(kRed);
    sigma->SetMarkerSize(0.8);
    sigma->Draw("same E1P");

    TLine* line = new TLine(0, 0, n, 0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw("same");
}

void plotPearsonCorrelation(const std::map<std::string, std::pair<std::vector<double>, std::vector<double>>>& results) {

    auto mean = [&](const std::vector<double> values) {
        double sum = 0;
        for (auto value: values) sum += value;
        return sum / values.size();
    };

    auto cov = [&](const std::vector<double> values1, const std::vector<double> values2) {
        double sum = 0;
        double size = values1.size();
        for (int i = 0; i < size; i++) sum += values1.at(i) * values2.at(i);
        return sum / size - mean(values1) * mean(values2);
    };
    auto corr = [&](const std::vector<double> values1, const std::vector<double> values2) {
        auto covariance = cov(values1, values2);
        auto sigma1 = cov(values1, values1);
        auto sigma2 = cov(values2, values2);
        return covariance / std::sqrt(sigma1*sigma2);
    };

    auto size = results.size();
    TH2F* correlations = new TH2F("Pearson correlation matrix", "", size, 0, size, size, 0, size);

    int i = 1;
    for (const auto& entry: results) {
        int j = 1;
        correlations->GetXaxis()->SetBinLabel(i, entry.first.c_str());
        for (auto otherEntry = results.begin(); otherEntry != results.upper_bound(entry.first); otherEntry++) {
            int bin = correlations->GetBin(i, size + 1 - j);
            double correlation = corr(entry.second.first, otherEntry->second.first);
            correlations->SetBinContent(bin, correlation);
            correlations->GetYaxis()->SetBinLabel(size + 1 - j, otherEntry->first.c_str());
            j++;
        }
        i++;
    }
    correlations->SetMinimum(-1);
    correlations->SetMaximum(1);
    plotCorrelationMatrix(correlations);
}

void plotPulls(std::map<std::string, std::pair<std::vector<double>, std::vector<double>>>& results, TFile* outFile) {
    std::map<std::string, std::vector<std::string>> groups;
    outFile->cd();
    TDirectory* pullDir = outFile->mkdir("pulls");
    TDirectory* valueDir = outFile->mkdir("values");
    TCanvas* pullSummary = new TCanvas("Pull summary",
                                   "Pull summary",
                                   4800, 3200);
    pullSummary->Divide(2, 2);
    TCanvas* valueSummary = new TCanvas("Value summary",
                                       "Value summary",
                                       4800, 3200);
    valueSummary->Divide(2, 2);
    std::vector<std::string> removeKeys;

    for (const auto& [param, hist] : results) {
        std::string prefix(1, param[0]);
        if (std::isupper(prefix[0]) && !hist.first.empty()) {
            groups[prefix].push_back(param);
        }
        if (hist.first.empty()) {
            removeKeys.push_back(param);
        }
    }
    for (auto key: removeKeys) {
        results.erase(key);
    }
    std::cout << "Done removing" << std::endl;


    int summaryIndex = 1;
    for (const auto& [type, params]: groups) {
        std::vector<std::pair<std::pair<double, double>, std::pair<double, double>>> pullResults;
        std::vector<std::pair<double, double>> valueResults;
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
            auto value = makeHist(Form("Values for %s", param.c_str()), hist.first, testValues.at(param), true);
            valueResults.push_back(std::make_pair(value.first.first, value.second.first));
        }
        pullSummary->cd(summaryIndex);
        plotPullSummary(Form("%s", type.c_str()), pullResults, names);
        valueSummary->cd(summaryIndex++);
        plotValueSummary(Form("%s", type.c_str()), valueResults, names);
        pulls->Update();
        values->Update();
        pullDir->cd();
        pulls->Write();
        valueDir->cd();
        values->Write();
    }
    pullSummary->Update();
    valueSummary->Update();
    outFile->cd();
    pullSummary->Write();
    valueSummary->Write();
    plotPearsonCorrelation(results);
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


void plotCorrelationMatrix(TH2* corrHist) {
    TCanvas* c = new TCanvas("Correlation Matrix", "", 1800, 1600);
    corrHist->SetTitle("");
    corrHist->GetXaxis()->SetLabelSize(0.025);
    corrHist->GetYaxis()->SetLabelSize(0.025);
    corrHist->GetZaxis()->SetLabelSize(0.03);
    corrHist->SetMarkerSize(0.5);
    gPad->SetRightMargin(0.15);
    corrHist->Draw("COLZ");
    int nBinsX = corrHist->GetNbinsX();
    int nBinsY = corrHist->GetNbinsY();

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            double val = corrHist->GetBinContent(i, j);
            TString text = Form("%.3f", val);
            double x = corrHist->GetXaxis()->GetBinCenter(i);
            double y = corrHist->GetYaxis()->GetBinCenter(j);
            TLatex* latex = new TLatex(x, y, text);
            latex->SetTextSize(0.01);
            latex->SetTextAlign(22);
            latex->Draw();
        }
    }

    TLegend* legend = new TLegend(0.15, 0.85, 0.5, 0.9);
    legend->SetHeader("", "C");
    legend->AddEntry(corrHist, "", "f");
    legend->SetBorderSize(1);
    legend->SetFillStyle(0);
    legend->Draw();

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
    gStyle->SetOptTitle(0);
    plotCorrelationMatrix(fitResult->correlationHist());
    gStyle->SetOptTitle(1);
}

