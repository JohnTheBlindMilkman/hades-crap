#include "TFile.h"
#include "TH2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"

#include <tuple>
#include <array>
#include <iostream>

std::pair<double,double> GetEdges(const TH2D *hist, double count)
{
    const int xBins = hist->GetNbinsX();
    const int yBins = hist->GetNbinsY();
    int binX = 0, binY = 0;

    // we assume we are working with CDFs -> bin content is not decreasing with each bin
    for (int i = 1; i <= xBins; ++i)
    {
        double binCont = hist->GetBinContent(i,yBins);
        binX = i;
        if (binCont > count) break;
    }
    for (int i = 1; i <= yBins; ++i)
    {
        double binCont = hist->GetBinContent(xBins,i);
        binY = i;
        if (binCont > count) break;
    }

    return std::make_pair(hist->GetXaxis()->GetBinCenter(binX),hist->GetYaxis()->GetBinCenter(binY));
}

double CalcEntries(const TH2D *hist, double xMin, double yMin, double xMax, double yMax)
{
    int xMinBin = hist->GetXaxis()->FindBin(xMin);
    int xMaxBin = hist->GetXaxis()->FindBin(xMax);
    int yMinBin = hist->GetYaxis()->FindBin(yMin);
    int yMaxBin = hist->GetYaxis()->FindBin(yMax);

    double totEntries = 0;
    for (int x = xMinBin; x <= xMaxBin; ++x)
        for (int y = yMinBin; y <= yMaxBin; ++y)
        {
            totEntries += hist->GetBinContent(x,y);
        }

    return totEntries;
}

void calcEqualBins()
{
    const TString cdfFilePath = "../output/pairDistCDFs.root";
    const TString pdfFilePath = "../slurmOutput/apr12sim_all_25_07_22.root";
    const TString cdfHistName = "hKtRapGoodCent1";
    const TString pdfHistName = "hKtRapGoodCent1";
    constexpr int nIntervals = 8;

    std::array<TLine*,nIntervals - 1> ktLines, rapLines;

    TFile *inpFileCdf = TFile::Open(cdfFilePath);
    TFile *inpFilePdf = TFile::Open(pdfFilePath);

    TH2D *cdfHist = inpFileCdf->Get<TH2D>(cdfHistName);
    double zMax = cdfHist->GetMaximum();
    double step = zMax / nIntervals;

    TH2D *pdfHist = inpFilePdf->Get<TH2D>(pdfHistName);

    double ktMin = pdfHist->GetXaxis()->GetXmin();
    double ktMax = pdfHist->GetXaxis()->GetXmax();
    double rapMin = pdfHist->GetYaxis()->GetXmin();
    double rapMax = pdfHist->GetYaxis()->GetXmax();

    TCanvas *c = new TCanvas("c");
    pdfHist->Draw("colz");

    for (int i = 1; i < nIntervals; ++i)
    {
        auto [ktVal,rapVal] = GetEdges(cdfHist,step * i);
        std::cout << "N = " << step * i << "\tkt = " << ktVal << "\ty = " << rapVal << "\n";
        ktLines.at(i - 1) = new TLine(ktVal,rapMin,ktVal,rapMax);
        rapLines.at(i - 1) = new TLine(ktMin,rapVal,ktMax,rapVal);

        ktLines.at(i - 1)->Draw("same");
        rapLines.at(i - 1)->Draw("same");
    }

    for (int i = 0; i < nIntervals; ++i)
    {
        if (i == 0)
        {
            ktMin = pdfHist->GetXaxis()->GetXmin();
            ktMax = ktLines.at(i)->GetX1();
            rapMin = pdfHist->GetYaxis()->GetXmin();
            rapMax = rapLines.at(i)->GetY1();
        }
        else if (i == nIntervals - 1)
        {
            ktMin = ktLines.at(i - 1)->GetX1();
            ktMax = pdfHist->GetXaxis()->GetXmax();
            rapMin = rapLines.at(i - 1)->GetY1();
            rapMax = pdfHist->GetYaxis()->GetXmax();
        }
        else
        {
            ktMin = ktLines.at(i - 1)->GetX1();
            ktMax = ktLines.at(i)->GetX1();
            rapMin = rapLines.at(i - 1)->GetY1();
            rapMax = rapLines.at(i)->GetY1();
        }

        std::cout << "kT in [" << ktMin << "," << ktMax << "]\ty in [" << rapMin << "," << rapMax << "]\tN = " << CalcEntries(pdfHist,ktMin,rapMin,ktMax,rapMax) << "\n";
    }
    
    TFile *otpFile = TFile::Open("../output/dividedKtRap.root","recreate");
    c->Write();
    pdfHist->Write();
    cdfHist->Write("hKtRapGoodCumulative");
    for (const auto &line : ktLines) 
        line->Write();

    for (const auto &line : rapLines)
        line->Write();
}