#include "TFile.h"
#include "TH2D.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TString.h"

#include <tuple>
#include <array>
#include <iostream>

template <typename T>
std::pair<T,T> Factorise(T n)
{
    if (n <= 0) // this algorithm doesn't work with negative numbers
        return {};

    T a, b, b1;
    if ((n % 2) == 0) // if n is even
    {
        a = n / 2;
        b = 2;
        while (a > b) // try to minimise the difference between the two numbers, e.g. instead of 98 -> (49,2) do 98 -> (14,7)
        {
            std::tie(a,b1) = Factorise(a);
            b *= b1;
        }
        return std::make_pair(std::max(a,b),std::min(a,b));
    }

    a = std::ceil(std::sqrt(n));
    if (a * a == n) // if n is a square of some number
    {
        return std::make_pair(a,a);
    }

    while (true) // Fermat's factorisation (only for odd numbers!)
    {
        b1 = a * a - n;
        b = static_cast<T>(std::sqrt(b1));

        if (b * b == b1)
        {
            break;
        }
        else
        {
            ++a;
        }
    }
    
    return std::make_pair(a + b, a - b);
}

std::pair<double,double> GetEdges(const TH2D *hist, double count)
{
    const int xBins = hist->GetNbinsX();
    const int yBins = hist->GetNbinsY();
    std::pair<int,double> binX{0,0.}, binY{0,0.};

    // we assume we are working with CDFs -> bin content is not decreasing with each bin
    for (int i = 1; i <= xBins; ++i)
    {
        std::pair<int,double> binXnew{i,hist->GetBinContent(i,yBins)};
        
        if (std::abs(binX.second - count) < std::abs(binXnew.second - count)) 
            break;
        binX = binXnew;
    }
    for (int i = 1; i <= yBins; ++i)
    {
        std::pair<int,double> binYnew{i,hist->GetBinContent(xBins,i)};
        
        if (std::abs(binY.second - count) < std::abs(binYnew.second - count)) 
            break;
        binY = binYnew;
    }   

    return std::make_pair(hist->GetXaxis()->GetBinCenter(binX.first),hist->GetYaxis()->GetBinCenter(binY.first));
}

double GetBinCenterByConent(const TH1D *hist, double content)
{
    const int xBins = hist->GetNbinsX();
    std::pair<int,double> bin{0,0.};

    // we assume we are working with CDFs -> bin content is not decreasing with each bin
    for (int i = 1; i <= xBins; ++i)
    {
        std::pair<int,double> binNew{i,hist->GetBinContent(i)};
        
        if (std::abs(bin.second - content) < std::abs(binNew.second - content)) 
            break;
        bin = binNew;
    }

    return hist->GetBinLowEdge(bin.first + 1);
}

std::pair<std::vector<TLine*>,std::vector<TLine*> > CreateLines(const TH2D *hist, std::size_t xIntervals, std::size_t yIntervals)
{
    const double zMax = hist->GetMaximum();
    const double ktStep = zMax / xIntervals;
    const double rapStep = zMax / yIntervals;
    const double ktMin = hist->GetXaxis()->GetXmin();
    const double ktMax = hist->GetXaxis()->GetXmax();
    const double rapMin = hist->GetYaxis()->GetXmin();
    const double rapMax = hist->GetYaxis()->GetXmax();

    std::vector<TLine*> ktLines, yLines;

    for (std::size_t i = 1; i < xIntervals; ++i)
    {
        double ktVal = GetBinCenterByConent(hist->ProjectionX("_px",hist->GetNbinsY(),hist->GetNbinsY()),ktStep * i);
        ktLines.push_back(new TLine(ktVal,rapMin,ktVal,rapMax));
    }
    for (std::size_t i = 1; i < yIntervals; ++i)
    {
        double rapVal = GetBinCenterByConent(hist->ProjectionY("_py",hist->GetNbinsX(),hist->GetNbinsX()),rapStep * i);
        yLines.push_back(new TLine(ktMin,rapVal,ktMax,rapVal));
    }

    return std::make_pair(ktLines,yLines);
}

double CalcEntries(const TH2D *hist, double xMin, double yMin, double xMax, double yMax)
{
    int xMinBin = hist->GetXaxis()->FindFixBin(xMin);
    int xMaxBin = hist->GetXaxis()->FindFixBin(xMax);
    int yMinBin = hist->GetYaxis()->FindFixBin(yMin);
    int yMaxBin = hist->GetYaxis()->FindFixBin(yMax);

    return hist->GetBinContent(xMaxBin,yMaxBin) - hist->GetBinContent(xMinBin,yMaxBin) - hist->GetBinContent(xMaxBin,yMinBin) + hist->GetBinContent(xMinBin,yMinBin);
}

void calcEqualBins()
{
    const TString cdfFilePath = "../output/pairDistCDFs.root";
    const TString pdfFilePath = "../slurmOutput/apr12ana_all_25_07_25.root";
    const TString cdfHistName = "hKtRapGoodCent1";
    const TString pdfHistName = "hKtRapGoodCent1";
    constexpr double entries = 2e9;

    std::vector<TLine*> ktLines, rapLines;

    TFile *inpFileCdf = TFile::Open(cdfFilePath);
    TFile *inpFilePdf = TFile::Open(pdfFilePath);

    TH2D *cdfHist = inpFileCdf->Get<TH2D>(cdfHistName);
    TH2D *pdfHist = inpFilePdf->Get<TH2D>(pdfHistName);
    double zMax = cdfHist->GetMaximum();
    std::size_t nSections = static_cast<std::size_t>(std::floor(zMax / entries));
    auto [ktIntervals,rapIntervals] = Factorise(nSections);
    std::cout << "Sections: " << nSections << "\tkt: " << ktIntervals << "\ty: " << rapIntervals << "\n";

    TCanvas *c = new TCanvas("c");
    pdfHist->Draw("colz");

    std::tie(ktLines, rapLines) = CreateLines(cdfHist,ktIntervals,rapIntervals);

    for (const auto &hist : ktLines)
        hist->Draw("same");
    for (const auto &hist : rapLines)
        hist->Draw("same");

    for (int i = 0; i < ktIntervals; ++i)
    {
        for (int j = 0; j < rapIntervals; ++j)
        {
            const double ktMin = (i == 0) ? pdfHist->GetXaxis()->GetXmin() : ktLines.at(i - 1)->GetX1();
            const double ktMax = (i == ktIntervals - 1) ? pdfHist->GetXaxis()->GetXmax() - pdfHist->GetXaxis()->GetBinWidth(1) : ktLines.at(i)->GetX1();
            const double rapMin = (j == 0) ? pdfHist->GetYaxis()->GetXmin() : rapLines.at(j - 1)->GetY1();
            const double rapMax = (j == rapIntervals - 1) ? pdfHist->GetYaxis()->GetXmax() - pdfHist->GetYaxis()->GetBinWidth(1) : rapLines.at(j)->GetY1();

            std::cout << "kT in [" << ktMin << "," << ktMax << "]\ty in [" << rapMin + 0.74 << "," << rapMax + 0.74 << "]\tN = " << CalcEntries(cdfHist, ktMin, rapMin, ktMax, rapMax) << "\n";
        }
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