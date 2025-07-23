#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"

#include <iostream>

#include "../Externals/indicators.hpp"

double GetBinCenterByConent(const TH1D *hist, double content)
{
    const int xBins = hist->GetNbinsX();
    int bin = 0;

    // we assume we are working with CDFs -> bin content is not decreasing with each bin
    for (int i = 1; i <= xBins; ++i)
    {
        double binCont = hist->GetBinContent(i);
        bin = i;
        if (binCont > content) break;
    }

    return hist->GetBinCenter(bin);
}

TH2D* CalcCDF(const TH2D *hist)
{
    TH2D *hCumulative = dynamic_cast<TH2D*>(hist->Clone());
    const int xBins = hCumulative->GetXaxis()->GetNbins();
    const int yBins = hCumulative->GetYaxis()->GetNbins();

    indicators::BlockProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::ShowPercentage{true},
                    indicators::option::PrefixText{"Calculating CDF"},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{xBins * yBins}};

    indicators::show_console_cursor(false);

    for (int xBin = 1; xBin <= xBins; ++xBin)
    {
        for (int yBin = 1; yBin <= yBins; ++yBin)
        {
            double binConent = 0;
            for (int i = 1; i <= xBin; ++i)
            {
                for (int j = 1; j <= yBin; ++j)
                {
                    binConent += hist->GetBinContent(i,j);
                }
            }
            hCumulative->SetBinContent(xBin,yBin,binConent);
            bar.tick();
        }
    }

    bar.mark_as_completed();
    indicators::show_console_cursor(true);

    return hCumulative;
}

void generateCdfs()
{
    const TString filePath = "../slurmOutput/apr12sim_all_25_07_22.root";
    const TString ktRapName = "hKtRapGoodCent";
    const TString ktCentNameBase = "ktDist_cent";
    const TString rapCentNameBase = "rapDist_cent";
    const TString nameEnd = "_good";
    const std::vector<int> centralities = {1,2,3,4};
    constexpr double entries = 6e8;

    std::vector<TH1D*> ktCdfs, rapCdfs;
    std::vector<TH2D*> ktRapCdfs;
    double maximum = 0, xValue = 0;
    int intervals = 0;

    TFile *inpFile = TFile::Open(filePath);
    

    for (const auto &cent : centralities)
    {
        TH2D *hPdfKtRap = inpFile->Get<TH2D>(TString::Format("%s%d",ktRapName.Data(),cent));
        TH2D *hCdfKtRap = CalcCDF(hPdfKtRap);
        ktRapCdfs.push_back(hCdfKtRap);

        TH1D *hKtCentPdf = inpFile->Get<TH1D>(TString::Format("%s%d%s",ktCentNameBase.Data(),cent,nameEnd.Data()));
        TH1D *hKtCentCdf = dynamic_cast<TH1D*>(hKtCentPdf->GetCumulative());

        maximum = hKtCentCdf->GetMaximum();
        intervals = static_cast<int>(maximum / entries);

        std::cout << "Minimum entries required: " << entries << "\n";
        std::cout << "Generating " << intervals << " kT intervals...\n";

        for (int i = 0; i <= intervals; ++i)
        {
            xValue = GetBinCenterByConent(hKtCentCdf,entries * i);
            std::cout << "kT = " << xValue << " MeV/c\n";
        }
        ktCdfs.push_back(hKtCentCdf);

        TH1D *hRapCentPdf = inpFile->Get<TH1D>(TString::Format("%s%d%s",rapCentNameBase.Data(),cent,nameEnd.Data()));
        TH1D *hRapCentCdf = dynamic_cast<TH1D*>(hRapCentPdf->GetCumulative());

        maximum = hRapCentCdf->GetMaximum();
        intervals = static_cast<int>(maximum / entries);

        std::cout << "Minimum entries required: " << entries << "\n";
        std::cout << "Generating " << intervals << " y intervals...\n";

        for (int i = 0; i <= intervals; ++i)
        {
            xValue = GetBinCenterByConent(hRapCentCdf,entries * i);
            std::cout << "y12 = " << xValue << "\n";
        }
        rapCdfs.push_back(hRapCentCdf);
    }

    TFile *otpFile = TFile::Open("../output/pairDistCDFs.root","recreate");

    for (const auto &hist : ktCdfs) 
        hist->Write();

    for (const auto &hist : rapCdfs) 
        hist ->Write();

    for (const auto &hist : ktRapCdfs) 
        hist ->Write();
}