#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"
#include "../FemtoMixer/PairUtils.hxx"

double getNorm(const TH1D *hInp, double xMin, double xMax)
{
    int nBins = 0;
    double val = 0., xVal;
    for (int i = 1; i < hInp->GetNbinsX(); i++)
    {
        xVal = hInp->GetBinCenter(i);
        if (xVal >= xMin && xVal <= xMax)
        {
            val += hInp->GetBinContent(i);
            nBins++;
        } 
    }
    if(nBins > 0)
        return val/nBins;
    else    
        return 0.;
}

void prepareGraph(TH1D* hist, int col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    //hist->GetXaxis()->SetLabelSize();
    //hist->GetXaxis()->SetLabelOffset();
    //hist->GetXaxis()->SetTitleSize();
    //hist->GetXaxis()->SetTitleOffset();

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
}

void drawProton1DMultiDiff()
{
    const TString fileName = "../slurmOutput/apr12ana_all_25_09_24_processed.root";
    const TString outputFile = "../output/1Dcorr_30_40_cent.root";
    const auto ktArr = Mixing::PairGrouping{}.GetKtIndexIntervalPairs();
    const auto yArr = Mixing::PairGrouping{}.GetRapIndexIntervalPairs();
    const int rebin = 1;
    constexpr int minX = 0, maxX = 300;

    float norm;    
    TLine *line = new TLine(minX,1,maxX,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.15,0.05,0.15,0.05);
    for (const auto &kt : ktArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignKt%ld",kt.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%ld",kt.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            hRat->Sumw2();
            norm = JJUtils::CF::GetNormByRange(hRat,200,300);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatKt%ld",kt.first));
            hRat->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF(q_{inv})",kt.second.Data()));
            prepareGraph(hRat,JJColor::fWutSecondaryColors11[kt.first-1]);

            hRat->GetXaxis()->SetRangeUser(minX,maxX);
            hRat->GetYaxis()->SetRangeUser(0,1.9);

            hRat->Write();
            if (kt.first -1 == 0)
                hRat->Draw("hist pe");
            else
                hRat->Draw("same hist pe");

            hSign->Write();
            hBckg->Write();
        }
    }
    
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvKt->Write();

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    canvY->SetMargin(0.15,0.05,0.15,0.05);
    for (const auto &y : yArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignY%ld",y.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgY%ld",y.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            hRat->Sumw2();
            norm = JJUtils::CF::GetNormByRange(hRat,200,300);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatY%ld",y.first));
            hRat->SetTitle(TString::Format("y_{lab}^{1,2} #in %s",y.second.Data()));
            prepareGraph(hRat,JJColor::fWutSecondaryColors11[y.first]);

            hRat->GetXaxis()->SetRangeUser(minX,maxX);
            hRat->GetYaxis()->SetRangeUser(0,1.9);

            hRat->Write();
            if (y.first -1 == 0)
                hRat->Draw("hist pe");
            else
                hRat->Draw("same hist pe");

            hSign->Write();
            hBckg->Write();
        }
    }
    
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvY->Write();

    for (const auto &kt : ktArr)
        for (const auto &y : yArr)
        {
            TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignKt%ldY%ld",kt.first,y.first));
            TH1D *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%ldY%ld",kt.first,y.first));
            if (hSign != nullptr && hBckg != nullptr)
            {
                TH1D *hRat = new TH1D(*hSign);
                hRat->Divide(hBckg);
                hRat->Sumw2();
                norm = JJUtils::CF::GetNormByRange(hRat,200,300);
                hRat->Rebin(rebin);
                norm *= rebin;
                hRat->Scale(1./norm);
                hRat->SetName(TString::Format("hQinvRatKt%ldY%ld",kt.first,y.first));
                hRat->SetTitle(TString::Format("k_{T} #in %s and y_{lab}^{1,2} #in %s;q_{inv} [MeV/c];CF(q_{inv})",kt.second.Data(),y.second.Data()));
                prepareGraph(hRat,JJColor::fWutSecondaryColors11[kt.first-1]);

                hRat->GetXaxis()->SetRangeUser(minX,maxX);
                hRat->GetYaxis()->SetRangeUser(0,1.9);

                hRat->Write();
                hSign->Write();
                hBckg->Write();
            }
        }
}