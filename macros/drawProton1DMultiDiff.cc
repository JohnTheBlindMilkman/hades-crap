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
    const TString fileName = "../slurmOutput/apr12ana_all_24_10_24_processed.root";
    const TString outputFile = "../output/1Dcorr_50_60_cent.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(200,400)"},{2,"(400,600)"},{3,"(600,800)"},{4,"(800,1000)"},{5,"(1000,1200)"},{6,"(1200,1400)"},{7,"(1400,1600)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.58,-0.35)"},{2,"(-0.35,-0.12)"},{3,"(-0.12,0.12)"},{4,"(0.12,0.35)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 5;

    float norm;    
    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &kt : ktArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignKt%d",kt.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%d",kt.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            norm = getNorm(hRat,300,600);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatKt%d",kt.first));
            hRat->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF(q_{inv})",kt.second.Data()));
            //hRat->SetMarkerStyle(20);
            //hRat->SetMarkerColor(JJColor::fWutAllColors[kt.first-1]);
            prepareGraph(hRat,JJColor::fWutSecondaryColors11[kt.first-1]);

            hRat->GetXaxis()->SetRangeUser(0,490);
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
    canvY->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &y : yArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignY%d",y.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgY%d",y.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            norm = getNorm(hRat,300,900);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatY%d",y.first));
            hRat->SetTitle(TString::Format("y #in %s",y.second.Data()));
            //hRat->SetMarkerStyle(20);
            //hRat->SetMarkerColor(JJColor::fWutAllColors[y.first-1]);
            prepareGraph(hRat,JJColor::fWutSecondaryColors[y.first]);

            hRat->GetXaxis()->SetRangeUser(0,490);
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

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    canvPsi->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &psi : psiArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignPsi%d",psi.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgPsi%d",psi.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            norm = getNorm(hRat,300,900);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatPsi%d",psi.first));
            hRat->SetTitle(TString::Format("#phi - #Psi_{E.P.} #in %s",psi.second.Data()));
            //hRat->SetMarkerStyle(20);
            //hRat->SetMarkerColor(JJColor::fWutAllColors[psi.first-1]);
            prepareGraph(hRat,JJColor::fWutSecondaryColors11[psi.first-1]);

            hRat->GetXaxis()->SetRangeUser(0,490);
            hRat->GetYaxis()->SetRangeUser(0,1.9);

            hRat->Write();
            if (psi.first -1 == 0)
                hRat->Draw("hist pe");
            else
                hRat->Draw("same hist pe");

            hSign->Write();
            hBckg->Write();
        }
    }
    
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();
}