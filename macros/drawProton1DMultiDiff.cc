#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "Palettes.hxx"

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

void drawProton1DMultiDiff()
{
    const TString fileName = "/home/jedkol/lxpool/hades-crap/slurmOutput/apr12ana_all_24_02_19_processed.root";
    const TString outputFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_forceEP.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(150,450)"},{2,"(450,750)"},{3,"(750,1050)"},{4,"(1050,1350)"},{5,"(1350,1650)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.75,-0.25)"},{2,"(-0.25,0.25)"},{3,"(0.25,0.75)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 1;

    float norm;    
    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    for (const auto &kt : ktArr)
    {
        TH1D *hSign = inpFile->Get<TH1D>(TString::Format("hQinvSignKt%d",kt.first));
        TH1D  *hBckg = inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%d",kt.first));
        if (hSign != nullptr && hBckg != nullptr)
        {
            TH1D *hRat = new TH1D(*hSign);
            hRat->Divide(hBckg);
            norm = getNorm(hRat,300,900);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->SetName(TString::Format("hQinvRatKt%d",kt.first));
            hRat->SetTitle(TString::Format("k_{T} #in %s",kt.second.Data()));
            hRat->SetMarkerStyle(20);
            hRat->SetMarkerColor(JJColor::fWutAllColors[kt.first-1]);

            hRat->Write();
            if (kt.first -1 == 0)
                hRat->Draw("hist pe pmc plc");
            else
                hRat->Draw("same hist pe pmc plc");
        }
    }
    
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvKt->Write();

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
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
            hRat->SetMarkerStyle(20);
            hRat->SetMarkerColor(JJColor::fWutAllColors[y.first-1]);

            hRat->Write();
            if (y.first -1 == 0)
                hRat->Draw("hist pe pmc plc");
            else
                hRat->Draw("same hist pe pmc plc");
        }
    }
    
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvY->Write();

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
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
            hRat->SetMarkerStyle(20);
            hRat->SetMarkerColor(JJColor::fWutAllColors[psi.first-1]);

            hRat->Write();
            if (psi.first -1 == 0)
                hRat->Draw("hist pe pmc plc");
            else
                hRat->Draw("same hist pe pmc plc");
        }
    }
    
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();
}