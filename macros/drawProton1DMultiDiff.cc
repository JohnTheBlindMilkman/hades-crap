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
    const TString fileName = "/home/jedkol/lxpool/hades-crap/slurmOutput/apr12ana_all_24_02_13.root";
    const TString outputFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_90deg.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(150,450)"},{2,"(450,750)"},{3,"(750,1050)"},{4,"(1050,1350)"},{5,"(1350,1650)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.75,-0.25)"},{2,"(-0.25,0.25)"},{3,"(0.25,0.75)"}};
    //const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-125,-115)"},{2,"(-115,-105)"},{3,"(-105,-95)"},{4,"(-95,-85)"},{5,"(-85,-75)"},{6,"(-75,-65)"},{7,"(-65,-55)"}};
    const int rebin = 1;

    float norm;
    std::vector<std::vector<std::vector<TH1D*> > > 
    hSign(ktArr.size(),std::vector<std::vector<TH1D*> >(yArr.size(),std::vector<TH1D*>(psiArr.size(),nullptr))), 
    hBckg(ktArr.size(),std::vector<std::vector<TH1D*> >(yArr.size(),std::vector<TH1D*>(psiArr.size(),nullptr)));

    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);

    for (const auto &kt : ktArr)
        for (const auto &y : yArr)
            for (const auto &psi : psiArr)
            {
                hSign[kt.first-1][y.first-1][psi.first-1] = inpFile->Get<TH1D>(TString::Format("hQinvSign_%d%d%d",kt.first,y.first,psi.first));
                if (hSign[kt.first-1][y.first-1][psi.first-1] != nullptr)
                    hSign[kt.first-1][y.first-1][psi.first-1]->Sumw2();

                hBckg[kt.first-1][y.first-1][psi.first-1] = inpFile->Get<TH1D>(TString::Format("hQinvBckg_%d%d%d",kt.first,y.first,psi.first));
                if (hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                    hBckg[kt.first-1][y.first-1][psi.first-1]->Sumw2();
            }

    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    std::vector<TH1D*> hRatKt(ktArr.size(),nullptr), hBckgKt(ktArr.size(),nullptr);
    for (const auto &kt : ktArr)
    {
        for (const auto &y : yArr)
            for (const auto &psi : psiArr)
            {
                if (hRatKt[kt.first-1] == nullptr && hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckgKt[kt.first-1] == nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatKt[kt.first-1] = new TH1D(*hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgKt[kt.first-1] = new TH1D(*hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
                else if(hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatKt[kt.first-1]->Add(hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgKt[kt.first-1]->Add(hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
            }

        if (hRatKt[kt.first-1] != nullptr && hBckgKt[kt.first-1] != nullptr)
        {
            hRatKt[kt.first-1]->Divide(hBckgKt[kt.first-1]);
            norm = getNorm(hRatKt[kt.first-1],300,900);
            hRatKt[kt.first-1]->Rebin(rebin);
            norm *= rebin;
            hRatKt[kt.first-1]->Scale(1./norm);
            hRatKt[kt.first-1]->SetName(TString::Format("hQinvRatKt%d",kt.first));
            hRatKt[kt.first-1]->SetTitle(TString::Format("k_{T} #in %s",kt.second.Data()));
            hRatKt[kt.first-1]->SetMarkerStyle(20);
            hRatKt[kt.first-1]->SetMarkerColor(JJColor::fWutAllColors[kt.first-1]);

            hRatKt[kt.first-1]->Write();
            if (kt.first -1 == 0)
                hRatKt[kt.first-1]->Draw("hist pe pmc plc");
            else
                hRatKt[kt.first-1]->Draw("same hist pe pmc plc");
        }
    }
    
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvKt->Write();

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    std::vector<TH1D*> hRatY(yArr.size(),nullptr), hBckgY(yArr.size(),nullptr);
    for (const auto &y : yArr)
    {
        for (const auto &kt : ktArr)
            for (const auto &psi : psiArr)
            {
                if (hRatY[y.first-1] == nullptr && hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckgY[y.first-1] == nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatY[y.first-1] = new TH1D(*hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgY[y.first-1] = new TH1D(*hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
                else if(hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatY[y.first-1]->Add(hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgY[y.first-1]->Add(hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
            }

        if (hRatY[y.first-1] != nullptr && hBckgY[y.first-1] != nullptr)
        {
            hRatY[y.first-1]->Divide(hBckgY[y.first-1]);
            norm = getNorm(hRatY[y.first-1],300,900);
            hRatY[y.first-1]->Rebin(rebin);
            norm *= rebin;
            hRatY[y.first-1]->Scale(1./norm);
            hRatY[y.first-1]->SetName(TString::Format("hQinvRatY%d",y.first));
            hRatY[y.first-1]->SetTitle(TString::Format("y_{c.m.} #in %s",y.second.Data()));
            hRatY[y.first-1]->SetMarkerStyle(20);

            hRatY[y.first-1]->Write();
            if (y.first -1 == 0)
                hRatY[y.first-1]->Draw("hist pe pmc plc");
            else
                hRatY[y.first-1]->Draw("same hist pe pmc plc");
        }
    }
    
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvY->Write();

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    std::vector<TH1D*> hRatPsi(psiArr.size(),nullptr), hBckgPsi(psiArr.size(),nullptr);
    for (const auto &psi : psiArr)
    {
        for (const auto &kt : ktArr)
            for (const auto &y : yArr)
            {
                if (hRatPsi[psi.first-1] == nullptr && hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckgPsi[psi.first-1] == nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatPsi[psi.first-1] = new TH1D(*hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgPsi[psi.first-1] = new TH1D(*hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
                else if(hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hRatPsi[psi.first-1]->Add(hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgPsi[psi.first-1]->Add(hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
            }

        if (hRatPsi[psi.first-1] != nullptr && hBckgPsi[psi.first-1] != nullptr)
        {
            hRatPsi[psi.first-1]->Divide(hBckgPsi[psi.first-1]);
            norm = getNorm(hRatPsi[psi.first-1],100,300);
            hRatPsi[psi.first-1]->Rebin(rebin);
            norm *= rebin;
            hRatPsi[psi.first-1]->Scale(1./norm);
            hRatPsi[psi.first-1]->SetName(TString::Format("hQinvRatPsi%d",psi.first));
            hRatPsi[psi.first-1]->SetTitle(TString::Format("#phi - #Psi_{E.P.} #in %s",psi.second.Data()));
            hRatPsi[psi.first-1]->SetMarkerStyle(20);

            hRatPsi[psi.first-1]->Write();
            if (psi.first -1 == 0)
                hRatPsi[psi.first-1]->Draw("hist pe pmc plc");
            else
                hRatPsi[psi.first-1]->Draw("same hist pe pmc plc");
        }
    }
    
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();
}