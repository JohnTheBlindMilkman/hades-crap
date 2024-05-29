#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
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

void prepareGraph(TH1D* hist, int col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    //hist->GetXaxis()->SetLabelSize();
    //hist->GetXaxis()->SetLabelOffset();
    //hist->GetXaxis()->SetTitleSize();
    //hist->GetXaxis()->SetTitleOffset();

    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelOffset(0.005);
    hist->GetYaxis()->SetLabelFont(62);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleFont(62);
    hist->GetYaxis()->SetNdivisions(510);
}

void drawProton3DMultiDiff()
{
    gStyle->SetOptStat(0);

    const TString fileName = "/u/kjedrzej/hades-crap/slurmOutput/apr12ana_all_24_05_29_processed.root";
    const TString outputFile = "/u/kjedrzej/hades-crap/output/3Dcorr_0_10_cent.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(150,450)"},{2,"(450,750)"},{3,"(750,1050)"},{4,"(1050,1350)"},{5,"(1350,1650)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.75,-0.25)"},{2,"(-0.25,0.25)"},{3,"(0.25,0.75)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};
    const int rebin = 1;
    const int wbin = 1; 

    float norm;
    int binc, binmn, binmx;

    TLine *line = new TLine(-250,1,250,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1800,600);
    canvKt->Divide(3,1);
    for (const auto &kt : ktArr)
    {
        TH3D *hSign = inpFile->Get<TH3D>(TString::Format("hQoslSignKt%d",kt.first));
        TH3D *hBckg = inpFile->Get<TH3D>(TString::Format("hQoslBckgKt%d",kt.first));
        TH3D *hRat3D = new TH3D(*hSign);
        hRat3D->Divide(hBckg);
        hRat3D->SetName(TString::Format("hQoslRatKt%d",kt.first));
        hRat3D->Write();

        if (hRat3D != nullptr)
        {
            for (const int &i : {0,1,2})
            {
                binc = hRat3D->GetXaxis()->FindFixBin(0.0);
                binmx = binc + wbin;
                binmn = binc - wbin;

                hRat3D->GetXaxis()->SetRange((i == 0) ? 1 : binmn, (i == 0) ? hSign->GetNbinsX() : binmx);
                hRat3D->GetYaxis()->SetRange((i == 1) ? 1 : binmn, (i == 1) ? hSign->GetNbinsY() : binmx);
                hRat3D->GetZaxis()->SetRange((i == 2) ? 1 : binmn, (i == 2) ? hSign->GetNbinsZ() : binmx);

                TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
                norm = getNorm(hRat,150,250);
                hRat->Rebin(rebin);
                norm *= rebin;
                hRat->Scale(1./norm);
                hRat->GetYaxis()->SetRangeUser(0.,2.);
                hRat->SetName(TString::Format("hQ%sRatKt%d",sProjName[i].Data(),kt.first));
                hRat->SetTitle(TString::Format("k_{T} #in %s;q_{%s} [MeV/c];CF(q_{%s})",kt.second.Data(),sProjName[i].Data(),sProjName[i].Data()));
                hRat->SetMarkerStyle(20);
                hRat->SetMarkerColor(JJColor::fWutAllColors[kt.first-1]);

                hRat->Write();

                canvKt->cd(i+1);
                if (kt.first -1 == 0)
                    hRat->Draw("hist pe pmc plc");
                else
                    hRat->Draw("same hist pe pmc plc");
            }
        }
    }
    
    for (const int &i : {0,1,2})
    {
        TVirtualPad *tvp = canvKt->cd(i+1);
        if (i == 1)
            tvp->BuildLegend(.3,.21,.3,.21,"","p");
        line->Draw("same");
    } 
    canvKt->Write();

    TCanvas *canvY = new TCanvas("canvY","",1800,600);
    canvY->Divide(3,1);
    for (const auto &y : yArr)
    {
        TH3D *hSign = inpFile->Get<TH3D>(TString::Format("hQoslSignY%d",y.first));
        TH3D *hBckg = inpFile->Get<TH3D>(TString::Format("hQoslBckgY%d",y.first));
        TH3D *hRat3D = new TH3D(*hSign);
        hRat3D->Divide(hBckg);
        hRat3D->SetName(TString::Format("hQoslRatY%d",y.first));
        hRat3D->Write();

        if (hRat3D != nullptr)
        {
            for (const int &i : {0,1,2})
            {
                binc = hRat3D->GetXaxis()->FindFixBin(0.0);
                binmx = binc + wbin;
                binmn = binc - wbin;

                hRat3D->GetXaxis()->SetRange((i == 0) ? 1 : binmn, (i == 0) ? hSign->GetNbinsX() : binmx);
                hRat3D->GetYaxis()->SetRange((i == 1) ? 1 : binmn, (i == 1) ? hSign->GetNbinsY() : binmx);
                hRat3D->GetZaxis()->SetRange((i == 2) ? 1 : binmn, (i == 2) ? hSign->GetNbinsZ() : binmx);

                TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
                norm = getNorm(hRat,150,250);
                hRat->Rebin(rebin);
                norm *= rebin;
                hRat->Scale(1./norm);
                hRat->GetYaxis()->SetRangeUser(0.,2.);
                hRat->SetName(TString::Format("hQ%sRatY%d",sProjName[i].Data(),y.first));
                hRat->SetTitle(TString::Format("y #in %s;q_{%s} [MeV/c];CF(q_{%s})",y.second.Data(),sProjName[i].Data(),sProjName[i].Data()));
                hRat->SetMarkerStyle(20);
                hRat->SetMarkerColor(JJColor::fWutAllColors[y.first-1]);

                hRat->Write();

                canvY->cd(i+1);
                if (y.first -1 == 0)
                    hRat->Draw("hist pe pmc plc");
                else
                    hRat->Draw("same hist pe pmc plc");
            }
        }
    }
    
    for (const int &i : {0,1,2})
    {
        TVirtualPad *tvp = canvY->cd(i+1);
        if (i == 1)
            tvp->BuildLegend(.3,.21,.3,.21,"","p");
        line->Draw("same");
    } 
    canvY->Write();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1800,600);
    canvPsi->Divide(3,1);
    for (const auto &psi : psiArr)
    {
        TH3D *hSign = inpFile->Get<TH3D>(TString::Format("hQoslSignPsi%d",psi.first));
        TH3D *hBckg = inpFile->Get<TH3D>(TString::Format("hQoslBckgPsi%d",psi.first));
        TH3D *hRat3D = new TH3D(*hSign);
        hRat3D->Divide(hBckg);
        hRat3D->SetName(TString::Format("hQoslRatPsi%d",psi.first));
        hRat3D->Write();

        if (hRat3D != nullptr)
        {
            for (const int &i : {0,1,2})
            {
                binc = hRat3D->GetXaxis()->FindFixBin(0.0);
                binmx = binc + wbin;
                binmn = binc - wbin;

                hRat3D->GetXaxis()->SetRange((i == 0) ? 1 : binmn, (i == 0) ? hSign->GetNbinsX() : binmx);
                hRat3D->GetYaxis()->SetRange((i == 1) ? 1 : binmn, (i == 1) ? hSign->GetNbinsY() : binmx);
                hRat3D->GetZaxis()->SetRange((i == 2) ? 1 : binmn, (i == 2) ? hSign->GetNbinsZ() : binmx);

                TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
                norm = getNorm(hRat,150,250);
                hRat->Rebin(rebin);
                norm *= rebin;
                hRat->Scale(1./norm);
                hRat->GetYaxis()->SetRangeUser(0.,2.);
                hRat->SetName(TString::Format("hQ%sRatPsi%d",sProjName[i].Data(),psi.first));
                hRat->SetTitle(TString::Format("#phi - #Psi_{E.P.} #in %s;q_{%s} [MeV/c];CF(q_{%s})",psi.second.Data(),sProjName[i].Data(),sProjName[i].Data()));
                hRat->SetMarkerStyle(20);
                hRat->SetMarkerColor(JJColor::fWutAllColors[psi.first-1]);

                hRat->Write();

                canvPsi->cd(i+1);
                if (psi.first -1 == 0)
                    hRat->Draw("hist pe pmc plc");
                else
                    hRat->Draw("same hist pe pmc plc");
            }
        }
    }
    
    for (const int &i : {0,1,2})
    {
        TVirtualPad *tvp = canvPsi->cd(i+1);
        if (i == 1)
            tvp->BuildLegend(.3,.21,.3,.21,"","p");
        line->Draw("same");
    } 
    canvPsi->Write();
}