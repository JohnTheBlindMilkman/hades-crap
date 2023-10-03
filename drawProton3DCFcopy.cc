#include <iostream>
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH3.h"
#include "TLine.h"

double getNorm(TH1D *hInp, double xMin, double xMax)
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
        return (double)(val/nBins);
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

void drawProton3DCFcopy()
{
    gStyle->SetOptStat(0);
    const TString filePath = "testOutFile.root";
    const TString sProj[] = {"x","y","z"};
    const TString sProjName[] = {"out","side","long"};
    const int cent = 1;
    const double drawMin = 0, drawMax = 1000;
    const int rebin = 2;
    const int wbin = 20;

    double norm = 1.;
    int binc = 0, binmn = 0, binmx = 0;
    TFile *impFile,*outFile;
    TCanvas *canv;
    TH3D *hNum3D, *hDen3D;
    TH1D *hNumPP[3],*hDenPP[3],*hRatPP[3];
    TLine *lline;

    canv = new TCanvas("canv","",1800,600);
    canv->Divide(3,1);

    lline = new TLine(drawMin,1.,drawMax,1.);
    lline->SetLineStyle(kDashed);

    impFile = TFile::Open(filePath);

    hNum3D = (TH3D*) impFile->Get("hQoslSign");
    hDen3D = (TH3D*) impFile->Get("hQoslBckg");

    for(int i : {0,1,2}) // loop over projections: out, side, long
    {
        binc = hNum3D->GetXaxis()->FindFixBin(0.0);
        binmx = binc + wbin;
        binmn = binc - wbin;

        hNum3D->GetXaxis()->SetRange(1, (i == 0) ? hNum3D->GetNbinsX() : binmx);
        hNum3D->GetYaxis()->SetRange(1, (i == 1) ? hNum3D->GetNbinsY() : binmx);
        hNum3D->GetZaxis()->SetRange(1, (i == 2) ? hNum3D->GetNbinsZ() : binmx);
        hDen3D->GetXaxis()->SetRange(1, (i == 0) ? hDen3D->GetNbinsX() : binmx);
        hDen3D->GetYaxis()->SetRange(1, (i == 1) ? hDen3D->GetNbinsY() : binmx);
        hDen3D->GetZaxis()->SetRange(1, (i == 2) ? hDen3D->GetNbinsZ() : binmx);

        hNumPP[i] = (TH1D*) hNum3D->Project3D(sProj[i].Data());
        hDenPP[i] = (TH1D*) hDen3D->Project3D(sProj[i].Data());

        hRatPP[i] = new TH1D(*hNumPP[i]);
        hRatPP[i]->SetName(Form("hQ%s",sProjName[i].Data()));
        hRatPP[i]->Divide(hDenPP[i]);
        //hRatPP[i]->Sumw2();
        //hRatPP[i]->Reset("ICE");
        hRatPP[i]->Rebin(rebin);
        norm = getNorm(hRatPP[i],150,300);
        hRatPP[i]->GetXaxis()->SetRangeUser(drawMin,drawMax);
        //hRatPP[i]->GetYaxis()->SetRangeUser(0.,2.2);
        hRatPP[i]->SetTitle(Form("Projection of 3D p-p correlation;q_{%s} [MeV/c];CF(q_{%s})",sProjName[i].Data(),sProjName[i].Data()));
        prepareGraph(hRatPP[i],TColor::GetColor("#965f77"));
        if (norm > 0)
            hRatPP[i]->Scale(1./norm);
        canv->cd(i+1)->SetMargin(0.15,0.01,0.15,0.1);
        hRatPP[i]->Draw("p histo");
        lline->Draw("same");
    }
    
    canv->SaveAs("./output/CF3D_cent0_10_kTInt.png");
    //outFile = TFile::Open("./output/ProtonProton1BPLCMS3DCF.root","recreate");
    //canv->Write();
    //for( int i : {0,1,2})
        //hRatPP[i]->Write();

    //outFile->Close();
    //impFile->Close();
}