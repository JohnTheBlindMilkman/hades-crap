#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "Palettes.hxx"
#include "TLine.h"

#include "hadesifyPlot.cc"

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
        return val/nBins;
    else    
        return 0.;
}

void drawProton1DJJFM()
{
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12ana_all_24_02_19_processed.root";
    const TString outputFile = "../output/1Dcorr_0_10_cent_Integ.root";
    const int rebin = 2;

    float norm;
    TH1D *hRat, *hSign, *hBckg;

    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    hSign = inpFile->Get<TH1D>("hQinvSignInteg");
    hBckg = inpFile->Get<TH1D>("hQinvBckgInteg");
    
    hSign->Sumw2();
    hBckg->Sumw2();
    hRat = new TH1D(*hSign);

    hRat->Divide(hBckg);
    hRat->Sumw2();
    norm = getNorm(hRat,300,900);
    hRat->Rebin(rebin);
    norm *= rebin;
    hRat->Scale(1./norm);
    hRat->SetName("hQinvRat");
    hRat->SetTitle(";q_{inv} [MeV/c];CF(q_{inv})");
    hRat->SetMarkerStyle(20);
    hRat->SetMarkerColor(TColor::GetColor(JJColor::navyWut.AsHexString()));

    hRat->GetYaxis()->SetTitleSize(0.06);
    hRat->GetYaxis()->SetLabelSize(0.06);
    hRat->GetYaxis()->SetNdivisions(506);

    hRat->GetXaxis()->SetTitleSize(0.06);
    hRat->GetXaxis()->SetLabelSize(0.06);
    hRat->GetXaxis()->SetNdivisions(506);

    TFile *outFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canv = new TCanvas("canv","",1600,900);
    canv->SetMargin(0.15,0.01,0.15,0.01);
    hRat->Write();
    hSign->Write();
    hBckg->Write();

    hRat->Draw("hist pe pmc plc");
    line->Draw("same");
    //canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canv->Write();
}