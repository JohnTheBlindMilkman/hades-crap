#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"
#include "TLine.h"

void drawProton1DInteg()
{
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12sim_all_25_07_22_processed.root";
    const TString outputFile = "../output/1Dcorr_0_10_cent_Integ.root";
    const int rebin = 1;
    constexpr int minX = 0, maxX = 300;

    float norm;
    TH1D *hRat, *hSign, *hBckg;

    TLine *line = new TLine(minX,1,maxX,1);
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
    norm = JJUtils::CF::GetNormByRange(hRat,200,600);
    hRat->Rebin(rebin);
    norm *= rebin;
    hRat->Scale(1./norm);
    hRat->SetName("hQinvRatInteg");
    hRat->SetTitle(";q_{inv} [MeV/c];CF(q_{inv})");
    hRat->SetMarkerStyle(20);
    hRat->SetMarkerColor(TColor::GetColor(JJColor::navyWut.AsHexString()));

    hRat->GetYaxis()->SetTitleSize(0.06);
    hRat->GetYaxis()->SetLabelSize(0.06);
    hRat->GetYaxis()->SetNdivisions(506);

    hRat->GetXaxis()->SetTitleSize(0.06);
    hRat->GetXaxis()->SetLabelSize(0.06);
    hRat->GetXaxis()->SetNdivisions(506);

    hRat->GetXaxis()->SetRangeUser(minX,maxX);
    hRat->GetYaxis()->SetRangeUser(0,1.9);

    TFile *outFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canv = new TCanvas("canv","",1600,900);
    canv->SetMargin(0.15,0.05,0.15,0.05);
    hRat->Write();
    hSign->Write();
    hBckg->Write();

    hRat->Draw("hist pe pmc plc");
    line->Draw("same");
    //canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canv->Write();
}
