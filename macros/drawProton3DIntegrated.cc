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

void drawProton3DIntegrated()
{
    gStyle->SetOptStat(0);
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12ana_quarter_24_06_05_processed.root";
    const TString outputFile = "../output/3Dcorr_0_10_cent_Integ.root";
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};
    const int rebin = 1;
    const int wbin = 1; 

    float norm;
    int binc, binmn, binmx;

    TLine *line = new TLine(-500,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canvInteg = new TCanvas("canvInteg","",1800,600);
    canvInteg->Divide(3,1);

    TH3D *hSign = inpFile->Get<TH3D>("hQoslSignInteg");
    TH3D *hBckg = inpFile->Get<TH3D>("hQoslBckgInteg");
    TH3D *hRat3D = new TH3D(*hSign);
    hRat3D->Divide(hBckg);
    hRat3D->SetName("hQoslRatInteg");
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
            norm = getNorm(hRat,150,500);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->GetYaxis()->SetRangeUser(0.,1.9);
            hRat->SetName(TString::Format("hQ%sRatInteg",sProjName[i].Data()));
            hRat->SetTitle(TString::Format(";q_{%s} [MeV/c];CF(q_{%s})",sProjName[i].Data(),sProjName[i].Data()));
            hRat->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
            hRat->GetXaxis()->SetTitleSize(0.06);
            hRat->GetXaxis()->SetLabelSize(0.06);
            hRat->GetXaxis()->SetNdivisions(506);
            hRat->GetYaxis()->SetTitleSize(0.06);
            hRat->GetYaxis()->SetLabelSize(0.06);
            hRat->GetYaxis()->SetNdivisions(506);
            hRat->SetMarkerStyle(20);
            hRat->SetMarkerColor(JJColor::fWutAllColors[1]); // navy WUT

            hRat->Write();

            canvInteg->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);
            hRat->Draw("hist pe pmc plc");
            line->Draw("same");
        }

        canvInteg->Write();
    }
}