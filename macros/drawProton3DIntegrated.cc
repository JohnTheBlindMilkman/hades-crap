#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"
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

double getNorm(const TH2D *hInp, double xMin, double xMax, double yMin, double yMax)
{
    const int binsX = hInp->GetNbinsX();
    const int binxY = hInp->GetNbinsY();
    int div = 0;
    double sum = 0., xVal, yVal;
    for (int i = 1; i <= binsX; ++i)
        for (int j = 1; j < binxY; ++j)
        {
            xVal = hInp->GetXaxis()->GetBinCenter(i);
            yVal = hInp->GetYaxis()->GetBinCenter(j);
            if (xVal >= xMin && xVal <= xMax && yVal >= yMin && yVal <= yMax)
            {
                sum += hInp->GetBinContent(i,j);
                ++div;
            } 
        }
    if(div > 0)
        return sum/div;
    else    
        return 0.;
}

double getNorm(const TH3D *num, const TH3D *den)
{
    return den->GetEntries()/num->GetEntries();
}

void drawProton3DIntegrated()
{
    gStyle->SetOptStat(0);
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12ana_all_24_09_13_processed.root";
    const TString outputFile = "../output/3Dcorr_0_10_cent_Integ_tmp3.root";
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};
    const int rebin = 2;
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

    TPaveText *ptInfo = new TPaveText(.35,.7,.8,.95,"NB");
    ptInfo->SetBorderSize(0);
    ptInfo->SetTextAlign(12);
    ptInfo->SetTextFont(102);
    ptInfo->SetFillStyle(1); // transparemt background
    ptInfo->AddText("proton - proton");
    ptInfo->AddText("0-10%");
    ptInfo->AddText("y_{c.m.} #in (-0.75,0.75)");
    ptInfo->AddText("k_{T} #in (150,1650) [MeV/c]");

    TH3D *hSign = inpFile->Get<TH3D>("hQoslSignInteg");
    TH3D *hBckg = inpFile->Get<TH3D>("hQoslBckgInteg");
    TH3D *hRat3D = new TH3D(*hSign);
    hRat3D->Divide(hBckg);
    hRat3D->SetName("hQoslRatInteg");
    //norm = getNorm(hSign,hBckg);
    //hRat3D->Scale(norm);
    hRat3D->Write();

    if (hRat3D != nullptr)
    {
        for (const int &i : {0,1,2})
        {
            binc = hRat3D->GetXaxis()->FindFixBin(0.0);
            binmx = binc + wbin;
            //binmn = binc - wbin ; // I don't have a central bin at 0...

            hRat3D->GetXaxis()->SetRange(1, (i == 0) ? hSign->GetNbinsX() : binmx);
            hRat3D->GetYaxis()->SetRange(1, (i == 1) ? hSign->GetNbinsY() : binmx);
            hRat3D->GetZaxis()->SetRange(1, (i == 2) ? hSign->GetNbinsZ() : binmx);

            TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
            norm = getNorm(hRat,350,500);
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

            TPaveText *pt = new TPaveText(0.3,0.4,0.7,0.6,"NB");
            pt->SetBorderSize(0);
            pt->SetFillStyle(1);
            pt->SetTextAlign(22);
            pt->SetTextFont(102);
            pt->AddText(TString::Format("q_{%s} #in (0,%.2f) [MeV/c]",sProjName[(i+1)%3].Data(),hRat->GetXaxis()->GetBinUpEdge(binmx)));
            pt->AddText(TString::Format("q_{%s} #in (0,%.2f) [MeV/c]",sProjName[(i+2)%3].Data(),hRat->GetXaxis()->GetBinUpEdge(binmx)));

            canvInteg->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);
            hRat->Draw("hist pe pmc plc");
            line->Draw("same");
            pt->Draw("same");
            if (i == 0)
                ptInfo->Draw("same");
        }

        canvInteg->Write();

        TCanvas *canv2DInteg = new TCanvas("canv2DInteg","",1800,600);
        canv2DInteg->Divide(3,1);
        //JJColor::CreatePrimaryWutGradient();
        gStyle->SetPalette(kPastel);
        for (const int &i : {0,1,2})
        {
            binc = hRat3D->GetXaxis()->FindFixBin(0.0);
            binmx = binc + wbin;

            hRat3D->GetXaxis()->SetRange(1, (i == 0 || (i+1)%3 == 0) ? hRat3D->GetNbinsX() : binmx);
            hRat3D->GetYaxis()->SetRange(1, (i == 1 || (i+1)%3 == 1) ? hRat3D->GetNbinsY() : binmx);
            hRat3D->GetZaxis()->SetRange(1, (i == 2 || (i+1)%3 == 2) ? hRat3D->GetNbinsZ() : binmx);

            TH2D *hRat = static_cast<TH2D*>(hRat3D->Project3D((sProj[i%3]+sProj[(i+1)%3]).Data()));
            norm = getNorm(hRat,350,500,350,500);
            hRat->Rebin2D(rebin,rebin);
            norm *= rebin*rebin;
            hRat->Scale(1./norm);
            hRat->GetZaxis()->SetRangeUser(0.5,1.5);
            hRat->SetName(TString::Format("hQ%s%sRatInteg",sProjName[i%3].Data(),sProjName[(i+1)%3].Data()));
            hRat->SetTitle(TString::Format(";q_{%s} [MeV/c];q_{%s} [MeV/c];",sProjName[i].Data(),sProjName[(i+1)%3].Data()));
            hRat->GetXaxis()->SetTitleOffset();
            hRat->GetYaxis()->SetTitleOffset(1.8);
            hRat->GetXaxis()->SetTitleSize(0.06);
            hRat->GetXaxis()->SetLabelSize(0.06);
            hRat->GetXaxis()->SetNdivisions(506);
            hRat->GetYaxis()->SetTitleSize(0.06);
            hRat->GetYaxis()->SetLabelSize(0.06);
            hRat->GetYaxis()->SetNdivisions(506);
            //hRat->SetMarkerStyle(20);
            //hRat->SetMarkerColor(JJColor::fWutAllColors[1]); // navy WUT

            hRat->Write();

            canv2DInteg->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);
            hRat->Draw("colz");
        }

        canv2DInteg->Write();
    }
}
