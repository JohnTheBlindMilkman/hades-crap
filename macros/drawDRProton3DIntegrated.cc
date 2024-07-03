#include "TStyle.h"
#include "TLine.h"
#include "TF3.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "../Externals/Palettes.hxx"

double getNorm(const TH1D *hInp, double xMin, double xMax)
{
    int nBins = 0;
    double val = 0., xVal;
    for (int i = 1; i <= hInp->GetNbinsX(); i++)
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

void SetErrors(TH3D *hout, TH3D *hNum, TH3D *hDen)
{
    const int binsX = hout->GetNbinsX();
    const int binsY = hout->GetNbinsY();
    const int binsZ = hout->GetNbinsZ();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= binsX; ++i)
        for (int j = 1; j <= binsY; ++j)
            for (int k = 1; k <= binsZ; ++k)
            {
                vErr = 0;
                vNum = hNum->GetBinContent(i,j,k);
                eNum = hNum->GetBinError(i,j,k);
                vDen = hDen->GetBinContent(hDen->GetXaxis()->FindBin(hNum->GetXaxis()->GetBinCenter(i)),hDen->GetYaxis()->FindBin(hNum->GetYaxis()->GetBinCenter(j)),hDen->GetZaxis()->FindBin(hNum->GetZaxis()->GetBinCenter(k)));
                eDen = hDen->GetBinError(hDen->GetXaxis()->FindBin(hNum->GetXaxis()->GetBinCenter(i)),hDen->GetYaxis()->FindBin(hNum->GetYaxis()->GetBinCenter(j)),hDen->GetZaxis()->FindBin(hNum->GetZaxis()->GetBinCenter(k)));

                // propagation of uncertainty for a function num/den with inclusion of the correlation between the constituents
                if (fabs(vDen) > 0.)
                    vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
                hout->SetBinError(i,j,k,vErr);
            }
}

void AnalyticalDivide(TH3D *hNum, TF3 *fDen)
{
    const int binsX = hNum->GetNbinsX();
    const int binsY = hNum->GetNbinsY();
    const int binsZ = hNum->GetNbinsZ();
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    fDen->GetRange(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);

    for (int i = 1; i <= binsX; ++i)
        for (int j = 1; j <= binsY; ++j)
            for (int k = 1; k <= binsZ; ++k)
            {
                if ((hNum->GetXaxis()->GetBinCenter(i) < Xmax && hNum->GetXaxis()->GetBinCenter(i) > Xmin) ||
                (hNum->GetYaxis()->GetBinCenter(j) < Ymax && hNum->GetYaxis()->GetBinCenter(j) > Ymin) || 
                (hNum->GetZaxis()->GetBinCenter(k) < Zmax && hNum->GetZaxis()->GetBinCenter(k) > Zmin))
                    hNum->SetBinContent(i,j,k,hNum->GetBinContent(i,j,k) / fDen->Eval(hNum->GetXaxis()->GetBinCenter(i),hNum->GetYaxis()->GetBinCenter(j),hNum->GetZaxis()->GetBinCenter(k)));
                else
                    hNum->SetBinContent(i,j,k,0.);
            }
}

void drawDRProton3DIntegrated()
{
    gStyle->SetOptStat(0);
    const TString inpFile = "../output/3Dcorr_0_10_cent_Integ.root";
    const TString simFile = "../output/3Dcorr_0_10_cent_HGeant_Integ_fit.root";
    const TString otpFile = "../output/3Dcorr_0_10_cent_DR_Integ.root";
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};
    const int rebin = 1;
    const int wbin = 1;

    float binc,binmx,binmn,norm;
    TLine line(0,1,500,1);
    line.SetLineColor(kGray);
    line.SetLineStyle(kDashed);

    TFile *fInpData = TFile::Open(inpFile);
    TFile *fInpSim = TFile::Open(simFile);

    TH3D *hExp3D = fInpData->Get<TH3D>("hQoslRatInteg");
    TF3 *fFit3D = fInpSim->Get<TF3>("fitFunc");
    TH3D *hErr3D = fInpSim->Get<TH3D>("hQoslErrInteg");

    TH3D *hRat3D = new TH3D(*hExp3D);
    //hRat3D->Divide(fFit3D);
    AnalyticalDivide(hRat3D,fFit3D);
    SetErrors(hRat3D,hExp3D,hErr3D);

    TCanvas *canvInteg = new TCanvas("canvInteg","",1800,600);
    canvInteg->Divide(3,1);

    TFile *outputFile = TFile::Open(otpFile,"RECREATE");
    hRat3D->Write();

    for (const int &i : {0,1,2})
    {
        binc = hRat3D->GetXaxis()->FindFixBin(0.0);
        binmx = binc + wbin;
        binmn = binc - wbin;

        hRat3D->GetXaxis()->SetRange(binc, (i == 0) ? hExp3D->GetNbinsX() : binmx);
        hRat3D->GetYaxis()->SetRange(binc, (i == 1) ? hExp3D->GetNbinsY() : binmx);
        hRat3D->GetZaxis()->SetRange(binc, (i == 2) ? hExp3D->GetNbinsZ() : binmx);

        TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
        norm = getNorm(hRat,350,500);
        hRat->Rebin(rebin);
        norm *= rebin;
        hRat->Scale(1./norm);
        hRat->GetYaxis()->SetRangeUser(0.,1.9);
        hRat->GetXaxis()->SetRangeUser(0,500);
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
        line.Draw("same");
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
        binmn = binc - wbin;

        hRat3D->GetXaxis()->SetRange(binc, (i == 0 || (i+1)%3 == 0) ? hExp3D->GetNbinsX() : binmx);
        hRat3D->GetYaxis()->SetRange(binc, (i == 1 || (i+1)%3 == 1) ? hExp3D->GetNbinsY() : binmx);
        hRat3D->GetZaxis()->SetRange(binc, (i == 2 || (i+1)%3 == 2) ? hExp3D->GetNbinsZ() : binmx);

        TH2D *hRat = static_cast<TH2D*>(hRat3D->Project3D((sProj[i%3]+sProj[(i+1)%3]).Data()));
        norm = getNorm(hRat,350,500,350,500);
        hRat->Rebin2D(rebin,rebin);
        norm *= rebin;
        hRat->Scale(1./norm);
        hRat->GetYaxis()->SetRangeUser(0,500);
        hRat->GetXaxis()->SetRangeUser(0,500);
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

        canv2DInteg->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);;
        hRat->Draw("colz");
    }

    canv2DInteg->Write();
}