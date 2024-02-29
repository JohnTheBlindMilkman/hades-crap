#include "TStyle.h"
#include "TLine.h"
#include "TF3.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "Palettes.hxx"

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
                vDen = hDen->GetBinContent(hDen->GetXaxis()->FindBin(hNum->GetBinCenter(i)),hDen->GetYaxis()->FindBin(hNum->GetBinCenter(j)),hDen->GetZaxis()->FindBin(hNum->GetBinCenter(k)));
                eDen = hDen->GetBinError(hDen->GetXaxis()->FindBin(hNum->GetBinCenter(i)),hDen->GetYaxis()->FindBin(hNum->GetBinCenter(j)),hDen->GetZaxis()->FindBin(hNum->GetBinCenter(k)));

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
    const TString inpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent.root";
    const TString simFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_HGeant_fit.root";
    const TString otpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR.root";
    const int rebin = 1;
    const int wbin = 1;

    TLine line(-250,1,250,1);
    line.SetLineColor(kGray);
    line.SetLineStyle(kDashed);

    TFile *fInpData = TFile::Open(inpFile);
    TFile *fInpSim = TFile::Open(simFile);
}