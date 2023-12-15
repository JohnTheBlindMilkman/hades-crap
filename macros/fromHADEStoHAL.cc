#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

TH1D* ConvertXaxisUnits(TH1D *hInp)
{
    const double MeVtoGeV = 1./1000.;
    const int nBins = hInp->GetNbinsX();
    const double newMin = hInp->GetXaxis()->GetXmin() * MeVtoGeV;
    const double newMax = hInp->GetXaxis()->GetXmax() * MeVtoGeV; // /2 to convert from qinv to k*

    TH1D *hOtp = new TH1D(hInp->GetName(),hInp->GetTitle(),nBins,newMin,newMax);
    //hOtp->GetXaxis()->SetTitle("k* [MeV/c]");

    for (int i = 1; i <= nBins; ++i)
    {
        hOtp->SetBinContent(i,hInp->GetBinContent(i));
        hOtp->SetBinError(i,hInp->GetBinError(i));
    }

    return hOtp;
}

void fromHADEStoHAL()
{
    const TString inpFilePath = "/home/jedkol/Downloads/HADES/HADES-CrAP/output/1Dcorr_0_10_cent.root";
    const TString inpHistBase = "hQinvRat";
    const std::vector<int> ktBins = {0,1,2,3,4};

    std::vector<TH1D*> histArray;
    TFile *inpFile,*otpFile;

    inpFile = TFile::Open(inpFilePath);
    for (const int &ktval : ktBins)
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("%s%d",inpHistBase.Data(),ktval))));

    // create output file which has "_forHAL" added to its name
    TString otpFilePath = inpFilePath;
    otpFilePath.Insert(otpFilePath.First('.'),"_forHAL");

    otpFile = TFile::Open(otpFilePath,"RECREATE");

    for (const auto elem : histArray)
        elem->Write();
}