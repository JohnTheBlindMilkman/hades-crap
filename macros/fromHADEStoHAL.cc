#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TString.h"

TH1D* ConvertXaxisUnits(TH1D *hInp)
{
    const double MeVtoGeV = 1./1000.;
    const int nBins = hInp->GetNbinsX();
    const double newMin = hInp->GetXaxis()->GetXmin() /* * MeVtoGeV */;
    const double newMax = hInp->GetXaxis()->GetXmax() / 2. /* * MeVtoGeV */; // /2 to convert from qinv to k*

    TH1D *hOtp = new TH1D(hInp->GetName(),hInp->GetTitle(),nBins,newMin,newMax);
    //hOtp->GetXaxis()->SetTitle("k* [MeV/c]");

    for (int i = 1; i <= nBins; ++i)
    {
        hOtp->SetBinContent(i,hInp->GetBinContent(i));
        hOtp->SetBinError(i,hInp->GetBinError(i));
    }

    return hOtp;
}

TH3D* ConvertXaxisUnits(TH3D *hInp)
{
    const double MeVtoGeV = 1./1000.;
    const int nBinsX = hInp->GetNbinsX();
    const int nBinsY = hInp->GetNbinsY();
    const int nBinsZ = hInp->GetNbinsZ();
    const double newMinX = hInp->GetXaxis()->GetXmin();
    const double newMaxX = hInp->GetXaxis()->GetXmax() / 2.; // /2 to convert from qinv to k*
    const double newMinY = hInp->GetYaxis()->GetXmin();
    const double newMaxY = hInp->GetYaxis()->GetXmax() / 2.; // /2 to convert from qinv to k*
    const double newMinZ = hInp->GetZaxis()->GetXmin();
    const double newMaxZ = hInp->GetZaxis()->GetXmax() / 2.; // /2 to convert from qinv to k*

    TH3D *hOtp = new TH3D(hInp->GetName(),hInp->GetTitle(),nBinsX,newMinX,newMaxX,nBinsY,newMinY,newMaxY,nBinsZ,newMinZ,newMaxZ);
    //hOtp->GetXaxis()->SetTitle("k* [MeV/c]");

    for (int i = 1; i <= nBinsX; ++i)
        for (int j = 1; j <= nBinsY; ++j)
            for (int k = 1; k <= nBinsZ; ++k)
            {
                hOtp->SetBinContent(i,j,k,hInp->GetBinContent(i,j,k));
                hOtp->SetBinError(i,j,k,hInp->GetBinError(i,j,k));
            }

    return hOtp;
}

void fromHADEStoHAL()
{
    const TString inpFilePath = "../output/3Dcorr_0_10_cent.root";
    const TString inpFilePathInteg = "../output/3Dcorr_0_10_cent_Integ.root";
    const TString inpHistBaseKt = "hQoslRatKt";
    const TString inpHistBaseRap = "hQoslRatY";
    const TString inpHistBasePsi = "hQoslRatPsi";
    const TString inpHitsInteg = "hQoslRatInteg";
    const std::vector<int> ktBins = {1,2,3,4,5};
    const std::vector<int> rapBins = {1,2,3};
    const std::vector<int> psiBins = {1,2,3,4,5,6,7,8};
    const std::vector<TString> sProjName{"out","side","long"};

    std::vector<TH3D*> histArray;
    std::vector<TH1D*> histProjArray;
    TFile *inpFile,*otpFile;

    inpFile = TFile::Open(inpFilePath);
    for (const int &ktval : ktBins)
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH3D>(TString::Format("%s%d",inpHistBaseKt.Data(),ktval))));
    for (const int &rapval : rapBins)
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH3D>(TString::Format("%s%d",inpHistBaseRap.Data(),rapval))));
    for (const int &psival : psiBins)
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH3D>(TString::Format("%s%d",inpHistBasePsi.Data(),psival))));
    for (const auto &proj : sProjName)
    {
        for (const int &ktval : ktBins)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatKt%d",proj.Data(),ktval))));
        for (const int &rapval : rapBins)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatY%d",proj.Data(),rapval))));
        for (const int &psival : psiBins)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatPsi%d",proj.Data(),psival))));
    }

    inpFile = TFile::Open(inpFilePathInteg);
    histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH3D>(inpHitsInteg)));
    histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQoutRatInteg")));
    histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQsideRatInteg")));
    histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQlongRatInteg")));

    // create output file which has "_forHAL" added to its name
    TString otpFilePath = inpFilePath;
    otpFilePath.Insert(otpFilePath.Last('.'),"_forHAL");

    otpFile = TFile::Open(otpFilePath,"RECREATE");

    for (const auto elem : histArray)
        elem->Write();

    for (const auto elem : histProjArray)
        elem->Write();
}