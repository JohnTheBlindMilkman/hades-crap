#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TString.h"
#include "../FemtoMixer/PairUtils.hxx"

TH1D* ConvertXaxisUnits(TH1D *hInp)
{
    if (hInp == nullptr)
        return nullptr;

    constexpr double GeVtoMeV = 1000.;
    constexpr double MeVtoGeV = 1. / GeVtoMeV;
    const int nBins = hInp->GetNbinsX();
    const double newMin = hInp->GetXaxis()->GetXmin() /* * MeVtoGeV */;
    const double newMax = hInp->GetXaxis()->GetXmax() / 2. /* * MeVtoGeV */; // /2 to convert from qinv to k*

    TH1D *hOtp = new TH1D(hInp->GetName(),hInp->GetTitle(),nBins,newMin,newMax);
    hOtp->GetXaxis()->SetTitle("k* [MeV/c]");

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
    const TString inpFilePath = "../output/1Dcorr_30_40_cent_Purity_MomRes.root";
    const auto ktArr = Mixing::PairGrouping{}.GetKtIndexSequence();
    const auto yArr = Mixing::PairGrouping{}.GetRapIndexSequence();
    //const std::vector<TString> sProjName{"out","side","long"};

    std::vector<TH1D*> histArray;
    std::vector<TH1D*> histProjArray;
    TFile *inpFile,*otpFile;

    inpFile = TFile::Open(inpFilePath);
    for (const auto &ktval : ktArr)
    {
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvRatKt%ld",ktval))));
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvSignKt%ld",ktval))));
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%ld",ktval))));
    }
    for (const auto &rapval : yArr)
    {
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvRatY%ld",rapval))));
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvSignY%ld",rapval))));
        histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvBckgY%ld",rapval))));
    }
    for (const auto &ktval : ktArr)
        for (const auto &rapval : yArr)
        {
            histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvRatKt%ldY%ld",ktval,rapval))));
            histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvSignKt%ldY%ld",ktval,rapval))));
            histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQinvBckgKt%ldY%ld",ktval,rapval))));
        }
    /* for (const auto &proj : sProjName)
    {
        for (const int &ktval : ktArr)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatKt%d",proj.Data(),ktval))));
        for (const int &rapval : yArr)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatY%d",proj.Data(),rapval))));
        for (const int &psival : psiBins)
            histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>(TString::Format("hQ%sRatPsi%d",proj.Data(),psival))));
    } */

    histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQinvRatInteg")));
    histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQinvSignInteg")));
    histArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQinvBckgInteg")));
    //histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQoutRatInteg")));
    //histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQsideRatInteg")));
    //histProjArray.push_back(ConvertXaxisUnits(inpFile->Get<TH1D>("hQlongRatInteg")));

    // create output file which has "_forHAL" added to its name
    TString otpFilePath = inpFilePath;
    otpFilePath.Insert(otpFilePath.Last('.'),"_forHAL");

    otpFile = TFile::Open(otpFilePath,"RECREATE");

    for (const auto elem : histArray)
        if (elem != nullptr)
            elem->Write();

    /* for (const auto elem : histProjArray)
        elem->Write(); */
}