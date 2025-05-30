#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

void makeMomResSigma()
{
    const TString fileName = "../slurmOutput/apr12sim_all_25_05_16.root";
    const TString inpHistName = "hQinvResolution";
    const TString outputFile = "../output/1DMomResSigma_0_40_cent.root";

    TFile *inpFile = TFile::Open(fileName);
    TH2D *hQinvResolution = inpFile->Get<TH2D>(inpHistName);
    TH1D *hSigma = new TH1D(
        "hQinvSigma",
        "Distribution of #sigma^{reco}; q_{inv}^{kine} [MeV/c]; #sigma^{reco} [MeV/c]",
        hQinvResolution->GetNbinsY(),
        hQinvResolution->GetYaxis()->GetXmin(),
        hQinvResolution->GetYaxis()->GetXmax()
    );

    TH1D *hQinvProj = nullptr;
    for (int bin = 1; bin <= hQinvResolution->GetNbinsX(); ++bin)
    {
        hQinvProj = hQinvResolution->ProjectionY("_py",bin,bin);
        if (hQinvProj != nullptr)
        {
            hSigma->SetBinContent(bin,hQinvProj->GetStdDev());
            hSigma->SetBinError(bin,hQinvProj->GetStdDevError());
            delete hQinvProj;
        }
    }
    TFile *otpFile = TFile::Open(outputFile,"recreate");
    hSigma->Write();
}