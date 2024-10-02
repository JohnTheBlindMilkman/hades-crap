#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphErrors.h"

void fitMomentumResolution()
{
    const TString fileName = "../slurmOutput/apr12momres_all_24_10_02.root";
    const TString outputFile = "../output/momentum_resolution.root";
    const TString inpMomHist = "hMomResolution";
    const TString inpPhiHist = "hPhiResolution";
    const TString inpThetaHist = "hThetaResolution";
    constexpr int entriesCutoff{10000};
    constexpr bool fitSlices{true};
    constexpr bool fitGaussianParams{true};

    TFile *inpFile = TFile::Open(fileName);

    TH2D *hMomHist = inpFile->Get<TH2D>(inpMomHist);
    TH2D *hPsiHist = inpFile->Get<TH2D>(inpPhiHist);
    TH2D *hThetaHist = inpFile->Get<TH2D>(inpThetaHist);

    TF1 *fGauss = new TF1("fGauss","gaus",-500,500,3);
    fGauss->SetParameters(1,0,1);

    // I will now assume that all three histograms have the same X axis range and binning!
    // why? because I'm lazy

    const int nBins = hMomHist->GetNbinsX();
    std::vector<float> momMuCollection, momSigmaCollection, momMuErrCollection, momSigmaErrCollection,
    phiMuCollection, phiSigmaCollection, phiMuErrCollection, phiSigmaErrCollection,
    thetaMuCollection, thetaSigmaCollection, thetaMuErrCollection, thetaSigmaErrCollection, 
    pRecoArr;
    
    for (int bin = 1; bin <= nBins; ++bin)
    {
        TH1D *hProjMom = hMomHist->ProjectionY("MomProjY",bin,bin);
        if (hProjMom->GetEntries() > entriesCutoff)
        {
            pRecoArr.push_back(hMomHist->GetXaxis()->GetBinCenter(bin));
            hProjMom->Fit(fGauss);
            momMuCollection.push_back(fGauss->GetParameter(1));
            momMuErrCollection.push_back(fGauss->GetParError(1));
            momSigmaCollection.push_back(fGauss->GetParameter(2));
            momSigmaErrCollection.push_back(fGauss->GetParError(2));
        }
    }

    TFile *outFile = TFile::Open(outputFile,"recreate");
    TGraphErrors *geMom = new TGraphErrors(momMuCollection.size(),pRecoArr.data(),momMuCollection.data(),nullptr,momMuErrCollection.data());
    geMom->Write();
}