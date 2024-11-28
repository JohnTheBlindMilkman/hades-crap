#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphErrors.h"

void fitMomentumResolution()
{
    const TString fileName = "../slurmOutput/apr12momres_all_24_10_08.root";
    const TString outputFile = "../output/momentum_resolution.root";
    const TString inpMomHist = "hMomResolution";
    const TString inpPhiHist = "hPhiResolution";
    const TString inpThetaHist = "hThetaResolution";
    constexpr int entriesCutoff{10000};
    constexpr bool fitSlices{true};
    constexpr bool fitGaussianParams{true};

    TFile *inpFile = TFile::Open(fileName);

    TH2D *hMomHist = inpFile->Get<TH2D>(inpMomHist);
    TH2D *hPhiHist = inpFile->Get<TH2D>(inpPhiHist);
    TH2D *hThetaHist = inpFile->Get<TH2D>(inpThetaHist);

    TF1 *fGauss = new TF1("fGauss","gaus",hMomHist->GetYaxis()->GetXmin(),hMomHist->GetYaxis()->GetXmax());

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

        TH1D *hProjPhi = hPhiHist->ProjectionY("PhiProjY",bin,bin);
        if (hProjPhi->GetEntries() > entriesCutoff)
        {
            hProjPhi->Fit(fGauss);
            phiMuCollection.push_back(fGauss->GetParameter(1));
            phiMuErrCollection.push_back(fGauss->GetParError(1));
            phiSigmaCollection.push_back(fGauss->GetParameter(2));
            phiSigmaErrCollection.push_back(fGauss->GetParError(2));
        }

        TH1D *hProjTheta = hThetaHist->ProjectionY("ThetaProjY",bin,bin);
        if (hProjTheta->GetEntries() > entriesCutoff)
        {
            hProjTheta->Fit(fGauss);
            thetaMuCollection.push_back(fGauss->GetParameter(1));
            thetaMuErrCollection.push_back(fGauss->GetParError(1));
            thetaSigmaCollection.push_back(fGauss->GetParameter(2));
            thetaSigmaErrCollection.push_back(fGauss->GetParError(2));
        }
    }

    TFile *outFile = TFile::Open(outputFile,"recreate");
    TGraphErrors *geMomMu = new TGraphErrors(momMuCollection.size(),pRecoArr.data(),momMuCollection.data(),nullptr,momMuErrCollection.data());
    geMomMu->SetName("geMomMu");
    TGraphErrors *geMomSig = new TGraphErrors(momSigmaCollection.size(),pRecoArr.data(),momSigmaCollection.data(),nullptr,momSigmaErrCollection.data());
    geMomSig->SetName("geMomSig");
    TGraphErrors *gePhiMu = new TGraphErrors(phiMuCollection.size(),pRecoArr.data(),phiMuCollection.data(),nullptr,phiMuErrCollection.data());
    gePhiMu->SetName("gePhiMu");
    TGraphErrors *gePhiSig = new TGraphErrors(phiSigmaCollection.size(),pRecoArr.data(),phiSigmaCollection.data(),nullptr,phiSigmaErrCollection.data());
    gePhiSig->SetName("gePhiSig");
    TGraphErrors *geThetaMu = new TGraphErrors(thetaMuCollection.size(),pRecoArr.data(),thetaMuCollection.data(),nullptr,thetaMuErrCollection.data());
    geThetaMu->SetName("geThetaMu");
    TGraphErrors *geThetaSig = new TGraphErrors(thetaSigmaCollection.size(),pRecoArr.data(),thetaSigmaCollection.data(),nullptr,thetaSigmaErrCollection.data());
    geThetaSig->SetName("geThetaSig");

    TF1 *fMomMuFit = new TF1("fMomMuFit","[0] + [1]/x + [2]/(x**2) + pol3(3)",0,3500);
    fMomMuFit->SetParameters(1.,1.,1.,1.,1.,1.,1.);
    TF1 *fMomSigFit = new TF1("fMomSigFit","pol7",0,3500);
    TF1 *fPhiMuFit = new TF1("fPhiMuFit","[0] + [1]/x",0,3500);
    fPhiMuFit->SetParameters(1.,1.);
    TF1 *fPhiSigFit = new TF1("fPhiSigFit","[0] + [1]/x + [2]/(x**2)",0,3500);
    fPhiMuFit->SetParameters(1.,1.,1.);
    TF1 *fThetaMuFit = new TF1("fThetaMuFit","[0] + [1]/x + [2]/(x**2) + pol3(3)",0,3500);
    fThetaMuFit->SetParameters(1.,1.,1.,1.,1.,1.,1.);
    TF1 *fThetaSigFit = new TF1("fThetaSigFit","[0] + [1]/x + [2]/(x**2)",0,3500);
    fThetaSigFit->SetParameters(1.,1.,1.);

    geMomMu->Fit(fMomMuFit,"ME");
    geMomSig->Fit(fMomSigFit,"ME");
    gePhiMu->Fit(fPhiMuFit,"ME");
    gePhiSig->Fit(fPhiSigFit,"ME");
    geThetaMu->Fit(fThetaMuFit,"ME");
    geThetaSig->Fit(fThetaSigFit,"ME");

    geMomMu->Write();
    geMomSig->Write();
    gePhiMu->Write();
    gePhiSig->Write();
    geThetaMu->Write();
    geThetaSig->Write();

    fMomMuFit->Write();
    fMomSigFit->Write();
    fPhiMuFit->Write();
    fPhiSigFit->Write();
    fThetaMuFit->Write();
    fThetaSigFit->Write();
}