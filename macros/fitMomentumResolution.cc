#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

void fitMomentumResolution()
{
    const TString fileName = "../slurmOutput/apr12sim_all_25_06_03.root";
    // const TString fileName = "../slurmOutput/apr12sim_all_25_09_11.root";
    const TString outputFile = "../output/momentum_resolution_old.root";
    const TString inpMomHist = "hMomResolution";
    const TString inpPhiHist = "hPhiResolution";
    const TString inpThetaHist = "hThetaResolution";
    constexpr int entriesCutoff{10000};
    constexpr double nSig = 1.5;

    TFile *inpFile = TFile::Open(fileName);

    TH2D *hMomHist = inpFile->Get<TH2D>(inpMomHist);
    TH2D *hPhiHist = inpFile->Get<TH2D>(inpPhiHist);
    TH2D *hThetaHist = inpFile->Get<TH2D>(inpThetaHist);

    TF1 *fGauss = new TF1("fGauss","gaus",-5,5);

    const int nBins = hMomHist->GetNbinsX();
    std::vector<double> momMuCollection, momSigmaCollection, momMuErrCollection, momSigmaErrCollection,
    phiMuCollection, phiSigmaCollection, phiMuErrCollection, phiSigmaErrCollection,
    thetaMuCollection, thetaSigmaCollection, thetaMuErrCollection, thetaSigmaErrCollection, pRecoArr;
    std::vector<TH1D*> momSlice, phiSlice, thetaSlice;
    double mean = 0, stdev = 0, mode = 0;

    TH1D *hChi2Mom = new TH1D("hChi2Mom",";p_{kine};#chi^{2}/NDF",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hChi2Phi = new TH1D("hChi2Phi",";#phi_{kine};#chi^{2}/NDF",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hChi2Theta = new TH1D("hChi2Theta",";#theta_{kine};#chi^{2}/NDF",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    
    TH1D *hSigDiffMom = new TH1D("hSigDiffMom",";p_{kine};#frac{|#sigma_{fit} - #sigma_{hist}|}{#sigma_{hist}} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hSigDiffPhi = new TH1D("hSigDiffPhi",";#phi_{kine};#frac{|#sigma_{fit} - #sigma_{hist}|}{#sigma_{hist}} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hSigDiffTheta = new TH1D("hSigDiffTheta",";#theta_{kine};#frac{|#sigma_{fit} - #sigma_{hist}|}{#sigma_{hist}} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hMeanDiffMom = new TH1D("hMeanDiffMom",";p_{kine};#frac{|#mu_{fit} - #mu_{hist}|}{|#mu_{hist}|} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hMeanDiffPhi = new TH1D("hMeanDiffPhi",";#phi_{kine};#frac{|#mu_{fit} - #mu_{hist}|}{|#mu_{hist}|} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    TH1D *hMeanDiffTheta = new TH1D("hMeanDiffTheta",";#theta_{kine};#frac{|#mu_{fit} - #mu_{hist}|}{|#mu_{hist}|} [%]",nBins,hMomHist->GetXaxis()->GetXmin(),hMomHist->GetXaxis()->GetXmax());
    
    for (int bin = 1; bin <= nBins; ++bin)
    {
        TH1D *hProjMom = hMomHist->ProjectionY(TString::Format("MomProjY_%d",bin),bin,bin);
        if (hProjMom->GetEntries() > entriesCutoff)
        {
            pRecoArr.push_back(hMomHist->GetXaxis()->GetBinCenter(bin));
            mean = hProjMom->GetMean();
            stdev = hProjMom->GetStdDev();
            mode = hProjMom->GetXaxis()->GetBinCenter(hProjMom->GetMaximumBin());
            fGauss->SetRange(mode - (nSig * stdev),mode + (nSig * stdev));
            TFitResultPtr fitPtr = hProjMom->Fit(fGauss,"SME");
            momMuCollection.push_back(fitPtr->Parameter(1));
            momMuErrCollection.push_back(fitPtr->ParError(1));
            momSigmaCollection.push_back(fitPtr->Parameter(2));
            momSigmaErrCollection.push_back(fitPtr->ParError(2));
            momSlice.push_back(hProjMom);

            hChi2Mom->SetBinContent(bin,fitPtr->Chi2() / fitPtr->Ndf());
            hSigDiffMom->SetBinContent(bin,std::abs(fitPtr->Parameter(2) - stdev) / stdev * 100);
            hMeanDiffMom->SetBinContent(bin,std::abs(fitPtr->Parameter(1) - mean) / std::abs(mean) * 100);
        }

        TH1D *hProjPhi = hPhiHist->ProjectionY(TString::Format("PhiProjY_%d",bin),bin,bin);
        if (hProjPhi->GetEntries() > entriesCutoff)
        {
            mean = hProjPhi->GetMean();
            stdev = hProjPhi->GetStdDev();
            mode = hProjPhi->GetXaxis()->GetBinCenter(hProjPhi->GetMaximumBin());
            fGauss->SetRange(mode - (nSig * stdev),mode + (nSig * stdev));
            TFitResultPtr fitPtr = hProjPhi->Fit(fGauss,"SME");
            phiMuCollection.push_back(fitPtr->Parameter(1));
            phiMuErrCollection.push_back(fitPtr->ParError(1));
            phiSigmaCollection.push_back(fitPtr->Parameter(2));
            phiSigmaErrCollection.push_back(fitPtr->ParError(2));
            phiSlice.push_back(hProjPhi);

            hChi2Phi->SetBinContent(bin,fitPtr->Chi2() / fitPtr->Ndf());
            hSigDiffPhi->SetBinContent(bin,std::abs(fitPtr->Parameter(2) - stdev) / stdev * 100);
            hMeanDiffPhi->SetBinContent(bin,std::abs(fitPtr->Parameter(1) - mean) / std::abs(mean) * 100);
        }

        TH1D *hProjTheta = hThetaHist->ProjectionY(TString::Format("ThetaProjY_%d",bin),bin,bin);
        if (hProjTheta->GetEntries() > entriesCutoff)
        {
            mean = hProjTheta->GetMean();
            stdev = hProjTheta->GetStdDev();
            mode = hProjTheta->GetXaxis()->GetBinCenter(hProjTheta->GetMaximumBin());
            fGauss->SetRange(mode - (nSig * stdev),mode + (nSig * stdev));
            TFitResultPtr fitPtr = hProjTheta->Fit(fGauss,"SME");
            thetaMuCollection.push_back(fitPtr->Parameter(1));
            thetaMuErrCollection.push_back(fitPtr->ParError(1));
            thetaSigmaCollection.push_back(fitPtr->Parameter(2));
            thetaSigmaErrCollection.push_back(fitPtr->ParError(2));
            thetaSlice.push_back(hProjTheta);

            hChi2Theta->SetBinContent(bin,fitPtr->Chi2() / fitPtr->Ndf());
            hSigDiffTheta->SetBinContent(bin,std::abs(fitPtr->Parameter(2) - stdev) / stdev * 100);
            hMeanDiffTheta->SetBinContent(bin,std::abs(fitPtr->Parameter(1) - mean) / std::abs(mean) * 100);
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

    hChi2Mom->Write();
    hChi2Phi->Write();
    hChi2Theta->Write();

    hSigDiffMom->Write();
    hSigDiffPhi->Write();
    hSigDiffTheta->Write();
    hMeanDiffMom->Write();
    hMeanDiffPhi->Write();
    hMeanDiffTheta->Write();

    for (const auto &hist : momSlice) hist->Write();
    for (const auto &hist : phiSlice) hist->Write();
    for (const auto &hist : thetaSlice) hist->Write();
}