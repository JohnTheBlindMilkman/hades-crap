#include <fstream>

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "Palettes.hxx"

void SetErrors(TH1D *hout, TH1D *hNum, TH1D *hDen)
{
    const int iterMax = hout->GetNbinsX();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= iterMax; i++)
    {
        vErr = 0;
        vNum = hNum->GetBinContent(i);
        eNum = hNum->GetBinError(i);
        vDen = hDen->GetBinContent(hDen->FindBin(hNum->GetBinCenter(i)));
        eDen = hDen->GetBinError(hDen->FindBin(hNum->GetBinCenter(i)));

        // propagation of uncertainty for a function num/den with onclusion of the correlation between the constituents
        if (fabs(vDen) > 0.)
            vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
}

void drawDRProton1D()
{
    gStyle->SetPalette(kPastel);

    const TString inpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent.root";
    const TString simFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_HGeant_fit.root";
    const TString otpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR.root";
    constexpr std::array<int,5> ktArr{1,2,3,4,5};
    constexpr std::array<int,3> yArr{1,2,3};
    constexpr std::array<int,8> psiArr{1,2,3,4,5,6,7,8};
    const int rebin = 1;

    std::vector<TH1D*> hCFkt(ktArr.size(),nullptr), hCFy(yArr.size(),nullptr), hCFpsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hSimkt(ktArr.size(),nullptr), hSimy(yArr.size(),nullptr), hSimpsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hErrkt(ktArr.size(),nullptr), hErry(yArr.size(),nullptr), hErrpsi(psiArr.size(),nullptr);
    std::vector<TF1*> fFitkt(ktArr.size(),nullptr), fFity(yArr.size(),nullptr), fFitpsi(psiArr.size(),nullptr);
    std::vector<TCanvas*> canvkt(ktArr.size(),nullptr), canvy(yArr.size(),nullptr), canvpsi(psiArr.size(),nullptr);
    TFile *fInpData,*fInpSim,*fOtpFile;
    TLine line(0,1,3000,1);
    line.SetLineColor(kGray);
    line.SetLineStyle(kDashed);

    fInpData = TFile::Open(inpFile);
    fInpSim = TFile::Open(simFile);

    for(const int &kt : ktArr)
    {
        TString histName = TString::Format("hQinvRatKt%d",kt);
        TString errName = TString::Format("hQinvErrKt%d",kt);
        TString fitName = TString::Format("fQinvFitKt%d",kt);

        fInpSim->cd();
        hSimkt[kt-1] = fInpSim->Get<TH1D>(histName);
        hSimkt[kt-1]->SetMarkerColor(TColor::GetColor(JJColor::greenWut.AsHexString()));
        hSimkt[kt-1]->SetLineColor(TColor::GetColor(JJColor::greenWut.AsHexString()));

        hErrkt[kt-1] = fInpSim->Get<TH1D>(errName);
        hErrkt[kt-1]->SetFillColorAlpha(TColor::GetColor(JJColor::fireRedWut.AsHexString()),.75);
        //hErrkt[kt-1]->SetFillStyle(3001);
        hErrkt[kt-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fFitkt[kt-1] = fInpSim->Get<TF1>(fitName);
        fFitkt[kt-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fInpData->cd();
        hCFkt[kt-1] = fInpData->Get<TH1D>(histName);
        hCFkt[kt-1]->SetMarkerColor(TColor::GetColor(JJColor::navyWut.AsHexString()));
        hCFkt[kt-1]->SetLineColor(TColor::GetColor(JJColor::navyWut.AsHexString()));

        hCFkt[kt-1]->Divide(hSimkt[kt-1]);
        SetErrors(hCFkt[kt-1],hCFkt[kt-1],hErrkt[kt-1]);

        canvkt[kt-1] = new TCanvas(TString::Format("cQinvRatKt%d",kt),"",1600,900);
        hCFkt[kt-1]->Draw("hist pe");
        hSimkt[kt-1]->Draw("same hist pe");
        hErrkt[kt-1]->Draw("same e3");
        fFitkt[kt-1]->Draw("same c");

        canvkt[kt-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
        line.Draw("same");
    }
    for(const int &y : yArr)
    {
        TString histName = TString::Format("hQinvRatY%d",y);
        TString errName = TString::Format("hQinvErrY%d",y);
        TString fitName = TString::Format("fQinvFitY%d",y);

        fInpSim->cd();
        hSimy[y-1] = fInpSim->Get<TH1D>(histName);
        hSimy[y-1]->SetMarkerColor(TColor::GetColor(JJColor::greenWut.AsHexString()));
        hSimy[y-1]->SetLineColor(TColor::GetColor(JJColor::greenWut.AsHexString()));

        hErry[y-1] = fInpSim->Get<TH1D>(errName);
        hErry[y-1]->SetFillColorAlpha(TColor::GetColor(JJColor::fireRedWut.AsHexString()),.75);
        //hErry[y-1]->SetFillStyle(3013);
        hErry[y-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fFity[y-1] = fInpSim->Get<TF1>(fitName);
        fFity[y-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fInpData->cd();
        hCFy[y-1] = fInpData->Get<TH1D>(histName);
        hCFy[y-1]->SetMarkerColor(TColor::GetColor(JJColor::navyWut.AsHexString()));
        hCFy[y-1]->SetLineColor(TColor::GetColor(JJColor::navyWut.AsHexString()));

        hCFy[y-1]->Divide(hSimy[y-1]);
        SetErrors(hCFy[y-1],hCFy[y-1],hErry[y-1]);

        canvy[y-1] = new TCanvas(TString::Format("cQinvRatY%d",y),"",1600,900);
        hCFy[y-1]->Draw("hist pe");
        hSimy[y-1]->Draw("same hist pe");
        hErry[y-1]->Draw("same e3");
        fFity[y-1]->Draw("same c");

        canvy[y-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
        line.Draw("same");
    }
    for(const int &psi : psiArr)
    {
        TString histName = TString::Format("hQinvRatPsi%d",psi);
        TString errName = TString::Format("hQinvErrPsi%d",psi);
        TString fitName = TString::Format("fQinvFitPsi%d",psi);

        fInpSim->cd();
        hSimpsi[psi-1] = fInpSim->Get<TH1D>(histName);
        hSimpsi[psi-1]->SetMarkerColor(TColor::GetColor(JJColor::greenWut.AsHexString()));
        hSimpsi[psi-1]->SetLineColor(TColor::GetColor(JJColor::greenWut.AsHexString()));

        hErrpsi[psi-1] = fInpSim->Get<TH1D>(errName);
        hErrpsi[psi-1]->SetFillColorAlpha(TColor::GetColor(JJColor::fireRedWut.AsHexString()),0.75);
        //hErrpsi[psi-1]->SetFillStyle(3013);
        hErrpsi[psi-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fFitpsi[psi-1] = fInpSim->Get<TF1>(fitName);
        fFitpsi[psi-1]->SetLineColor(TColor::GetColor(JJColor::fireRedWut.AsHexString()));

        fInpData->cd();
        hCFpsi[psi-1] = fInpData->Get<TH1D>(histName);
        hCFpsi[psi-1]->SetMarkerColor(TColor::GetColor(JJColor::navyWut.AsHexString()));
        hCFpsi[psi-1]->SetLineColor(TColor::GetColor(JJColor::navyWut.AsHexString()));

        hCFpsi[psi-1]->Divide(hSimpsi[psi-1]);
        SetErrors(hCFpsi[psi-1],hCFpsi[psi-1],hErrpsi[psi-1]);

        canvpsi[psi-1] = new TCanvas(TString::Format("cQinvRatPsi%d",psi),"",1600,900);
        hCFpsi[psi-1]->Draw("hist pe");
        hSimpsi[psi-1]->Draw("same hist pe");
        hErrpsi[psi-1]->Draw("same e3");
        fFitpsi[psi-1]->Draw("same c");

        canvpsi[psi-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
        line.Draw("same");
    }

    fOtpFile = TFile::Open(otpFile,"RECREATE");

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    for (auto hist : hCFkt)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFkt.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvKt->Write();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    for (auto hist : hCFy)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFy.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvY->Write();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    for (auto hist : hCFpsi)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFpsi.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvPsi->Write();

    for (const auto &kt : ktArr)
        canvkt[kt-1]->Write();
    for (const auto &y : yArr)
        canvy[y-1]->Write();
    for (const auto &psi : psiArr)
        canvpsi[psi-1]->Write();
}