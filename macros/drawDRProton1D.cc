#include <fstream>

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
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

void AnalyticalDivide(TH1D *hNum, TF1 *fDen)
{
    const int iterMax = hNum->GetNbinsX();
    double min,max;
    fDen->GetRange(min,max);

    for (int i = 1; i <= iterMax; i++)
    {
        if (hNum->GetBinCenter(i) < max)
            hNum->SetBinContent(i,hNum->GetBinContent(i) / fDen->Eval(hNum->GetBinCenter(i)));
        else
            hNum->SetBinContent(i,0.);
    }
}

void drawDRProton1D()
{
    //gStyle->SetPalette(kPastel);
    //JJColor::CreateSecondaryWutGradient();

    const TString inpFile = "../output/1Dcorr_0_10_cent.root";
    const TString simFile = "../output/1Dcorr_0_10_cent_HGeant_fit.root";
    const TString otpFile = "../output/1Dcorr_0_10_cent_DR.root";
    constexpr std::array<int,5> ktArr{1,2,3,4,5};
    constexpr std::array<int,3> yArr{1,2,3};
    constexpr std::array<int,8> psiArr{1,2,3,4,5,6,7,8};
    const int rebin = 1;

    std::vector<TH1D*> hDRkt(ktArr.size(),nullptr), hDRy(yArr.size(),nullptr), hDRpsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hCFkt(ktArr.size(),nullptr), hCFy(yArr.size(),nullptr), hCFpsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hSimkt(ktArr.size(),nullptr), hSimy(yArr.size(),nullptr), hSimpsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hErrkt(ktArr.size(),nullptr), hErry(yArr.size(),nullptr), hErrpsi(psiArr.size(),nullptr);
    std::vector<TF1*> fFitkt(ktArr.size(),nullptr), fFity(yArr.size(),nullptr), fFitpsi(psiArr.size(),nullptr);
    std::vector<TCanvas*> canvkt(ktArr.size(),nullptr), canvy(yArr.size(),nullptr), canvpsi(psiArr.size(),nullptr);
    std::vector<TLegend*> legkt(ktArr.size(),nullptr), legy(yArr.size(),nullptr), legpsi(psiArr.size(),nullptr);
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

        hDRkt[kt-1] = static_cast<TH1D*>(hCFkt[kt-1]->Clone(TString::Format("hQinvDRKt%d",kt)));
        AnalyticalDivide(hDRkt[kt-1],fFitkt[kt-1]);
        SetErrors(hDRkt[kt-1],hCFkt[kt-1],hErrkt[kt-1]);
        hDRkt[kt-1]->GetXaxis()->SetTitleOffset();
        hDRkt[kt-1]->GetXaxis()->SetTitleSize(0.06);
        hDRkt[kt-1]->GetXaxis()->SetLabelSize(0.06);
        hDRkt[kt-1]->GetXaxis()->SetNdivisions(506);
        hDRkt[kt-1]->GetYaxis()->SetTitleSize(0.06);
        hDRkt[kt-1]->GetYaxis()->SetLabelSize(0.06);
        hDRkt[kt-1]->GetYaxis()->SetNdivisions(506);

        canvkt[kt-1] = new TCanvas(TString::Format("cQinvRatKt%d",kt),"",1600,900);
        canvkt[kt-1]->SetMargin(0.15,0.01,0.15,0.01);
        hDRkt[kt-1]->Draw("hist pe");
        hSimkt[kt-1]->Draw("same hist pe");
        hErrkt[kt-1]->Draw("same e4");
        fFitkt[kt-1]->Draw("same c");

        legkt[kt-1] = new TLegend(0.3,0.21,0.3,0.21,"","NB");
        legkt[kt-1]->AddEntry(hDRkt[kt-1],"Corrected result","p");
        legkt[kt-1]->AddEntry(hSimkt[kt-1],"Simulation","p");
        legkt[kt-1]->AddEntry(hErrkt[kt-1],"Fit to simulation","l");
        legkt[kt-1]->Draw("same");

        //canvkt[kt-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
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

        hDRy[y-1] = static_cast<TH1D*>(hCFy[y-1]->Clone(TString::Format("hQinvDRy%d",y)));
        AnalyticalDivide(hDRy[y-1],fFity[y-1]);
        SetErrors(hDRy[y-1],hCFy[y-1],hErry[y-1]);
        hDRy[y-1]->GetXaxis()->SetTitleOffset();
        hDRy[y-1]->GetXaxis()->SetTitleSize(0.06);
        hDRy[y-1]->GetXaxis()->SetLabelSize(0.06);
        hDRy[y-1]->GetXaxis()->SetNdivisions(506);
        hDRy[y-1]->GetYaxis()->SetTitleSize(0.06);
        hDRy[y-1]->GetYaxis()->SetLabelSize(0.06);
        hDRy[y-1]->GetYaxis()->SetNdivisions(506);

        canvy[y-1] = new TCanvas(TString::Format("cQinvRatY%d",y),"",1600,900);
        canvy[y-1]->SetMargin(0.15,0.01,0.15,0.01);
        hDRy[y-1]->Draw("hist pe");
        hSimy[y-1]->Draw("same hist pe");
        hErry[y-1]->Draw("same e3");
        fFity[y-1]->Draw("same c");

        legy[y-1] = new TLegend(0.3,0.21,0.3,0.21,"","NB");
        legy[y-1]->AddEntry(hDRy[y-1],"Corrected result","p");
        legy[y-1]->AddEntry(hSimy[y-1],"Simulation","p");
        legy[y-1]->AddEntry(hErry[y-1],"Fit to simulation","l");
        legy[y-1]->Draw("same");

        //canvy[y-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
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

        hDRpsi[psi-1] = static_cast<TH1D*>(hCFpsi[psi-1]->Clone(TString::Format("hQinvDRpsi%d",psi)));
        AnalyticalDivide(hDRpsi[psi-1],fFitpsi[psi-1]);
        SetErrors(hDRpsi[psi-1],hCFpsi[psi-1],hErrpsi[psi-1]);
        hDRpsi[psi-1]->GetXaxis()->SetTitleOffset();
        hDRpsi[psi-1]->GetXaxis()->SetTitleSize(0.06);
        hDRpsi[psi-1]->GetXaxis()->SetLabelSize(0.06);
        hDRpsi[psi-1]->GetXaxis()->SetNdivisions(506);
        hDRpsi[psi-1]->GetYaxis()->SetTitleSize(0.06);
        hDRpsi[psi-1]->GetYaxis()->SetLabelSize(0.06);
        hDRpsi[psi-1]->GetYaxis()->SetNdivisions(506);

        canvpsi[psi-1] = new TCanvas(TString::Format("cQinvRatPsi%d",psi),"",1600,900);
        canvpsi[psi-1]->SetMargin(0.15,0.01,0.15,0.01);
        hCFpsi[psi-1]->Draw("hist pe");
        hSimpsi[psi-1]->Draw("same hist pe");
        hErrpsi[psi-1]->Draw("same e3");
        fFitpsi[psi-1]->Draw("same c");

        legpsi[psi-1] = new TLegend(0.3,0.21,0.3,0.21,"","NB");
        legpsi[psi-1]->AddEntry(hDRpsi[psi-1],"Corrected result","p");
        legpsi[psi-1]->AddEntry(hSimpsi[psi-1],"Simulation","p");
        legpsi[psi-1]->AddEntry(hErrpsi[psi-1],"Fit to simulation","l");
        legpsi[psi-1]->Draw("same");

        //canvpsi[psi-1]->BuildLegend(0.3,0.21,0.3,0.21,"","p");
        line.Draw("same");
    }

    fOtpFile = TFile::Open(otpFile,"RECREATE");

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    //TLegend *legendKt = new TLegend(.3,.21,.3,.21,"","nb");
    for (const int &kt : ktArr)
    {
        if (hDRkt[kt-1] != nullptr)
        {
            hDRkt[kt-1]->Write();
            hDRkt[kt-1]->SetMarkerColor(JJColor::fWutSecondaryColors[kt-1]);
            hDRkt[kt-1]->SetLineColor(JJColor::fWutSecondaryColors[kt-1]);
            //legendKt->AddEntry(hist,hist->GetTitle(),"p");
            //hist->SetTitle("");

            if (kt == 1)
            {
                hDRkt[kt-1]->GetYaxis()->SetRangeUser(0,1.95);
                hDRkt[kt-1]->GetXaxis()->SetRangeUser(0,495);
                hDRkt[kt-1]->Draw("hist pe");
            }
            else
                hDRkt[kt-1]->Draw("hist pe same");
        }
    }
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canvKt->SetMargin(0.15,0.01,0.15,0.01);
    //legendKt->Draw("same");
    line.Draw("same");
    canvKt->Write();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    canvY->SetMargin(0.15,0.01,0.15,0.01);
    for (const int &y : yArr)
    {
        if (hDRy[y-1] != nullptr)
        {
            hDRy[y-1]->Write();
            hDRy[y-1]->SetMarkerColor(JJColor::fWutSecondaryColors[y-1]);
            hDRy[y-1]->SetLineColor(JJColor::fWutSecondaryColors[y-1]);
            
            if (y == 1)
            {
                hDRy[y-1]->GetYaxis()->SetRangeUser(0,1.95);
                hDRy[y-1]->GetXaxis()->SetRangeUser(0,495);
                hDRy[y-1]->Draw("hist pe");
            }
            else
                hDRy[y-1]->Draw("hist pe same");
        }
    }
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvY->Write();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    canvPsi->SetMargin(0.15,0.01,0.15,0.01);
    for (auto hist : hDRpsi)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFpsi.front())
            {
                hist->GetYaxis()->SetRangeUser(0,1.95);
                hist->GetXaxis()->SetRangeUser(0,495);
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