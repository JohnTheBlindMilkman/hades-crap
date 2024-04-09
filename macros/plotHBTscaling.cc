#include <array>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"

#include "Palettes.hxx"

void MakeNice(TH1D* hist)
{
    hist->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(JJColor::fWutAllColors[1]); // navy WUT
    hist->SetFillColor(JJColor::fWutAllColors[1]); // navy WUT
}

void plotHBTscaling()
{
    gStyle->SetOptStat(0);

    // setting up constants
    // data from fit_proton.C, using D. Wielanek's HAL for fitting
    constexpr unsigned long ktArrSize{5}; 
    float ktVal[] = {150,450,750,1050,1350,1650};
    constexpr std::array<float,ktArrSize> ktRinvVal{3.093,2.456,2.342,1.88,1.695};
    constexpr std::array<float,ktArrSize> ktRinvErr{0.005,0.024,0.019,0.024,0.035};
    constexpr std::array<float,ktArrSize> ktLambdaVal{0.241,0.187,0.254,0.19,0.188};
    constexpr std::array<float,ktArrSize> ktLambdaErr{0.001,0.006,0.006,0.006,0.009};

    constexpr unsigned long yArrSize{3}; 
    float yVal[] = {-0.75,-0.25,0.25,0.75};
    constexpr std::array<float,yArrSize> yRinvVal{1.733,2.926,2.721};
    constexpr std::array<float,yArrSize> yRinvErr{0.015,0.014,0.014};
    constexpr std::array<float,yArrSize> yLambdaVal{0.089,0.522,0.496};
    constexpr std::array<float,yArrSize> yLambdaErr{0.002,0.008,0.008};

    TH1D *hKtRinv = new TH1D("hKtRinv",";k_{T} [MeV/c];R_{inv} [fm]",ktArrSize,ktVal);
    TH1D *hKtLambda = new TH1D("hKtLambda",";k_{T} [MeV/c];#lambda [a.u.]",ktArrSize,ktVal);

    for (std::size_t kt = 1; kt <= ktArrSize; ++kt)
    {
        hKtRinv->SetBinContent(kt,ktRinvVal[kt-1]);
        hKtRinv->SetBinError(kt,ktRinvErr[kt-1]);

        hKtLambda->SetBinContent(kt,ktLambdaVal[kt-1]);
        hKtLambda->SetBinError(kt,ktLambdaErr[kt-1]);
    }

    MakeNice(hKtRinv);
    MakeNice(hKtLambda);

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvKtRinv = new TCanvas("canvKtRinv","",800,450);
    canvKtRinv->SetMargin(0.2,0.02,0.15,0.02);
    hKtRinv->Draw("p e x0");

    TCanvas *canvKtLambda = new TCanvas("canvKtLambda","",800,450);
    canvKtLambda->SetMargin(0.2,0.02,0.15,0.02);
    hKtLambda->Draw("p e x0");

    TH1D *hYRinv = new TH1D("hYRinv",";y_{c.m.} [a.u.];R_{inv} [fm]",yArrSize,yVal);
    TH1D *hYLambda = new TH1D("hYLambda",";y_{c.m.} [a.u.];#lambda [a.u.]",yArrSize,yVal);

    for (std::size_t y = 1; y <= yArrSize; ++y)
    {
        hYRinv->SetBinContent(y,yRinvVal[y-1]);
        hYRinv->SetBinError(y,yRinvErr[y-1]);

        hYLambda->SetBinContent(y,yLambdaVal[y-1]);
        hYLambda->SetBinError(y,yLambdaErr[y-1]);
    }

    MakeNice(hYRinv);
    MakeNice(hYLambda);

    TCanvas *canvYRinv = new TCanvas("canvYRinv","",800,450);
    canvYRinv->SetMargin(0.2,0.02,0.15,0.02);
    hYRinv->Draw("p e x0");

    TCanvas *canvYLambda = new TCanvas("canvYLambda","",800,450);
    canvYLambda->SetMargin(0.2,0.02,0.15,0.02);
    hYLambda->Draw("p e x0");

    TFile *outputFile = TFile::Open("../output/1Dfit_0_10.root","RECREATE");
    canvKtRinv->Write();
    canvKtLambda->Write();
    canvYRinv->Write();
    canvYLambda->Write();

    hKtRinv->Write();
    hKtLambda->Write();
    hYRinv->Write();
    hYLambda->Write();
}