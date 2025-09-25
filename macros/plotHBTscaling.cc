#include <array>
#include <iostream>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"

#include "../Externals/Palettes.hxx"

double PowerLaw(double *x, double *par)
{
    constexpr double protonMass = 938.27208943;
    // return par[0] * std::pow(std::sqrt(x[0] * x[0] + protonMass * protonMass),-par[1]);
    return par[0] * std::pow(x[0],-par[1]);
}

void MakeNice(TH1D* hist, int col)
{
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(col);
    hist->SetLineColor(col);

    hist->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
}

void MakeNice(TH2D* hist, int col)
{
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(col);
    hist->SetLineColor(col);

    hist->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleOffset();
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
    hist->GetZaxis()->SetTitleSize(0.06);
    hist->GetZaxis()->SetLabelSize(0.06);
    hist->GetZaxis()->SetNdivisions(506);
}

void RemoveUnwantedBins(TH2D *hist, const std::vector<std::pair<int,int> > &unwantedBins)
{
    for (const auto &[xBin,yBin] : unwantedBins)
    {
        hist->SetBinContent(xBin,yBin,0.);
        hist->SetBinError(xBin,yBin,0.);
    }
}

void RemoveUnwantedBins(TH1D *hist, const std::vector<int> &unwantedBins)
{
    for (const auto &bin : unwantedBins)
    {
        hist->SetBinContent(bin,0.);
        hist->SetBinError(bin,0.);
    }
}

TH1D* CloneBoosted(const TH1D *hist, float midRapidity)
{
    const int nBins = hist->GetXaxis()->GetNbins();
    const std::vector<double> data(hist->GetXaxis()->GetXbins()->GetArray(), hist->GetXaxis()->GetXbins()->GetArray() + nBins + 1);
    std::vector<double> dataBoosted;

    std::transform(data.begin(),data.end(),std::back_inserter(dataBoosted),[midRapidity](double x){return x - midRapidity;});

    TH1D *boosted = new TH1D(TString::Format("%s_boosted",hist->GetName()),hist->GetTitle(),nBins,dataBoosted.data());

    for (int i = 1; i <= nBins; ++i)
    {
        boosted->SetBinContent(i,hist->GetBinContent(i));
        boosted->SetBinError(i,hist->GetBinError(i));
    }

    return boosted;
}

TH2D* CloneBoosted(const TH2D *hist, float midRapidity)
{
    const int nBinsX = hist->GetXaxis()->GetNbins();
    const int nBinsY = hist->GetYaxis()->GetNbins();
    const std::vector<double> data(hist->GetYaxis()->GetXbins()->GetArray(),hist->GetYaxis()->GetXbins()->GetArray() + nBinsY + 1);
    std::vector<double> dataBoosted;

    std::transform(data.begin(),data.end(),std::back_inserter(dataBoosted),[midRapidity](double x){return x - midRapidity;});

    TH2D *boosted = new TH2D(TString::Format("%s_boosted",hist->GetName()),hist->GetTitle(),nBinsX,hist->GetXaxis()->GetXbins()->GetArray(),nBinsY,dataBoosted.data());

    for (int i = 1; i <= nBinsX; ++i)
        for (int j = 1; j <= nBinsY; ++j)
        {
            boosted->SetBinContent(i,j,hist->GetBinContent(i,j));
            boosted->SetBinError(i,j,hist->GetBinError(i,j));
        }

    return boosted;
}

TH1D* ExtendXaxis(const TH1D *hist, double factor)
{
    const int nbins = hist->GetNbinsX();
    TH1D *newHist = new TH1D(
        hist->GetName(),
        hist->GetTitle(),
        nbins * factor,
        hist->GetXaxis()->GetXmin(),
        hist->GetXaxis()->GetXmax() * factor);

    for (std::size_t bin = 1; bin <= nbins; ++bin)
    {
        newHist->SetBinContent(bin,hist->GetBinContent(bin));
        newHist->SetBinError(bin,hist->GetBinError(bin));
    }

    return newHist;
}

TH1D* CloneMirrored(const TH1D *hist)
{
    const int nBins = hist->GetXaxis()->GetNbins();
    const std::vector<double> data(hist->GetXaxis()->GetXbins()->GetArray(), hist->GetXaxis()->GetXbins()->GetArray() + nBins + 1);    
    std::vector<double> dataMirrored;

    std::transform(data.begin(),data.end(),std::back_inserter(dataMirrored),[](double x){return -x;});
    std::reverse(dataMirrored.begin(),dataMirrored.end());

    TH1D *mirrored = new TH1D(TString::Format("%s_mirrored",hist->GetName()),hist->GetTitle(),nBins,dataMirrored.data());

    for (int i = 1; i <= nBins; ++i)
    {
        mirrored->SetBinContent(i,hist->GetBinContent((nBins + 1) - i));
        mirrored->SetBinError(i,hist->GetBinError((nBins + 1) - i));
    }

    return mirrored;
}

void plotHBTscaling()
{
    gStyle->SetOptStat(0);

    // setting up constants
    const std::vector<std::vector<int> > rmKts{
        {12},
        {11,12},
        {11,12},
        {10,11,12}
    };
    const std::vector<std::vector<int> > rmYs{
        {},
        {},
        {},
        {}
    };
    const std::vector<std::vector<std::vector<int> > > rmKtsYs{
    {
        {1,11,12},
        {1,11,12},
        {},
        {12},
        {9,12,11,12},
        {9,12,11,12},
        {1,10,11,12},
        {1,2,9,10,11,12},
        {1,2,7,8,9,10,11,12}
    },
    {
        {1,10,11,12},
        {1,10,11,12},
        {11,12},
        {12},
        {8,9,10,11,12},
        {1,10,11,12},
        {1,10},
        {1,10,11,12},
        {1,2,8,9,10,11,12}
    },
    {
        {1,9,10,11,12},
        {1,9,10,11,12},
        {9,10,11,12},
        {9,10,11,12},
        {7,8,9,10,11,12},
        {1,9,10,11,12},
        {1,9,10,11,12},
        {1,2,8,9,10,11,12},
        {1,2,7,8,9,10,11,12}
    },
    {
        {1,9,10,11,12},
        {1,9,10,11,12},
        {9,10,11,12},
        {8,9,10,11,12},
        {7,8,9,10,11,12},
        {9,10,11,12},
        {1,7,8,9,10,11,12},
        {1,7,8,9,10,11,12},
        {1,2,7,8,9,10,11,12}
    }};
    const std::vector<TString> rapIntervals = {"(-0.65,-0.55)","(-0.55,-0.45)","(-0.45,-0.35)","(-0.35,-0.25)","(-0.25,-0.15)","(-0.05,0.05)","(0.05,0.15)","(0.15,0.25)","(0.25,0.35)"};
    const std::vector<int> centrlaities{1,2,3,4};
    const std::vector<TString> centTitle = {"0% - 10%","10% - 20%","20% - 30%","30% - 40%"};

    // TF1 *fit = new TF1("fit",PowerLaw,0,1400,2);
    // fit->SetParameters(3,0.5);
    // fit->SetParLimits(0,1,100);
    // fit->SetParLimits(1,0,1);

    TFile *inpFile;
    TFile *outputFile = TFile::Open("../output/1Dfit.root","RECREATE");

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvKtRinv = new TCanvas("canvKtRinv","",800,450);
    canvKtRinv->SetMargin(0.15,0.05,0.15,0.05);
    for (const auto &cent : centrlaities)
    {
        inpFile = TFile::Open(TString::Format("../output/fitAll1DCF_%d.root",cent));
        TH1D *hRinvTmp = inpFile->Get<TH1D>("radiiKt");
        inpFile->Close();
        if (hRinvTmp != nullptr)
        {
            MakeNice(hRinvTmp, JJColor::fWutSecondaryColors[cent]);
            RemoveUnwantedBins(hRinvTmp,rmKts.at(cent - 1));
            hRinvTmp->SetMarkerSize(1.5);
            hRinvTmp->SetLineWidth(2);

            hRinvTmp->SetTitle(TString::Format("%s;k_{T} [MeV/c];R_{inv}",centTitle[cent-1].Data()));
            if (cent == centrlaities.front())
            {
                hRinvTmp->GetYaxis()->SetRangeUser(2,4);
                hRinvTmp->Draw("p e1x0");
            }
            else
            {
                hRinvTmp->Draw("p e1x0 same");
            }
            //hRinvTmp->Fit(fit,"R");
            outputFile->cd();
            hRinvTmp->Write(TString::Format("radiiKt_%d",cent));
        }
    }
    canvKtRinv->BuildLegend(0.3,0.21,0.3,0.21,"","p");
    outputFile->cd();
    canvKtRinv->Write();

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvYRinv = new TCanvas("canvYRinv","",800,450);
    canvYRinv->SetMargin(0.15,0.05,0.15,0.05);
    TLegend *leg = new TLegend(0.3,0.3,0.5,0.5,"","nbNDC");
    leg->SetNColumns(2);
    // leg->AddEntry(static_cast<TObject*>(nullptr),"Data","");
    // leg->AddEntry(static_cast<TObject*>(nullptr),"Mirrored","");
    // leg->AddEntry(static_cast<TObject*>(nullptr),"","");
    for (const auto &cent : centrlaities)
    {
        inpFile = TFile::Open(TString::Format("../output/fitAll1DCF_%d.root",cent));
        TH1D *hRinvTmp = inpFile->Get<TH1D>("radiiY");
        inpFile->Close();
        if (hRinvTmp != nullptr)
        {
            TH1D *hRinvTmpBoost = CloneBoosted(hRinvTmp,0.74);
            MakeNice(hRinvTmpBoost, JJColor::fWutSecondaryColors[cent]);
            hRinvTmpBoost->SetFillColor(JJColor::fWutSecondaryColors[cent]);
            RemoveUnwantedBins(hRinvTmpBoost,rmYs.at(cent - 1));
            if(cent == centrlaities.front())
            {
                TH1D *hXaxis = new TH1D("hXaxis",";y^{1,2}_{c.m.};R_{inv}",1,-0.7,0.7);
                MakeNice(hXaxis, JJColor::fWutSecondaryColors[cent]);
                hXaxis->GetYaxis()->SetRangeUser(2.2,3.8);
                hXaxis->Draw();
            }
            hRinvTmpBoost->SetMarkerSize(1.5);
            hRinvTmpBoost->SetLineWidth(2);
            hRinvTmpBoost->Draw("p e1x0 same");
            leg->AddEntry(hRinvTmpBoost,"","p");

            TH1D *mirrored = CloneMirrored(hRinvTmpBoost);
            mirrored->SetMarkerStyle(71);
            mirrored->SetMarkerSize(1.5);
            mirrored->SetLineWidth(2);
            mirrored->SetMarkerColor(JJColor::fWutSecondaryColors[cent]);
            mirrored->SetLineColor(JJColor::fWutSecondaryColors[cent]);
            mirrored->Draw("p e1x0 same");
            leg->AddEntry(mirrored,centTitle[cent-1],"p");
            // leg->AddEntry(static_cast<TObject*>(nullptr),centTitle[cent-1],"");

            outputFile->cd();
            hRinvTmpBoost->Write(TString::Format("radiiY_%d",cent));
            mirrored->Write(TString::Format("radiiY_mirrored_%d",cent));
        }    
    }
    // canvYRinv->BuildLegend();
    leg->Draw("same");
    outputFile->cd();
    canvYRinv->Write(); 

    for (const auto &cent : centrlaities)
    {
        JJColor::CreateSecondaryWutGradient();
        TCanvas *canvKtYRinv = new TCanvas(TString::Format("canvKtYRinv_%d",cent),"",800,450);
        canvKtYRinv->SetMargin(0.15,0.05,0.15,0.05);

        inpFile = TFile::Open(TString::Format("../output/fitAll1DCF_%d.root",cent));
        for (const auto &y : {1,2,3,4,5,6,7,8,9})
        {
            TH1D *hRinvTmp = inpFile->Get<TH1D>(TString::Format("radiiKtY_%d",y));
            if (hRinvTmp != nullptr)
            {
                MakeNice(hRinvTmp, JJColor::fWutSecondaryColors11[y]);
                RemoveUnwantedBins(hRinvTmp,rmKtsYs.at(cent - 1).at(y-1));

                hRinvTmp->SetTitle(TString::Format("y_{c.m.}^{1,2} #in %s;k_{T} [MeV/c];R_{inv}",rapIntervals.at(y-1).Data()));
                hRinvTmp->GetYaxis()->SetRangeUser(2.,4.2);
                (y == 1) ? hRinvTmp->Draw("p e1x0") : hRinvTmp->Draw("p e1x0 same");

                outputFile->cd();
                hRinvTmp->Write(TString::Format("radiiKtY_%d_%d",y,cent));
            }  
        }
        inpFile->Close();
        outputFile->cd();  
        canvKtYRinv->BuildLegend(0.3,0.21,0.3,0.21,"","p");
        canvKtYRinv->Write();
    }
}