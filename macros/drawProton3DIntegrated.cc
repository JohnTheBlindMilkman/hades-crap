#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph2DErrors.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"

#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"

#include <iostream>

enum class Projection{Out,Side,Long,OutSide,SideLong,OutLong};

struct Proj1D
{
    Projection proj;
    std::string arg;
    std::string name;
    std::size_t num;
};

std::string Get1DProjName(Projection proj)
{
    switch (proj)
    {
        case Projection::Out:
            return "out";
        case Projection::Side:
            return "side";
        case Projection::Long:
            return "long";
    
        default:
            return "";
    }
}

std::pair<std::string,std::string> Get2DProjName(Projection proj)
{
    switch (proj)
    {
        case Projection::OutSide:
            return {"out","side"};
        case Projection::SideLong:
            return {"side","long"};
        case Projection::OutLong:
            return {"out","long"};
    
        default:
            return {};
    }
}

void SetAxesRange(TH3D *data, Projection proj, int nBins)
{
    switch (proj)
    {
        case Projection::Out:
            data->GetXaxis()->SetRange(1,data->GetNbinsX());
            data->GetYaxis()->SetRange(1,nBins);
            data->GetZaxis()->SetRange(1,nBins);
            break;
        case Projection::Side:
            data->GetXaxis()->SetRange(1,nBins);
            data->GetYaxis()->SetRange(1,data->GetNbinsY());
            data->GetZaxis()->SetRange(1,nBins);
            break;
        case Projection::Long:
            data->GetXaxis()->SetRange(1,nBins);
            data->GetYaxis()->SetRange(1,nBins);
            data->GetZaxis()->SetRange(1,data->GetNbinsZ());
            break;
        case Projection::OutSide:
            data->GetXaxis()->SetRange(1,data->GetNbinsX());
            data->GetYaxis()->SetRange(1,data->GetNbinsY());
            data->GetZaxis()->SetRange(1,nBins);
            break;
        case Projection::SideLong:
            data->GetXaxis()->SetRange(1,nBins);
            data->GetYaxis()->SetRange(1,data->GetNbinsY());
            data->GetZaxis()->SetRange(1,data->GetNbinsZ());
            break;
        case Projection::OutLong:
            data->GetXaxis()->SetRange(1,data->GetNbinsX());
            data->GetYaxis()->SetRange(1,nBins);
            data->GetZaxis()->SetRange(1,data->GetNbinsZ());
            break;
        
        default:
            break;
    }
}

TGraph2DErrors* ConvertToGraph(const TH2D *hist)
{
    const int nBinsX = hist->GetNbinsX();
    const int nBinxY = hist->GetNbinsY();
    const std::size_t nBinxTot = nBinsX * nBinxY;
    std::vector<double> xVals(nBinxTot,0.), yVals(nBinxTot,0.), zVals(nBinxTot,0.), xErrs(nBinxTot,0.), yErrs(nBinxTot,0.), zErrs(nBinxTot,0.);

    std::size_t index = 0;
    for (int binX = 1; binX <= nBinsX; ++binX)
        for (int binY = 1; binY <= nBinxY; ++binY)
        {
            xVals[index] = hist->GetXaxis()->GetBinCenter(binX);
            yVals[index] = hist->GetYaxis()->GetBinCenter(binY);
            zVals[index] = hist->GetBinContent(binX,binY);
            zErrs[index] = hist->GetBinError(binX,binY);
            ++index;
        }

    auto graph = new TGraph2DErrors(nBinxTot,xVals.data(),yVals.data(),zVals.data(),zErrs.data(),yErrs.data(),zErrs.data());
    graph->SetName(TString::Format("%s_graph",hist->GetName()));
    graph->SetTitle(hist->GetTitle());

    return graph;
}

void prepareGraph(TH1D *hist, Color_t col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
}

void prepareGraph(TH2D *hist, Color_t col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
    hist->GetZaxis()->SetTitleSize(0.06);
    hist->GetZaxis()->SetLabelSize(0.06);
    hist->GetZaxis()->SetNdivisions(506);
}

void drawProton3DIntegrated()
{
    gStyle->SetOptStat(0);
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12ana_all_25_11_04_processed.root";
    const TString outputFile = "../output/3Dcorr_0_10_cent_Integ.root";
    const std::vector<TString> sProjName = {"out","side","long"};
    const std::array<Proj1D,3> oneDimProjections = {
        Proj1D{Projection::Out,"x",";k^{*}_{out} [MeV/c];C(k^{*}_{out})",0},
        Proj1D{Projection::Side,"y",";k^{*}_{side} [MeV/c];C(k^{*}_{side})",1},
        Proj1D{Projection::Long,"z",";k^{*}_{long} [MeV/c];C(k^{*}_{long})",2}
    };
    const std::array<Proj1D,3> twoDimProjections = {
        Proj1D{Projection::OutSide,"xy",";k^{*}_{side} [MeV/c];k^{*}_{out} [MeV/c];C(k^{*}_{out},k^{*}_{side})",0},
        Proj1D{Projection::SideLong,"yz",";k^{*}_{long} [MeV/c];k^{*}_{side} [MeV/c];C(k^{*}_{side},k^{*}_{long})",1},
        Proj1D{Projection::OutLong,"xz",";k^{*}_{long} [MeV/c];k^{*}_{out} [MeV/c];C(k^{*}_{out},k^{*}_{long})",2}
    };
    const int rebin = 1;
    const int wbin = 1; 

    float norm;
    int binc, binmn, binmx;

    TLine *line = new TLine(-500,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canvInteg = new TCanvas("canvInteg","",1800,600);
    canvInteg->Divide(3,1);
    JJColor::CreateSecondaryWutGradient();
    std::array<TH1D*,oneDimProjections.size()> dataProj1D;

    TPaveText *ptInfo = new TPaveText(.35,.7,.8,.95,"NB NDC");
    ptInfo->SetBorderSize(0);
    ptInfo->SetTextAlign(12);
    ptInfo->SetTextFont(102);
    ptInfo->SetFillStyle(1); // transparemt background
    ptInfo->AddText("proton - proton");
    ptInfo->AddText("0-10%");
    ptInfo->AddText("y_{c.m.} #in (-0.75,0.75)");
    ptInfo->AddText("k_{T} #in (150,1650) [MeV/c]");

    TH3D *hSign = inpFile->Get<TH3D>("hQoslSignInteg");
    TH3D *hBckg = inpFile->Get<TH3D>("hQoslBckgInteg");
    TH3D *hRat3D = new TH3D(*hSign);
    hRat3D->Divide(hBckg);
    hRat3D->SetName("hQoslRatInteg");
    hRat3D->Write();

    if (hRat3D != nullptr)
    {
        for (const auto &[proj,arg,name,num] : oneDimProjections)
        {
            binc = hRat3D->GetXaxis()->FindFixBin(0.0);
            binmx = binc + wbin;
            SetAxesRange(hRat3D,proj,binmx);
            dataProj1D[num] = static_cast<TH1D*>(hRat3D->Project3D(arg.c_str()));
            dataProj1D[num]->SetTitle(name.c_str());
            norm = JJUtils::CF::GetNormByRange(dataProj1D[num],350,500);
            dataProj1D[num]->Rebin(rebin);
            norm *= rebin;
            dataProj1D[num]->Scale(1./norm);
            prepareGraph(dataProj1D[num],JJColor::fWutAllColors[1]); // navy WUT
            dataProj1D[num]->Write();
            TPaveText *pt = new TPaveText(0.3,0.4,0.7,0.6,"NB");
            pt->SetBorderSize(0);
            pt->SetFillStyle(1);
            pt->SetTextAlign(22);
            pt->SetTextFont(102);
            pt->AddText(TString::Format("q_{%s} #in (0,%.2f) [MeV/c]",sProjName[(num + 1) % 3].Data(),dataProj1D[num]->GetXaxis()->GetBinUpEdge(binmx)));
            pt->AddText(TString::Format("q_{%s} #in (0,%.2f) [MeV/c]",sProjName[(num + 2) % 3].Data(),dataProj1D[num]->GetXaxis()->GetBinUpEdge(binmx)));

            canvInteg->cd(num + 1)->SetMargin(0.2,0.02,0.15,0.02);
            dataProj1D[num]->Draw("p e");
            line->Draw("same");
            pt->Draw("same");
            if (num == 0)
                ptInfo->Draw("same");
        }

        canvInteg->Write();

        TCanvas *canv2DInteg = new TCanvas("canv2DInteg","",1800,600);
        canv2DInteg->Divide(3,1);
        JJColor::CreatePrimaryWutGradient();
        std::array<TH2D*,oneDimProjections.size()> dataProj2D;
        for (const auto &[proj,arg,name,num] : twoDimProjections)
        {
            binc = hRat3D->GetXaxis()->FindFixBin(0.0);
            binmx = binc + wbin;

            SetAxesRange(hRat3D,proj,binmx);
            dataProj2D[num] = static_cast<TH2D*>(hRat3D->Project3D(arg.c_str()));
            dataProj2D[num]->SetTitle(name.c_str());
            norm = JJUtils::CF::GetNormByRange(dataProj2D[num],350,500,350,500);
            dataProj2D[num]->Rebin2D(rebin,rebin);
            norm *= rebin*rebin;
            dataProj2D[num]->Scale(1./norm);
            prepareGraph(dataProj2D[num],JJColor::fWutMainColors[1]);
            dataProj2D[num]->Write();

            canv2DInteg->cd(num + 1)->SetMargin(0.2,0.02,0.15,0.02);
            dataProj2D[num]->Draw("colz");
        }

        canv2DInteg->Write();
    }
}
