#include <iostream>
#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"

/* void SetErrors(TH1D *hout, const TH1D *hNum, TH1D *hDen)
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
            vErr = std::sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
} */

TH1D* CorrectForMomRes(const TH1D *data, TH1D *momRes)
{
    const int nbins = data->GetNbinsX();
    TH1D *hCorrected = new TH1D(*data);

    for (int i = 1; i <= nbins; ++i)
    {
        hCorrected->SetBinContent(i,data->GetBinContent(i) * momRes->GetBinContent(i));
        //std::cout << "PCC: " << JJUtils::Generic::Detail::GetPCC(data,momRes) << "\n";
    }

    JJUtils::Generic::SetErrorsMultiply(hCorrected,data,momRes);

    return hCorrected;
}

void makeMomResCorrection1D()
{
    const TString fileNameExp = "../output/1Dcorr_0_10_cent_Purity.root";
    const TString fileNameMomRes = "../output/1Dmomres_0_10_cent_SMASH.root";

    const std::vector<std::pair<int,TString> > ktArr{{1,"(0,200)"},{2,"(200,400)"},{3,"(400,600)"},{4,"(600,800)"},{5,"(800,1000)"},{6,"(1000,1200)"},{7,"(1200,1400)"},{8,"(1400,1600)"},{9,"(1600,1800)"},{10,"(1800,2000)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.65,-0.55)"},{2,"(-0.55,-0.45)"},{3,"(-0.45,-0.35)"},{4,"(-0.35,-0.25)"},{5,"(-0.25,-0.15)"},{6,"(-0.15,-0.05)"},{7,"(-0.05,0.05)"},{8,"(0.05,0.15)"},{9,"(0.15,0.25)"},{10,"(0.25,0.35)"},{11,"(0.35,0.45)"},{12,"(0.45,0.55)"},{13,"(0.55,0.65)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 1;

    bool isFirst = true;

    TLine *line = new TLine(0,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFileData = TFile::Open(fileNameExp);
    TFile *inpFileMomRes = TFile::Open(fileNameMomRes);
    TH1D *hMomResFactor = inpFileMomRes->Get<TH1D>("hRatio");

    TString otpFilePath = fileNameExp;
    otpFilePath.Insert(otpFilePath.Last('.'),"_MomRes");
    TFile *otpFile = TFile::Open(otpFilePath,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.2,0.02,0.15,0.02);
    isFirst = true;
    for (const auto &kt : ktArr)
    {
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatKt%d",kt.first));

        if (hData != nullptr)
        {
            TH1D *hCorrect = CorrectForMomRes(hData,hMomResFactor);
            if (isFirst)
            {
                hCorrect->Draw("hist pe");
                isFirst = false;
            }
            else
            {
                hCorrect->Draw("same hist pe");
            }

            hCorrect->Write();
        }
    }
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvKt->Write();

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    canvY->SetMargin(0.2,0.02,0.15,0.02);
    isFirst = true;
    for (const auto &y : yArr)
    {
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatY%d",y.first));

        if (hData != nullptr)
        {
            TH1D *hCorrect = CorrectForMomRes(hData,hMomResFactor);
            if (isFirst)
            {
                hCorrect->Draw("hist pe");
                isFirst = false;
            }
            else
            {
                hCorrect->Draw("same hist pe");
            }

            hCorrect->Write();
        }
    }
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvY->Write();

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    canvPsi->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &psi : psiArr)
    {
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatPsi%d",psi.first));

        if (hData != nullptr)
        {
            TH1D *hCorrect = CorrectForMomRes(hData,hMomResFactor);
            if (isFirst)
            {
                hCorrect->Draw("hist pe");
                isFirst = false;
            }
            else
            {
                hCorrect->Draw("same hist pe");
            }

            hCorrect->Write();
        }
    }
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();

    TH1D *hData = inpFileData->Get<TH1D>("hQinvRatInteg");
    TH1D *hCorrect = CorrectForMomRes(hData,hMomResFactor);
    hCorrect->Write();
}