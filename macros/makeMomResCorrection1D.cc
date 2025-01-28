#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "../Externals/Palettes.hxx"

void SetErrors(TH1D *hout, const TH1D *hNum, TH1D *hDen)
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
}

TH1D* CorrectForMomRes(const TH1D *data, TH1D *momRes)
{
    const int nbins = data->GetNbinsX();
    TH1D *hCorrected = new TH1D(*data);

    for (int i = 1; i <= nbins; ++i)
    {
        hCorrected->SetBinContent(i,data->GetBinContent(i) * momRes->GetBinContent(i));
        SetErrors(hCorrected,data,momRes);
    }

    return hCorrected;
}

void makeMomResCorrection1D()
{
    const TString fileNameExp = "../output/1Dcorr_0_10_cent_Purity.root";
    const TString fileNameMomRes = "../output/1Dmomres_0_10_cent.root";

    const std::vector<std::pair<int,TString> > ktArr{{1,"(200,400)"},{2,"(400,600)"},{3,"(600,800)"},{4,"(800,1000)"},{5,"(1000,1200)"},{6,"(1200,1400)"},{7,"(1400,1600)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.58,-0.35)"},{2,"(-0.35,-0.12)"},{3,"(-0.12,0.12)"},{4,"(0.12,0.35)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 1;

    bool isFirst = true;

    TLine *line = new TLine(0,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFileData = TFile::Open(fileNameExp);
    TFile *inpFileMomRes = TFile::Open(fileNameMomRes);

    TString otpFilePath = fileNameExp;
    otpFilePath.Insert(otpFilePath.Last('.'),"_MomRes");
    TFile *otpFile = TFile::Open(otpFilePath,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.2,0.02,0.15,0.02);
    isFirst = true;
    for (const auto &kt : ktArr)
    {
        TH1D *hMomRes = inpFileMomRes->Get<TH1D>(TString::Format("hQinvMomResKt%d",kt.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatKt%d",kt.first));

        if (hData != nullptr && hMomRes != nullptr)
        {
            TH1D *hCorrect = CorrectForMomRes(hData,hMomRes);
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
            hMomRes->Write();
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
        TH1D *hMomRes = inpFileMomRes->Get<TH1D>(TString::Format("hQinvMomResY%d",y.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatY%d",y.first));

        if (hData != nullptr && hMomRes != nullptr)
        {
            TH1D *hCorrect = CorrectForMomRes(hData,hMomRes);
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
            hMomRes->Write();
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
        TH1D *hMomRes = inpFileMomRes->Get<TH1D>(TString::Format("hQinvMomResPsi%d",psi.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatPsi%d",psi.first));

        if (hData != nullptr && hMomRes != nullptr)
        {
            TH1D *hPur = new TH1D(*hMomRes);
            TH1D *hCorrect = CorrectForMomRes(hData,hMomRes);
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
            hMomRes->Write();
        }
    }
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();

    TH1D *hMomRes = inpFileMomRes->Get<TH1D>("hQinvMomResInteg");
    TH1D *hData = inpFileData->Get<TH1D>("hQinvRatInteg");
    TH1D *hCorrect = CorrectForMomRes(hData,hMomRes);
    hCorrect->Write();
    hMomRes->Write();
}