#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"

TH1D* CorrectForPurity(const TH1D *hDataRaw, TH1D *hPurity)
{
    const int nbins = hDataRaw->GetNbinsX();
    TH1D *hCorrected = new TH1D(*hDataRaw);

    for (int i = 1; i <= nbins; ++i)
    {
        //hCorrected->SetBinContent(i,hDataRaw->GetBinContent(i) * hPurity->GetBinContent(i));
        // here we assume that the purity plot has 1% uncertainty
        //hCorrected->SetBinError(i,hCorrected->GetBinContent(i) * std::sqrt(std::pow(hDataRaw->GetBinError(i)/hDataRaw->GetBinContent(i),2) + std::pow(0.01/hPurity->GetBinContent(i),2)));
    }

    //JJUtils::Generic::SetErrorsMultiply(hCorrected,hDataRaw,hPurity);

    return hCorrected;
}

void makePurityCorrection1D()
{
    const TString fileNameExp = "../output/1Dcorr_0_10_cent.root";
    const TString fileNameExpInteg = "../output/1Dcorr_0_10_cent_Integ.root";
    const TString fileNamePur = "../slurmOutput/apr12pur_all_24_11_20_1_processed.root";

    const std::vector<std::pair<int,TString> > ktArr{{1,"(0,200)"},{2,"(200,400)"},{3,"(400,600)"},{4,"(600,800)"},{5,"(800,1000)"},{6,"(1000,1200)"},{7,"(1200,1400)"},{8,"(1400,1600)"},{9,"(1600,1800)"},{10,"(1800,2000)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.65,-0.55)"},{2,"(-0.55,-0.45)"},{3,"(-0.45,-0.35)"},{4,"(-0.35,-0.25)"},{5,"(-0.25,-0.15)"},{6,"(-0.15,-0.05)"},{7,"(-0.05,0.05)"},{8,"(0.05,0.15)"},{9,"(0.15,0.25)"},{10,"(0.25,0.35)"},{11,"(0.35,0.45)"},{12,"(0.45,0.55)"},{13,"(0.55,0.65)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 1;

    bool isFirst = true;

    TLine *line = new TLine(0,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFileData = TFile::Open(fileNameExp);
    TFile *inpFileDataInteg = TFile::Open(fileNameExpInteg);
    TFile *inpFilePur = TFile::Open(fileNamePur);

    TString otpFilePath = fileNameExp;
    otpFilePath.Insert(otpFilePath.Last('.'),"_Purity");
    TFile *otpFile = TFile::Open(otpFilePath,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.2,0.02,0.15,0.02);
    isFirst = true;
    for (const auto &kt : ktArr)
    {
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignKt%d",kt.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgKt%d",kt.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatKt%d",kt.first));

        if (hData != nullptr && hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            JJUtils::Generic::SetErrorsDivide(hPur,hPurNum,hPurDen);
            hPur->SetName(TString::Format("hQinvPurKt%d",kt.first));
            hPur->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF(q_{inv})",kt.second.Data()));
            TH1D *hCorrect = CorrectForPurity(hData,hPur);
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
            hPur->Write();
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
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignY%d",y.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgY%d",y.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatY%d",y.first));

        if (hData != nullptr && hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            JJUtils::Generic::SetErrorsDivide(hPur,hPurNum,hPurDen);
            hPur->SetName(TString::Format("hQinvPurY%d",y.first));
            hPur->SetTitle(TString::Format("y #in %s",y.second.Data()));
            TH1D *hCorrect = CorrectForPurity(hData,hPur);
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
            hPur->Write();
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
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignPsi%d",psi.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgPsi%d",psi.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatPsi%d",psi.first));

        if (hData != nullptr && hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            JJUtils::Generic::SetErrorsDivide(hPur,hPurNum,hPurDen);
            hPur->SetName(TString::Format("hQinvPurPsi%d",psi.first));
            hPur->SetTitle(TString::Format("#phi - #Psi_{E.P.} #in %s",psi.second.Data()));
            TH1D *hCorrect = CorrectForPurity(hData,hPur);
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
            hPur->Write();
        }
    }
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvPsi->Write();

    TH1D *hPurNum = inpFilePur->Get<TH1D>("hQinvSignInteg");
    TH1D *hPurDen = inpFilePur->Get<TH1D>("hQinvBckgInteg");
    TH1D *hData = inpFileDataInteg->Get<TH1D>("hQinvRatInteg");
    TH1D *hPur = new TH1D(*hPurNum);
    hPur->Divide(hPurDen);
    JJUtils::Generic::SetErrorsDivide(hPur,hPurNum,hPurDen);
    hPur->SetName("hQinvPurInterg");
    TH1D *hCorrect = CorrectForPurity(hData,hPur);
    hCorrect->Write();
}