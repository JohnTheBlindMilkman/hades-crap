#include "TString.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "../Externals/Palettes.hxx"

void makePurityCorrection1D()
{
    const TString fileNameExp = "../output/1Dcorr_0_10_cent.root";
    const TString fileNameExpInteg = "../output/1Dcorr_0_10_cent_Integ.root";
    const TString fileNamePur = "../slurmOutput/apr12pur_all_24_09_26_processed.root";

    const std::vector<std::pair<int,TString> > ktArr{{1,"(200,400)"},{2,"(400,600)"},{3,"(600,800)"},{4,"(800,1000)"},{5,"(1000,1200)"},{6,"(1200,1400)"},{7,"(1400,1600)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.58,-0.35)"},{2,"(-0.35,-0.12)"},{3,"(-0.12,0.12)"},{4,"(0.12,0.35)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const int rebin = 1;

    TFile *inpFileData = TFile::Open(fileNameExp);
    TFile *inpFileDataInteg = TFile::Open(fileNameExpInteg);
    TFile *inpFilePur = TFile::Open(fileNamePur);

    TString otpFilePath = fileNameExp;
    otpFilePath.Insert(otpFilePath.Last('.'),"_Purity");
    TFile *otpFile = TFile::Open(otpFilePath,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    //TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    //canvKt->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &kt : ktArr)
    {
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignKt%d",kt.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgKt%d",kt.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatKt%d",kt.first));

        if (/* hData != nullptr && */ hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            hPur->SetName(TString::Format("hQinvPurKt%d",kt.first));
            hPur->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF(q_{inv})",kt.second.Data()));
            // add purity correction of hData

            hPur->Write();
        }
    }

    JJColor::CreatePrimaryWutGradient();

    //TCanvas *canvY = new TCanvas("canvY","",1600,900);
    //canvY->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &y : yArr)
    {
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignY%d",y.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgY%d",y.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatY%d",y.first));

        if (/* hData != nullptr && */ hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            hPur->SetName(TString::Format("hQinvPurY%d",y.first));
            hPur->SetTitle(TString::Format("y #in %s",y.second.Data()));
            // add purity correction of hData

            hPur->Write();
        }
    }

    JJColor::CreateSecondaryWutGradient();

    //TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    //canvPsi->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &psi : psiArr)
    {
        TH1D *hPurNum = inpFilePur->Get<TH1D>(TString::Format("hQinvSignPsi%d",psi.first));
        TH1D *hPurDen = inpFilePur->Get<TH1D>(TString::Format("hQinvBckgPsi%d",psi.first));
        TH1D *hData = inpFileData->Get<TH1D>(TString::Format("hQinvRatPsi%d",psi.first));

        if (/* hData != nullptr && */ hPurNum != nullptr && hPurDen != nullptr)
        {
            TH1D *hPur = new TH1D(*hPurNum);
            hPur->Divide(hPurDen);
            hPur->SetName(TString::Format("hQinvPsiPsi%d",psi.first));
            hPur->SetTitle(TString::Format("#phi - #Psi_{E.P.} #in %s",psi.second.Data()));
            // add purity correction of hData

            hPur->Write();
        }
    }
}