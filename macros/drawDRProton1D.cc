#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

void SetErrors(TH1D *hout, TH1D *hNum, TH1D *hDen)
{
    const int iterMax = hout->GetNbinsX();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= iterMax; i++)
    {
        vErr = 0;
        vNum = hNum->GetBinContent(i);
        eNum = hNum->GetBinError(i);
        vDen = hDen->GetBinContent(i);
        eDen = hDen->GetBinError(i);

        // propagacja błędów dla funcji num/den z uwzględnieniem korelacji pomiędzy składowymi
        if (fabs(vDen) > std::numeric_limits<double>::epsilon())
            vErr = TMath::Sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
}

void drawDRProton1D()
{
    gStyle->SetPalette(kPastel);

    const TString dataFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent.root";
    const TString simFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_HGeant.root";
    const TString otpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR.root";
    const std::vector<std::pair<TString,TString> > kTbins = {{"150","450"},{"450","750"},{"750","1050"},{"1050","1350"},{"1350","1650"}};
    const TString histNameBase = "hQinvRat";
    const std::size_t kTlen = kTbins.size();
    const int rebin = 1;

    std::array<TH1D*,5> hAna, hSim, hRatio;
    TFile *fInpData,*fInpSim,*fOtpFile;

    fInpData = TFile::Open(dataFile);
    fInpSim = TFile::Open(simFile);

    for (const auto &iter : {0,1,2,3,4})
    {
        hAna.at(iter) = fInpData->Get<TH1D>(TString::Format("%s%d",histNameBase.Data(),iter));
        hSim.at(iter) = fInpSim->Get<TH1D>(TString::Format("%s%d",histNameBase.Data(),iter));

        hRatio.at(iter) = new TH1D(*hAna.at(iter));
        hRatio.at(iter)->Divide(hSim.at(iter));
        hRatio.at(iter)->SetName(TString::Format("hQinvRatDR%d",iter));
        SetErrors(hRatio.at(iter),hAna.at(iter),hSim.at(iter));

        if (rebin > 1)
        {
            hRatio.at(iter)->Rebin(rebin);
            hRatio.at(iter)->Scale(1./rebin);
        }
        
        hRatio.at(iter)->SetTitle(TString::Format("k_{T} #in (%s,%s) [MeV/c];q_{inv} [MeV/c];C(q_{inv})",std::get<0>(kTbins.at(iter)).Data(),std::get<1>(kTbins.at(iter)).Data()));
        hRatio.at(iter)->SetMarkerStyle(20);
    }

    fOtpFile = TFile::Open(otpFile,"RECREATE");

    TCanvas *canv = new TCanvas("canv","",1600,900);
    for (auto &hist : hRatio)
    {
        hist->Write();
        if (hist == hRatio.at(0))
        {
            hist->GetYaxis()->SetRangeUser(0,2);
            hist->Draw("hist pe pmc plc");
        }
        else
            hist->Draw("hist pe pmc plc same");
    }
    canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canv->Write();
}