#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

double getNorm(TH1D *hInp, double xMin, double xMax)
{
    int nBins = 0;
    double val = 0., xVal;
    for (int i = 1; i < hInp->GetNbinsX(); i++)
    {
        xVal = hInp->GetBinCenter(i);
        if (xVal >= xMin && xVal <= xMax)
        {
            val += hInp->GetBinContent(i);
            nBins++;
        } 
    }
    if(nBins > 0)
        return val/nBins;
    else    
        return 0.;
}

void drawProton1DJJFM()
{
    gStyle->SetPalette(kPastel);

    const TString fileName = "/home/jedkol/Downloads/HADES/HADES-CrAP/slurmOutput/apr12ana_all_24_02_19.root";
    const TString outputFile = "/home/jedkol/Downloads/HADES/HADES-CrAP/output/1Dcorr_0_10_cent_forceEP_tmp.root";
    const int rebin = 2;

    float norm;
    TH1D *hSign = new TH1D("hSign","Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
    TH1D *hBckg = new TH1D("hBckg","Background of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
    TH1D *hRat = new TH1D("hRat","CF of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);

    TH2D *hDpDtSign = new TH2D("hDpDtSign","#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45);
    TH2D *hDpDtBckg = new TH2D("hDpDtBckg","#Delta#phi vs #Delta#theta distribution of backgound 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45);
    TH2D *hDpDtRat = new TH2D("hDpDtRat","#Delta#phi vs #Delta#theta distribution of signal/backgound ratio 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45);

    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);

    for(const int i : {1,2,3,4,5})
        for(const int j : {1,2,3})
            for(const int k : {1,2,3,4,5,6,7,8})
            {
                TH1D *hS = inpFile->Get<TH1D>(TString::Format("hQinvSign_%d%d%d",i,j,k));
                TH2D *hDpDtS = inpFile->Get<TH2D>(TString::Format("hDphiDthetaSign_%d%d%d",i,j,k));
                if (hS != nullptr && hDpDtS != nullptr)
                {
                    hSign->Add(hS);
                    hDpDtSign->Add(hDpDtS);
                }
                
                TH1D *hB = inpFile->Get<TH1D>(TString::Format("hQinvBckg_%d%d%d",i,j,k));
                TH2D *hDpDtB = inpFile->Get<TH2D>(TString::Format("hDphiDthetaBckg_%d%d%d",i,j,k));
                if (hB != nullptr && hDpDtB != nullptr)
                {
                    hBckg->Add(hB);
                    hDpDtBckg->Add(hDpDtB);
                }
            }

    hSign->Sumw2();
    hBckg->Sumw2();
    hRat = new TH1D(*hSign);

    hRat->Divide(hBckg);
    hRat->Sumw2();
    norm = getNorm(hRat,300,900);
    hRat->Rebin(rebin);
    norm *= rebin;
    hRat->Scale(1./norm);
    hRat->SetName("hQinvRat");
    hRat->SetTitle("CF of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})");
    hRat->SetMarkerStyle(20);

    hDpDtRat = new TH2D(*hDpDtSign);
    hDpDtRat->Divide(hDpDtBckg);
    hDpDtRat->SetName("hDphiDthetaRat");
    hDpDtRat->SetTitle("#Delta#phi vs #Delta#theta distribution of signal/backgound ratio 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]");
    
    TFile *outFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canv = new TCanvas("canv","",1600,900);
    hRat->Write();
    hSign->Write();
    hBckg->Write();
    hDpDtRat->Write();
    hRat->Draw("hist pe pmc plc");
    line->Draw("same");
    //canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canv->Write();

}