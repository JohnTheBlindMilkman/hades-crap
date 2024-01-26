#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"

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

    const TString fileName = "/home/jedkol/lxpool/hades-crap/femtoOutFile.root";
    const TString outputFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_tmp.root";
    const int rebin = 2;

    float norm;
    TH1D *hSign = new TH1D("hSign","Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
    TH1D *hBckg = new TH1D("hBckg","Background of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
    TH1D *hRat = new TH1D("hRat","CF of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);

    TFile *inpFile = TFile::Open(fileName);

    for(const int i : {1,2,3,4,5})
        for(const int j : {1})
            for(const int k : {1,2,3,4})
            {
                TH1D *hS = inpFile->Get<TH1D>(TString::Format("hQinvSign_%d%d%d",i,j,k));
                if (hS != nullptr)
                {
                    hSign->Add(hS);
                }
                
                TH1D *hB = inpFile->Get<TH1D>(TString::Format("hQinvBckg_%d%d%d",i,j,k));
                if (hB != nullptr)
                {
                    hBckg->Add(hB);
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
    
    TFile *outFile = TFile::Open(outputFile,"RECREATE");

    TCanvas *canv = new TCanvas("canv","",1600,900);
    hRat->Write();
    hRat->Draw("hist pe pmc plc");
    //canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    canv->Write();
}