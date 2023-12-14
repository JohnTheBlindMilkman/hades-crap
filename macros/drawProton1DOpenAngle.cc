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
        return (double)(val/nBins);
    else    
        return 0.;
}

void drawProton1DOpenAngle()
{
    gStyle->SetPalette(kPastel);

    const TString fileName = "/home/jedkol/lxpool/hades-crap/testOutFile2.root";
    const std::vector<int> alphaList = {0,1,2,3,4,5};
    const int alphaLen = static_cast<int>(alphaList.size());
    const int rebin = 2;

    float norm;
    TH1D *hSign[alphaLen],*hBckg[alphaLen],*hRat[alphaLen];

    TFile *inpFile = TFile::Open(fileName);

    for (int i = 0; i < alphaLen; i++)
    {
        hSign[i] = inpFile->Get<TH1D>(TString::Format("hQinvSign%d",i));
        hSign[i]->Sumw2();
        hBckg[i] = inpFile->Get<TH1D>(TString::Format("hQinvBckg%d",i));
        hBckg[i]->Sumw2();
        hRat[i] = new TH1D(*hSign[i]);

        hRat[i]->Divide(hBckg[i]);
        hRat[i]->Sumw2();
        norm = getNorm(hRat[i],300,900);
        hRat[i]->Rebin(rebin);
        norm *= rebin;
        hRat[i]->Scale(1./norm);
        hRat[i]->SetName(TString::Format("hQinvRat%d",i));
        hRat[i]->SetTitle(TString::Format("#alpha = %d#circ;q_{inv} [MeV/c];C(q_{inv})",alphaList.at(i)));
        hRat[i]->SetMarkerStyle(20);
    }
    
    TCanvas *canv = new TCanvas("canv","",1600,900);
    for (auto hist : hRat)
    {
        if (hist == hRat[0])
        {
            hist->GetXaxis()->SetRangeUser(0,500);
            hist->GetYaxis()->SetRangeUser(0.,2.);
            hist->Draw("hist pe pmc plc");
        }
        else
            hist->Draw("hist pe pmc plc same");
    }
    canv->BuildLegend(0.2,0.2,0.5,0.5,"","p");
}