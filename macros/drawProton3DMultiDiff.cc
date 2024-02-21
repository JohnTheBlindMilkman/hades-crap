#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "Palettes.hxx"

double getNorm(const TH1D *hInp, double xMin, double xMax)
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

void prepareGraph(TH1D* hist, int col)
{
    hist->SetMarkerColor(col);
    hist->SetMarkerStyle(20);
    hist->SetLineColor(col);

    //hist->GetXaxis()->SetLabelSize();
    //hist->GetXaxis()->SetLabelOffset();
    //hist->GetXaxis()->SetTitleSize();
    //hist->GetXaxis()->SetTitleOffset();

    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetLabelOffset(0.005);
    hist->GetYaxis()->SetLabelFont(62);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetYaxis()->SetTitleFont(62);
    hist->GetYaxis()->SetNdivisions(510);
}

void drawProton3DMultiDiff()
{
    const TString fileName = "/home/jedkol/lxpool/hades-crap/slurmOutput/apr12ana_all_24_02_20_dp8dt2.root";
    const TString outputFile = "/home/jedkol/lxpool/hades-crap/output/3Dcorr_0_10_cent_dp8dp2.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(150,450)"},{2,"(450,750)"},{3,"(750,1050)"},{4,"(1050,1350)"},{5,"(1350,1650)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.75,-0.25)"},{2,"(-0.25,0.25)"},{3,"(0.25,0.75)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};
    const int rebin = 1;
    const int wbin = 0; 

    float norm;
    int binc, binmn, binmx;
    std::vector<std::vector<std::vector<TH3D*> > > 
    hSign(ktArr.size(),std::vector<std::vector<TH3D*> >(yArr.size(),std::vector<TH3D*>(psiArr.size(),nullptr))), 
    hBckg(ktArr.size(),std::vector<std::vector<TH3D*> >(yArr.size(),std::vector<TH3D*>(psiArr.size(),nullptr)));

    TLine *line = new TLine(0,1,3000,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);

    for (const auto &kt : ktArr)
        for (const auto &y : yArr)
            for (const auto &psi : psiArr)
            {
                hSign[kt.first-1][y.first-1][psi.first-1] = inpFile->Get<TH3D>(TString::Format("hQoslSign_%d%d%d",kt.first,y.first,psi.first));
                if (hSign[kt.first-1][y.first-1][psi.first-1] != nullptr)
                    hSign[kt.first-1][y.first-1][psi.first-1]->Sumw2();

                hBckg[kt.first-1][y.first-1][psi.first-1] = inpFile->Get<TH3D>(TString::Format("hQoslBckg_%d%d%d",kt.first,y.first,psi.first));
                if (hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                    hBckg[kt.first-1][y.first-1][psi.first-1]->Sumw2();
            }

    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->Divide(3,1);
    std::vector<TH3D*> hSignKt(ktArr.size(),nullptr), hBckgKt(ktArr.size(),nullptr);
    std::vector<std::vector<TH1D*> > hRatPPKt(ktArr.size(),std::vector<TH1D*>(sProj.size(),nullptr)),
                                     hNumPPKt(ktArr.size(),std::vector<TH1D*>(sProj.size(),nullptr)), 
                                     hDenPPKt(ktArr.size(),std::vector<TH1D*>(sProj.size(),nullptr));
    for (const auto &kt : ktArr)
    {
        for (const auto &y : yArr)
            for (const auto &psi : psiArr)
            {
                if (hSignKt[kt.first-1] == nullptr && hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckgKt[kt.first-1] == nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hSignKt[kt.first-1] = new TH3D(*hSign[kt.first-1][y.first-1][psi.first-1]);
                    hSignKt[kt.first-1]->SetName(TString::Format("hQoslSignKt%d",kt.first));
                    hBckgKt[kt.first-1] = new TH3D(*hBckg[kt.first-1][y.first-1][psi.first-1]);
                    hBckgKt[kt.first-1]->SetName(TString::Format("hQoslBckgKt%d",kt.first));
                }
                else if(hSign[kt.first-1][y.first-1][psi.first-1] != nullptr && hBckg[kt.first-1][y.first-1][psi.first-1] != nullptr)
                {
                    hSignKt[kt.first-1]->Add(hSign[kt.first-1][y.first-1][psi.first-1]);
                    hBckgKt[kt.first-1]->Add(hBckg[kt.first-1][y.first-1][psi.first-1]);
                }
            }

        if (hSignKt[kt.first-1] != nullptr && hBckgKt[kt.first-1] != nullptr)
        {
            for (const int &i : {0,1,2})
            {
                binc = hSignKt[kt.first-1]->GetXaxis()->FindFixBin(0.0);
                binmx = binc + wbin;
                binmn = binc - wbin;

                hSignKt[kt.first-1]->GetXaxis()->SetRange(1, (i == 0) ? hSignKt[kt.first-1]->GetNbinsX() : binmx);
                hSignKt[kt.first-1]->GetYaxis()->SetRange(1, (i == 1) ? hSignKt[kt.first-1]->GetNbinsY() : binmx);
                hSignKt[kt.first-1]->GetZaxis()->SetRange(1, (i == 2) ? hSignKt[kt.first-1]->GetNbinsZ() : binmx);
                hBckgKt[kt.first-1]->GetXaxis()->SetRange(1, (i == 0) ? hBckgKt[kt.first-1]->GetNbinsX() : binmx);
                hBckgKt[kt.first-1]->GetYaxis()->SetRange(1, (i == 1) ? hBckgKt[kt.first-1]->GetNbinsY() : binmx);
                hBckgKt[kt.first-1]->GetZaxis()->SetRange(1, (i == 2) ? hBckgKt[kt.first-1]->GetNbinsZ() : binmx);

                hNumPPKt[kt.first-1][i] = static_cast<TH1D*>(hSignKt[kt.first-1]->Project3D(sProj[i].Data()));
                hDenPPKt[kt.first-1][i] = static_cast<TH1D*>(hBckgKt[kt.first-1]->Project3D(sProj[i].Data()));

                hRatPPKt[kt.first-1][i] = new TH1D(*hNumPPKt[kt.first-1][i]);
                hRatPPKt[kt.first-1][i]->Divide(hDenPPKt[kt.first-1][i]);
                norm = getNorm(hRatPPKt[kt.first-1][i],300,900);
                hRatPPKt[kt.first-1][i]->Rebin(rebin);
                norm *= rebin;
                hRatPPKt[kt.first-1][i]->Scale(1./norm);
            }
            hRatKt[kt.first-1] = new TH1D(*hSignKt[kt.first-1]);
            hRatKt[kt.first-1]->Divide(hBckgKt[kt.first-1]);
            norm = getNorm(hRatKt[kt.first-1],300,900);
            hRatKt[kt.first-1]->Rebin(rebin);
            norm *= rebin;
            hRatKt[kt.first-1]->Scale(1./norm);
            hRatKt[kt.first-1]->SetName(TString::Format("hQinvRatKt%d",kt.first));
            hRatKt[kt.first-1]->SetTitle(TString::Format("k_{T} #in %s",kt.second.Data()));
            hRatKt[kt.first-1]->SetMarkerStyle(20);
            hRatKt[kt.first-1]->SetMarkerColor(JJColor::fWutAllColors[kt.first-1]);

            hRatKt[kt.first-1]->Write();
            hSignKt[kt.first-1]->Write();
            hBckgKt[kt.first-1]->Write();
            if (kt.first -1 == 0)
                hRatKt[kt.first-1]->Draw("hist pe pmc plc");
            else
                hRatKt[kt.first-1]->Draw("same hist pe pmc plc");
        }
    }
    
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line->Draw("same");
    canvKt->Write();
}