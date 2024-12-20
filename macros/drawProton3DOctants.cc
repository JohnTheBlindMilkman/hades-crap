#include <iostream>

#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"
#include "TPaveText.h"
#include "TPad.h"
#include "../Externals/Palettes.hxx"

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

void FillOctant(TH3D *hOctant, TH3D *hRat3D, int xFrom, int xTo, int yFrom, int yTo, int zFrom, int zTo)
{
    for (int xx = xFrom; xx <= xTo; ++xx)
        for (int yy = yFrom; yy <= yTo; ++yy)
            for (int zz = zFrom; zz <= zTo; ++zz)
                hOctant->SetBinContent(
                    hOctant->GetXaxis()->FindBin(abs(hRat3D->GetXaxis()->GetBinCenter(xx))),
                    hOctant->GetYaxis()->FindBin(abs(hRat3D->GetYaxis()->GetBinCenter(yy))),
                    hOctant->GetZaxis()->FindBin(abs(hRat3D->GetZaxis()->GetBinCenter(zz))),
                    hRat3D->GetBinContent(xx,yy,zz));
}

TH3D* CreateOctant(TH3D *hRat3D, TString id)
{
    int xMin = hRat3D->GetXaxis()->GetFirst();
    int xCent = hRat3D->GetXaxis()->FindBin(0.0);
    int xMax = hRat3D->GetXaxis()->GetLast();
    int yMin = hRat3D->GetYaxis()->GetFirst();
    int yCent = hRat3D->GetYaxis()->FindBin(0.0);
    int yMax = hRat3D->GetYaxis()->GetLast();
    int zMin = hRat3D->GetZaxis()->GetFirst();
    int zCent = hRat3D->GetZaxis()->FindBin(0.0);
    int zMax = hRat3D->GetZaxis()->GetLast();

    /* std::cout << xMin << "\t" << xCent << "\t" << xMax << "\n";
    std::cout << yMin << "\t" << yCent << "\t" << yMax << "\n";
    std::cout << zMin << "\t" << zCent << "\t" << zMax << "\n"; */

    TH3D *hOct = new TH3D(hRat3D->GetName() + TString("_") + id,hRat3D->GetTitle(),
    xMax-xCent+1,0,hRat3D->GetXaxis()->GetBinCenter(xMax),
    yMax-yCent+1,0,hRat3D->GetYaxis()->GetBinCenter(yMax),
    zMax-zCent+1,0,hRat3D->GetZaxis()->GetBinCenter(zMax));

    if (!id.CompareTo("ppp"))
    {
        FillOctant(hOct,hRat3D,xCent,xMax,yCent,yMax,zCent,zMax);
    }
    else if (!id.CompareTo("mpp"))
    {
        FillOctant(hOct,hRat3D,xMin,xCent,yCent,yMax,zCent,zMax);
    }
    else if (!id.CompareTo("pmp"))
    {
        FillOctant(hOct,hRat3D,xCent,xMax,yMin,yCent,zCent,zMax);
    }
    else if (!id.CompareTo("mmp"))
    {
        FillOctant(hOct,hRat3D,xMin,xCent,yMin,yCent,zCent,zMax);
    }
    else if (!id.CompareTo("ppm"))
    {
        FillOctant(hOct,hRat3D,xCent,xMax,yCent,yMax,zMin,zCent);
    }
    else if (!id.CompareTo("mpm"))
    {
        FillOctant(hOct,hRat3D,xMin,xCent,yCent,yMax,zMin,zCent);
    }
    else if (!id.CompareTo("pmm"))
    {
        FillOctant(hOct,hRat3D,xCent,xMax,yMin,yCent,zMin,zCent);
    }
    else if (!id.CompareTo("mmm"))
    {
        FillOctant(hOct,hRat3D,xMin,xCent,yMin,yCent,zMin,zCent);
    }

    std::cout << hOct->GetEntries() << "\n";
    hOct->Sumw2();

    return hOct;
}

void drawProjections(TH3D *hRat3D, TH3D *graphic3D, TString id, const int wbin, const int rebin)
{
    const std::vector<TString> sProj{"x","y","z"};
    const std::vector<TString> sProjName{"out","side","long"};

    TLine *line = new TLine(0,1,500,1);
    TCanvas *canvInteg = new TCanvas(TString::Format("canvInteg_%s",id.Data()),"",1800,600);
    canvInteg->Divide(3,1);

    TCanvas *canvWithGraphic = new TCanvas(TString::Format("canvWithGraphic_%s",id.Data()),"",1800,1200);
    canvWithGraphic->SetMargin(0.,0.,0.,0.);
    //canvWithGraphic->Divide(2,1);

    TPad *pLeft = new TPad("pLeft","",.0,.0,.6,1.);
    TPad *pRight = new TPad("pRight","",.6,.0,1.,1.);
    pRight->Divide(1,3);

    int binc,binmn,binmx;
    double norm;

    if (hRat3D != nullptr)
    {
        for (const int &i : {0,1,2})
        {
            binc = hRat3D->GetXaxis()->FindFixBin(0.0);
            binmx = binc + wbin;

            hRat3D->GetXaxis()->SetRange(1, (i == 0) ? hRat3D->GetNbinsX() : binmx);
            hRat3D->GetYaxis()->SetRange(1, (i == 1) ? hRat3D->GetNbinsY() : binmx);
            hRat3D->GetZaxis()->SetRange(1, (i == 2) ? hRat3D->GetNbinsZ() : binmx);

            TH1D *hRat = static_cast<TH1D*>(hRat3D->Project3D(sProj[i].Data()));
            norm = getNorm(hRat,350,500);
            hRat->Rebin(rebin);
            norm *= rebin;
            hRat->Scale(1./norm);
            hRat->GetYaxis()->SetRangeUser(0.,1.9);
            hRat->SetName(TString::Format("hQ%sRatInteg_%s",sProjName[i].Data(),id.Data()));
            hRat->SetTitle(TString::Format(";q_{%s} [MeV/c];CF(q_{%s})",sProjName[i].Data(),sProjName[i].Data()));
            hRat->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
            hRat->GetXaxis()->SetTitleSize(0.06);
            hRat->GetXaxis()->SetLabelSize(0.06);
            hRat->GetXaxis()->SetNdivisions(506);
            hRat->GetYaxis()->SetTitleSize(0.06);
            hRat->GetYaxis()->SetLabelSize(0.06);
            hRat->GetYaxis()->SetNdivisions(506);
            hRat->SetMarkerStyle(20);
            hRat->SetMarkerColor(JJColor::fWutAllColors[1]); // navy WUT

            hRat->Write();

            canvInteg->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);
            hRat->Draw("hist pe pmc plc");
            line->Draw("same");

            pRight->cd(i+1)->SetMargin(0.2,0.02,0.15,0.02);
            hRat->Draw("hist pe pmc plc");
            line->Draw("same");
        }

        pLeft->cd()->SetMargin(0.2,0.2,0.15,0.02);
        //graphic3D->SetMarkerColorAlpha(JJColor::fWutSecondaryColors[3],0.5);
        //graphic3D->SetFillColorAlpha(JJColor::fWutSecondaryColors[3],0.5);
        graphic3D->GetXaxis()->SetTitleSize(0.06);
        //graphic3D->GetXaxis()->SetLabelSize(0.06);
        graphic3D->GetXaxis()->SetNdivisions(206);
        graphic3D->GetYaxis()->SetTitleSize(0.06);
        //graphic3D->GetYaxis()->SetLabelSize(0.06);
        graphic3D->GetYaxis()->SetNdivisions(206);
        graphic3D->GetZaxis()->SetTitleSize(0.06);
        graphic3D->GetZaxis()->SetTitleOffset(1.2);
        //graphic3D->GetZaxis()->SetLabelSize(0.06);
        graphic3D->GetZaxis()->SetNdivisions(206);
        graphic3D->Draw("fb bb box3");

        canvWithGraphic->cd();
        pLeft->Draw();
        //canvWithGraphic->cd(2);
        pRight->Draw();

        canvInteg->Write();
        canvWithGraphic->Write();
        canvWithGraphic->SaveAs(TString::Format("../output/gif/canvWithGraphic_%s.gif",id.Data()));
    }
}

void drawProton3DOctants()
{
    gStyle->SetOptStat(0);
    JJColor::CreatePrimaryWutGradient();

    const TString fileName = "../slurmOutput/apr12ana_all_24_07_09_processed.root";
    const TString outputFile = "../output/3Dcorr_0_10_cent_Integ_Oct.root";
    const std::array<TString,2> signArr = {"p","m"};
    const int rebin = 1;
    const int wbin = 1; 

    float norm;
    int binc, binmn, binmx;
    TString octantHash = "";

    TLine *line = new TLine(-500,1,500,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    TH3D *hSign = inpFile->Get<TH3D>("hQoslSignInteg");
    TH3D *hBckg = inpFile->Get<TH3D>("hQoslBckgInteg");

    for (const int &x : {0,1})
        for (const int &y : {0,1})
            for (const int &z : {0,1})
            {
                octantHash = signArr[x] + signArr[y] + signArr[z];
                TH3D *hTmpSign = CreateOctant(hSign,octantHash);
                hTmpSign->Write();
                TH3D *hTmpBckg = CreateOctant(hBckg,octantHash);
                hTmpBckg->Write();

                TH3D *hRat3D = new TH3D(*hTmpSign);
                hRat3D->Divide(hTmpBckg);
                hRat3D->SetName("hQoslRatInteg_"+octantHash);
                hRat3D->Write();

                TH3D *hGraphic = new TH3D("hGraphic",";q_{out};q_{side};q_{long}",2,-500,500,2,-500,500,2,-500,500);
                //hGraphic->SetMarkerColor(JJColor::fWutSecondaryColors[3]);
                hGraphic->SetFillColor(JJColor::fWutMainColors[1]);
                hGraphic->SetBinContent(x+1,y+1,z+1,1);
                drawProjections(hRat3D,hGraphic,octantHash,wbin,rebin);
            }
}