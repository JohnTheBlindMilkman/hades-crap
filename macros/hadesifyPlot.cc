#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "../Externals/Palettes.hxx"

TPaveText* addTextBox(const TString &text)
{
    TPaveText *pt = new TPaveText(0.2,0.8,0.5,0.9,"NB");
    pt->AddText(text);
    pt->SetTextAlign(kHAlignCenter+kVAlignCenter);
    pt->SetFillStyle(1);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.05);

    return pt;
}

void addLegend(TCanvas *c)
{

}

TCanvas* addFrame(TCanvas *c, TString title)
{
    const float topMargin = 0.075;
    const float whModifier = 1 + topMargin;
    TCanvas *newCanv = new TCanvas(c->GetName(),"",c->GetWw(),c->GetWh()*whModifier);
    newCanv->SetMargin(0.,0.,0.,0.);

    TPad *top = new TPad("top","",0.,1.-topMargin,1.,1.);
    TPad *bottom = new TPad("bottom","",0.,0.,1.,1.-topMargin);

    TPaveText *pt = new TPaveText(c->GetLeftMargin(),0.01,(1. - c->GetRightMargin()),0.9,"NB");
    pt->AddText(title);
    pt->SetTextFont(62);
    pt->SetTextAlign(kHAlignRight+kVAlignCenter);
    pt->SetFillStyle(1);
    pt->SetBorderSize(0);

    top->cd();
    pt->Draw();

    bottom->cd();
    c->DrawClonePad();

    newCanv->cd();
    top->Draw();
    bottom->Draw();

    return newCanv;
}

// Use this macro when you want to add a canvas with the "HADES Au+Au 0-10%" text and additional text box / legend / HADES preliminary
void hadesifyPlot()
{
    JJColor::CreateSecondaryWutGradient();

    TString inputFilePath = "../output/1Dcorr_0_10_cent.root";
    /* const TString canvasName1 = "canvKt";
    const TString canvasName2 = "canvY";
    const TString canvasName3 = "canvPsi"; */
    const TString canvasName4 = "canvKt";

    TFile *inpFile = TFile::Open(inputFilePath);
    /* TCanvas *canv1 = inpFile->Get<TCanvas>(canvasName1);
    TCanvas *canv2 = inpFile->Get<TCanvas>(canvasName2);
    TCanvas *canv3 = inpFile->Get<TCanvas>(canvasName3); */
    TCanvas *canv4 = inpFile->Get<TCanvas>(canvasName4);

    TFile *outFile = TFile::Open(inputFilePath.Insert(inputFilePath.Last('.'),"_hadesified"),"recreate");
    /* addFrame(canv1,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (50-60%)")->Write();
    addFrame(canv2,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (50-60%)")->Write();
    addFrame(canv3,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (50-60%)")->Write(); */
    TCanvas *newCanv = addFrame(canv4,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (0-10%)");
    TPaveText *tpt = addTextBox("#font[62]{HADES} #font[42]{Work in Progress}");
    newCanv->cd(2);
    tpt->Draw("same");

    newCanv->Write();
}