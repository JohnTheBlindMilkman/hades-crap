#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TPaveText.h"

void addTextBox(TCanvas *c, TString text)
{

}

void addLegend(TCanvas *c)
{

}

TCanvas* addFrame(TCanvas *c, TString title)
{
    const float topMargin = 0.1;
    const float whModifier = 1 + topMargin;
    TCanvas *newCanv = new TCanvas(c->GetName(),"",c->GetWw(),c->GetWh()*whModifier);
    newCanv->SetMargin(0.,0.,0.,0.);

    TPad *top = new TPad("top","",0.,1.-topMargin,1.,1.);
    TPad *bottom = new TPad("bottom","",0.,0.,1.,1.-topMargin);

    TPaveText *pt = new TPaveText(0.01,0.,0.99,1.,"NB");
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
    TString inputFilePath = "../output/3Dcorr_0_10_cent_dp8dt2_Integ.root";
    const TString canvasName1 = "canvInteg";
    //const TString canvasName2 = "canv2DInteg";

    TFile *inpFile = TFile::Open(inputFilePath);
    TCanvas *canv1 = inpFile->Get<TCanvas>(canvasName1);
    //TCanvas *canv2 = inpFile->Get<TCanvas>(canvasName2);

    TFile *outFile = TFile::Open(inputFilePath.Insert(inputFilePath.Last('.'),"_hadesified"),"recreate");
    addFrame(canv1,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (0-10%)")->Write();
    //addFrame(canv2,"Au+Au #sqrt{s_{NN}} = 2.4 GeV (0-10%)")->Write();
}