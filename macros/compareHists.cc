#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "Palettes.hxx"

void compareHists()
{
    const TString fileName1 = "../output/3Dcorr_0_10_cent_Integ_tmp.root";
    const TString fileName2 = "../output/3Dcorr_0_10_cent_Integ.root";
    const TString fileNameOut = "../output/3Dcorr_0_10_cent_Integ_Comp.root";
    const std::vector<TString> histNames = {"hQoutRatInteg","hQsideRatInteg","hQlongRatInteg"};
    const std::vector<TString> titleNames = {"Uncorrected","Corrected"};
    const std::vector<TString> sProjName{"out","side","long"};

    std::array<std::array<TH1D*,3>,2> histArr;

    TFile *inpFile1 = TFile::Open(fileName1);
    TFile *inpFile2 = TFile::Open(fileName2);
    TFile *outFile = TFile::Open(fileNameOut,"RECREATE");

    TCanvas *c = new TCanvas("c","",1800,600);
    c->Divide(3,1);

    JJColor::CreatePrimaryWutGradient();

    for (const int &i : {0,1,2})
    {
        inpFile1->cd();
        TH1D *hInp1 = inpFile1->Get<TH1D>(histNames.at(i));
        hInp1->SetMarkerColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
        hInp1->SetName(TString::Format("%s_1",histNames.at(i).Data()));
        hInp1->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[i].Data(),sProjName[i].Data(),sProjName[i].Data()));
        outFile->cd();
        hInp1->Write();

        inpFile2->cd();
        TH1D *hInp2 = inpFile2->Get<TH1D>(histNames.at(i));
        hInp2->SetMarkerColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
        hInp2->SetName(TString::Format("%s_2",histNames.at(i).Data()));
        hInp2->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[i].Data(),sProjName[i].Data(),sProjName[i].Data()));
        outFile->cd();
        hInp2->Write();

        TVirtualPad *tvp = c->cd(i+1);

        hInp1->Draw("hist pe");
        hInp2->Draw("hist pe same");

        if (i == 1)
            tvp->BuildLegend(.3,.21,.3,.21,"","p");
    }
    
    outFile->cd();
    c->Write();
}