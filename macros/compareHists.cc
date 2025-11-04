#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TString.h"

#include "../Externals/Palettes.hxx"
#include "MacroUtils.hxx"

void compareHists()
{
    /* const TString fileName1 = "../output/3Dcorr_0_10_cent_Integ.root";
    const TString fileName2 = "../output/3Dcorr_0_10_cent_Integ_tmp.root";
    const TString fileName3 = "../output/3Dcorr_0_10_cent_Integ_tmp2.root";
    const TString fileName4 = "../output/3Dcorr_0_10_cent_Integ_tmp3.root";
    const TString fileNameOut = "../output/3Dcorr_0_10_cent_Integ_Comp.root";
    const std::vector<TString> histNames = {"hQoutRatInteg","hQsideRatInteg","hQlongRatInteg"};
    const std::vector<TString> titleNames = {"Old Implementation","New Implementation","MDC min HV #pm 5V","FCH = 75%"};
    const std::vector<TString> sProjName{"out","side","long"};

    TFile *inpFile1 = TFile::Open(fileName1);
    TFile *inpFile2 = TFile::Open(fileName2);
    TFile *inpFile3 = TFile::Open(fileName3);
    TFile *inpFile4 = TFile::Open(fileName4);
    TFile *outFile = TFile::Open(fileNameOut,"RECREATE");

    TLine *line = new TLine(0,1,480,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TCanvas *c = new TCanvas("c","",1800,600);
    c->Divide(3,1);

    JJColor::CreatePrimaryWutGradient();

    for (const int &i : {0,1,2})
    {
        inpFile1->cd();
        TH1D *hInp1 = inpFile1->Get<TH1D>(histNames.at(i));
        hInp1->SetMarkerColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
        hInp1->GetYaxis()->SetRangeUser(0.6,1.4);
        hInp1->GetXaxis()->SetRangeUser(-1,480);
        hInp1->SetName(TString::Format("%s_1",histNames.at(i).Data()));
        hInp1->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[0].Data(),sProjName[i].Data(),sProjName[i].Data()));
        hInp1->GetXaxis()->SetTitleSize(0.06);
        hInp1->GetXaxis()->SetLabelSize(0.06);
        hInp1->GetXaxis()->SetNdivisions(506);
        hInp1->GetYaxis()->SetTitleSize(0.06);
        hInp1->GetYaxis()->SetLabelSize(0.06);
        hInp1->GetYaxis()->SetNdivisions(506);
        outFile->cd();
        hInp1->Write();

        inpFile2->cd();
        TH1D *hInp2 = inpFile2->Get<TH1D>(histNames.at(i));
        hInp2->SetMarkerColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
        hInp2->SetName(TString::Format("%s_2",histNames.at(i).Data()));
        hInp2->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[1].Data(),sProjName[i].Data(),sProjName[i].Data()));
        outFile->cd();
        hInp2->Write();

        inpFile3->cd();
        TH1D *hInp3 = inpFile3->Get<TH1D>(histNames.at(i));
        hInp3->SetMarkerColor(JJColor::fWutSecondaryColors[3]); // secondary fire red WUT
        hInp3->SetName(TString::Format("%s_3",histNames.at(i).Data()));
        hInp3->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[2].Data(),sProjName[i].Data(),sProjName[i].Data()));
        outFile->cd();
        hInp3->Write();

        inpFile4->cd();
        TH1D *hInp4 = inpFile4->Get<TH1D>(histNames.at(i));
        hInp4->SetMarkerColor(JJColor::fWutSecondaryColors[5]); // secondary violet WUT
        hInp4->SetName(TString::Format("%s_4",histNames.at(i).Data()));
        hInp4->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[3].Data(),sProjName[i].Data(),sProjName[i].Data()));
        outFile->cd();
        hInp4->Write();

        TVirtualPad *tvp = c->cd(i+1);
        tvp->SetMargin(0.2,0.02,0.15,0.02);
        
        hInp1->Draw("hist pe");
        hInp2->Draw("hist pe same");
        hInp3->Draw("hist pe same");
        hInp4->Draw("hist pe same");
        line->Draw("same");

        if (i == 1)
            tvp->BuildLegend(.3,.21,.3,.21,"","p");
    }
    
    outFile->cd();
    c->Write(); */

    const TString fileName1 = "../output/1Dcorr_0_10_cent_Purity_MomRes_old.root";
    const TString fileName2 = "../output/1Dcorr_0_10_cent_Purity_MomRes.root";
    const TString fileName3 = "../output/1Dcorr_0_10_cent_Integ.root";
    // const TString fileName4 = "../output/1Dcorr_0_10_cent_Integ_tmp3.root";
    // const TString fileName5 = "../output/1Dcorr_0_10_cent_Integ_tmp4.root";
    const TString fileNameOut = "../output/1Dcorr_0_10_cent_Purity_MomRes_Comp_Ratio.root";
    const std::vector<TString> histNames = {"hQinvRatInteg","hQinvRatInteg","hQinvRatInteg"};
    const std::vector<TString> titleNames = 
    {
        "#frac{Mom. Res. Correction with 1/p}{Mom. Res. Correction with p}",
        "#frac{Mom. Res. Correction with 1/p}{Uncorrected Function}",
        "#frac{Mom. Res. Correction with p}{Uncorrected Function}"
    };

    TFile *inpFile1 = TFile::Open(fileName1);
    TFile *inpFile2 = TFile::Open(fileName2);
    TFile *inpFile3 = TFile::Open(fileName3);
    // TFile *inpFile4 = TFile::Open(fileName4);
    // TFile *inpFile5 = TFile::Open(fileName5);
    TFile *outFile = TFile::Open(fileNameOut,"RECREATE");

    TLine *line = new TLine(0,1,480,1);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kGray);

    TCanvas *c = new TCanvas("c","",1800,600);
    c->SetMargin(0.2,0.02,0.15,0.02);

    JJColor::CreatePrimaryWutGradient();

    inpFile1->cd();
    TH1D *hInp1 = inpFile1->Get<TH1D>(histNames.at(0));
    hInp1->SetMarkerColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
    hInp1->SetLineColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
    hInp1->SetTitle(titleNames.at(0));
    /* hInp1->GetYaxis()->SetRangeUser(0.6,1.4);
    hInp1->GetXaxis()->SetRangeUser(-1,480);
    hInp1->SetName(TString::Format("%s_1",histNames.at(i).Data()));
    hInp1->SetTitle(TString::Format("%s;q_{%s} [MeV/c];CF(q_{%s})",titleNames[0].Data(),sProjName[i].Data(),sProjName[i].Data()));
    hInp1->GetXaxis()->SetTitleSize(0.06);
    hInp1->GetXaxis()->SetLabelSize(0.06);
    hInp1->GetXaxis()->SetNdivisions(506);
    hInp1->GetYaxis()->SetTitleSize(0.06);
    hInp1->GetYaxis()->SetLabelSize(0.06);
    hInp1->GetYaxis()->SetNdivisions(506); */
    outFile->cd();
    hInp1->Write();

    inpFile2->cd();
    TH1D *hInp2 = inpFile2->Get<TH1D>(histNames.at(1));
    hInp2->SetMarkerColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
    hInp2->SetLineColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
    hInp2->SetTitle(titleNames.at(1));
    outFile->cd();
    hInp2->Write();

    inpFile3->cd();
    TH1D *hInp3 = inpFile3->Get<TH1D>(histNames.at(2));
    hInp3->SetMarkerColor(JJColor::fWutSecondaryColors[3]); // secondary red WUT
    hInp3->SetLineColor(JJColor::fWutSecondaryColors[3]); // secondary red WUT
    hInp3->SetTitle(titleNames.at(2));
    outFile->cd();
    hInp3->Write();

    /* inpFile4->cd();
    TH1D *hInp4 = inpFile4->Get<TH1D>(histNames.at(3));
    hInp4->SetMarkerColor(JJColor::fWutSecondaryColors[4]);
    hInp4->SetLineColor(JJColor::fWutSecondaryColors[4]);
    hInp4->SetTitle(titleNames.at(3));
    outFile->cd();
    hInp4->Write();

    inpFile5->cd();
    TH1D *hInp5 = inpFile5->Get<TH1D>(histNames.at(4));
    hInp5->SetMarkerColor(JJColor::fWutSecondaryColors[5]);
    hInp5->SetLineColor(JJColor::fWutSecondaryColors[5]);
    hInp5->SetTitle(titleNames.at(4));
    outFile->cd();
    hInp5->Write(); */

    TH1D *hRatio1 = static_cast<TH1D*>(hInp2->Clone("hRatio1"));
    hRatio1->Divide(hInp1);
    JJUtils::Generic::SetErrorsDivide(hRatio1,hInp2,hInp1);
    hRatio1->SetMarkerColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
    hRatio1->SetLineColor(JJColor::fWutSecondaryColors[1]); // secondary blue WUT
    hRatio1->SetTitle(titleNames.at(0));

    TH1D *hRatio2 = static_cast<TH1D*>(hInp2->Clone("hRatio2"));
    hRatio2->Divide(hInp3);
    JJUtils::Generic::SetErrorsDivide(hRatio2,hInp1,hInp3);
    hRatio2->SetMarkerColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
    hRatio2->SetLineColor(JJColor::fWutSecondaryColors[2]); // secondary gold WUT
    hRatio2->SetTitle(titleNames.at(1));

    TH1D *hRatio3 = static_cast<TH1D*>(hInp1->Clone("hRatio3"));
    hRatio3->Divide(hInp3);
    JJUtils::Generic::SetErrorsDivide(hRatio3,hInp1,hInp3);
    hRatio3->SetMarkerColor(JJColor::fWutSecondaryColors[3]); // secondary red WUT
    hRatio3->SetLineColor(JJColor::fWutSecondaryColors[3]); // secondary red WUT
    hRatio3->SetTitle(titleNames.at(2));

    hRatio1->Draw("pe");
    hRatio2->Draw("pe same");
    hRatio3->Draw("pe same");
    // hInp1->Draw("pe");
    // hInp2->Draw("pe same");
    // hInp3->Draw("pe same");
    // hInp4->Draw("pe same");
    // hInp5->Draw("pe same");

    c->BuildLegend(0.2,0.2,0.5,0.5,"","p");

    c->Write();
}
