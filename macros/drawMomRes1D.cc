#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TLegend.h"
#include "../Externals/Palettes.hxx"
#include "TLine.h"

void SetErrors(TH1D *hout, const TH1D *hNum, TH1D *hDen)
{
    const int iterMax = hout->GetNbinsX();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= iterMax; i++)
    {
        vErr = 0;
        vNum = hNum->GetBinContent(i);
        eNum = hNum->GetBinError(i);
        vDen = hDen->GetBinContent(hDen->FindBin(hNum->GetBinCenter(i)));
        eDen = hDen->GetBinError(hDen->FindBin(hNum->GetBinCenter(i)));

        // propagation of uncertainty for a function num/den with onclusion of the correlation between the constituents
        if (fabs(vDen) > 0.)
            vErr = std::sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
        hout->SetBinError(i,vErr);
    }
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

    hist->GetXaxis()->SetTitleOffset();
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
}

void drawMomRes1D()
{
    const TString fileName = "../slurmOutput/apr12momres_all_25_01_13_processed.root";
    const TString outputFile = "../output/1Dmomres_0_10_cent.root";
    const std::vector<std::pair<int,TString> > ktArr{{1,"(200,400)"},{2,"(400,600)"},{3,"(600,800)"},{4,"(800,1000)"},{5,"(1000,1200)"},{6,"(1200,1400)"},{7,"(1400,1600)"}};
    const std::vector<std::pair<int,TString> > yArr{{1,"(-0.58,-0.35)"},{2,"(-0.35,-0.12)"},{3,"(-0.12,0.12)"},{4,"(0.12,0.35)"}};
    const std::vector<std::pair<int,TString> > psiArr{{1,"(-202.5,-157.5)"},{2,"(-157.5,-112.5)"},{3,"(-112.5,-67.5)"},{4,"(-67.5,-22.5)"},{5,"(-22.5,22.5)"},{6,"(22.5,67.5)"},{7,"(67.5,112.5)"},{8,"(112.5,157.5)"}};

    TFile *inpFile = TFile::Open(fileName);
    TFile *otpFile = TFile::Open(outputFile,"RECREATE");

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    canvKt->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &kt : ktArr)
    {
        TH1D *hSignGK = inpFile->Get<TH1D>(TString::Format("hQinvSignGKKt%d",kt.first));
        TH1D *hBckgGK = inpFile->Get<TH1D>(TString::Format("hQinvBckgGKKt%d",kt.first));
        TH1D *hSignPC = inpFile->Get<TH1D>(TString::Format("hQinvSignPCKt%d",kt.first));
        TH1D *hBckgPC = inpFile->Get<TH1D>(TString::Format("hQinvBckgPCKt%d",kt.first));
        if (hSignGK != nullptr && hBckgGK != nullptr && hSignPC != nullptr && hBckgPC != nullptr)
        {
            TH1D *hRatGK = new TH1D(*hSignGK);
            hRatGK->Divide(hBckgGK);
            hRatGK->SetName(TString::Format("hQinvRatGKKt%d",kt.first));
            hRatGK->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})",kt.second.Data()));
            prepareGraph(hRatGK,JJColor::fWutMainColors[1]);

            TH1D *hRatPC = new TH1D(*hSignPC);
            hRatPC->Divide(hBckgPC);
            hRatPC->SetName(TString::Format("hQinvRatPCKt%d",kt.first));
            hRatPC->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF_{smeared}(q_{inv})",kt.second.Data()));
            prepareGraph(hRatPC,JJColor::fWutMainColors[1]);

            TH1D *hRat = new TH1D(*hRatGK);
            hRat->Divide(hRatPC);
            SetErrors(hRat,hRatGK,hRatPC);
            hRat->SetName(TString::Format("hQinvMomResKt%d",kt.first));
            hRat->SetTitle(TString::Format("k_{T} #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})/CF_{smeared}(q_{inv})",kt.second.Data()));
            prepareGraph(hRat,JJColor::fWutMainColors[1]);

            hRat->Write();
            hRatGK->Write();
            hRatPC->Write();
        }
    }
    canvKt->Write();

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    canvY->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &y : yArr)
    {
        TH1D *hSignGK = inpFile->Get<TH1D>(TString::Format("hQinvSignGKY%d",y.first));
        TH1D *hBckgGK = inpFile->Get<TH1D>(TString::Format("hQinvBckgGKY%d",y.first));
        TH1D *hSignPC = inpFile->Get<TH1D>(TString::Format("hQinvSignPCY%d",y.first));
        TH1D *hBckgPC = inpFile->Get<TH1D>(TString::Format("hQinvBckgPCY%d",y.first));
        if (hSignGK != nullptr && hBckgGK != nullptr && hSignPC != nullptr && hBckgPC != nullptr)
        {
            TH1D *hRatGK = new TH1D(*hSignGK);
            hRatGK->Divide(hBckgGK);
            hRatGK->SetName(TString::Format("hQinvRatGKY%d",y.first));
            hRatGK->SetTitle(TString::Format("y #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})",y.second.Data()));
            prepareGraph(hRatGK,JJColor::fWutMainColors[1]);

            TH1D *hRatPC = new TH1D(*hSignPC);
            hRatPC->Divide(hBckgPC);
            hRatPC->SetName(TString::Format("hQinvRatPCY%d",y.first));
            hRatPC->SetTitle(TString::Format("y #in %s;q_{inv} [MeV/c];CF_{smeared}(q_{inv})",y.second.Data()));
            prepareGraph(hRatPC,JJColor::fWutMainColors[1]);

            TH1D *hRat = new TH1D(*hRatGK);
            hRat->Divide(hRatPC);
            SetErrors(hRat,hRatGK,hRatPC);
            hRat->SetName(TString::Format("hQinvMomResY%d",y.first));
            hRat->SetTitle(TString::Format("y #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})/CF_{smeared}(q_{inv})",y.second.Data()));
            prepareGraph(hRat,JJColor::fWutMainColors[1]);

            hRat->Write();
            hRatGK->Write();
            hRatPC->Write();
        }
    }
    canvY->Write();

    JJColor::CreateSecondaryWutGradient();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    canvPsi->SetMargin(0.2,0.02,0.15,0.02);
    for (const auto &psi : psiArr)
    {
        TH1D *hSignGK = inpFile->Get<TH1D>(TString::Format("hQinvSignGKPsi%d",psi.first));
        TH1D *hBckgGK = inpFile->Get<TH1D>(TString::Format("hQinvBckgGKPsi%d",psi.first));
        TH1D *hSignPC = inpFile->Get<TH1D>(TString::Format("hQinvSignPCPsi%d",psi.first));
        TH1D *hBckgPC = inpFile->Get<TH1D>(TString::Format("hQinvBckgPCPsi%d",psi.first));
        if (hSignGK != nullptr && hBckgGK != nullptr && hSignPC != nullptr && hBckgPC != nullptr)
        {
            TH1D *hRatGK = new TH1D(*hSignGK);
            hRatGK->Divide(hBckgGK);
            hRatGK->SetName(TString::Format("hQinvRatGKPsi%d",psi.first));
            hRatGK->SetTitle(TString::Format("#phi #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})",psi.second.Data()));
            prepareGraph(hRatGK,JJColor::fWutMainColors[1]);

            TH1D *hRatPC = new TH1D(*hSignPC);
            hRatPC->Divide(hBckgPC);
            hRatPC->SetName(TString::Format("hQinvRatPCPsi%d",psi.first));
            hRatPC->SetTitle(TString::Format("#phi #in %s;q_{inv} [MeV/c];CF_{smeared}(q_{inv})",psi.second.Data()));
            prepareGraph(hRatPC,JJColor::fWutMainColors[1]);

            TH1D *hRat = new TH1D(*hRatGK);
            hRat->Divide(hRatPC);
            SetErrors(hRat,hRatGK,hRatPC);
            hRat->SetName(TString::Format("hQinvMomResPsi%d",psi.first));
            hRat->SetTitle(TString::Format("#phi #in %s;q_{inv} [MeV/c];CF_{raw}(q_{inv})/CF_{smeared}(q_{inv})",psi.second.Data()));
            prepareGraph(hRat,JJColor::fWutMainColors[1]);

            hRat->Write();
            hRatGK->Write();
            hRatPC->Write();
        }
    }
    canvPsi->Write();

    TH1D *hSignGK = inpFile->Get<TH1D>("hQinvSignGKInteg");
    TH1D *hBckgGK = inpFile->Get<TH1D>("hQinvBckgGKInteg");
    TH1D *hSignPC = inpFile->Get<TH1D>("hQinvSignPCInteg");
    TH1D *hBckgPC = inpFile->Get<TH1D>("hQinvBckgPCInteg");
    if (hSignGK != nullptr && hBckgGK != nullptr && hSignPC != nullptr && hBckgPC != nullptr)
    {
        TH1D *hRatGK = new TH1D(*hSignGK);
        hRatGK->Divide(hBckgGK);
        hRatGK->SetName("hQinvRatGKInteg");
        hRatGK->SetTitle(";q_{inv} [MeV/c];CF_{raw}(q_{inv})");
        prepareGraph(hRatGK,JJColor::fWutMainColors[1]);

        TH1D *hRatPC = new TH1D(*hSignPC);
        hRatPC->Divide(hBckgPC);
        hRatPC->SetName("hQinvRatPCInteg");
        hRatPC->SetTitle(";q_{inv} [MeV/c];CF_{smeared}(q_{inv})");
        prepareGraph(hRatPC,JJColor::fWutMainColors[1]);

        TH1D *hRat = new TH1D(*hRatGK);
        hRat->Divide(hRatPC);
        hRat->SetName("hQinvMomResInteg");
        hRat->SetTitle(";q_{inv} [MeV/c];CF_{raw}(q_{inv})/CF_{smeared}(q_{inv})");
        prepareGraph(hRat,JJColor::fWutMainColors[1]);

        hRat->Write();
        hRatGK->Write();
        hRatPC->Write();
    }
}