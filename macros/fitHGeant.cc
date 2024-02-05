#include <fstream>
#include <chrono>
#include <ctime>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TDatime.h"

double LogFunc(double *x,double *par)
{
    return par[0] * log(par[1]*x[0] + par[2]) + par[3];
}

double InvFunc(double *x, double *par)
{
    return par[0]/(par[1] * pow(x[0],par[2])) + par[3];
}

double TanhFunc(double *x, double *par) // this one works best
{
    return par[0] * tanh(par[1] * pow(x[0],par[2])) + par[3];
}

double ExpFunc(double *x, double *par)
{
    return par[0] * exp(par[1] * pow(x[0],par[2])) + par[3];
}

void fitHGeant()
{
    constexpr std::array<int,5> ktArr{1,2,3,4,5};
    constexpr std::array<int,3> yArr{1,2,3};
    constexpr std::array<int,8> psiArr{1,2,3,4,5,6,7,8};

    std::ofstream outFile("DRparams.txt");
    std::time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    outFile << "Fit result file, created on " << std::ctime(&currentTime) << "\n\n";

    TF1 *fitFunc = new TF1("fitFunc",TanhFunc,0,500,4);
    fitFunc->SetLineColor(EColor::kRed);
    fitFunc->SetParameters(0.876870,-14.6676,-1.11005,1.01467);
    fitFunc->SetParLimits(0,0.7,1.2);
    fitFunc->SetParLimits(1,-25,-14);
    fitFunc->SetParLimits(2,-1.5,-0.7);
    fitFunc->SetParLimits(3,0.9,1.1);

    TFile *inpFile = TFile::Open("/home/jedkol/Downloads/HADES/HADES-CrAP/output/1Dcorr_0_10_cent_HGeant.root");

    for (const int &kt : ktArr)
    {
        TString histName = TString::Format("hQinvRatKt%d",kt);
        TH1D *hSim = inpFile->Get<TH1D>(histName);
        hSim->GetYaxis()->SetRangeUser(0,1.2);
        hSim->GetXaxis()->SetRangeUser(0,500);

        hSim->Fit(fitFunc,"EMR");

        outFile <<  histName << "\n";
        outFile << fitFunc->GetParameter(0) << "\t" << fitFunc->GetParError(0) << "\n";
        outFile << fitFunc->GetParameter(1) << "\t" << fitFunc->GetParError(1) << "\n";
        outFile << fitFunc->GetParameter(2) << "\t" << fitFunc->GetParError(2) << "\n";
        outFile << fitFunc->GetParameter(3) << "\t" << fitFunc->GetParError(3) << "\n";
    }
    for (const int &y : yArr)
    {
        TString histName = TString::Format("hQinvRatY%d",y);
        TH1D *hSim = inpFile->Get<TH1D>(histName);
        hSim->GetYaxis()->SetRangeUser(0,1.2);
        hSim->GetXaxis()->SetRangeUser(0,500);

        hSim->Fit(fitFunc,"EMR");

        outFile <<  histName << "\n";
        outFile << fitFunc->GetParameter(0) << "\t" << fitFunc->GetParError(0) << "\n";
        outFile << fitFunc->GetParameter(1) << "\t" << fitFunc->GetParError(1) << "\n";
        outFile << fitFunc->GetParameter(2) << "\t" << fitFunc->GetParError(2) << "\n";
        outFile << fitFunc->GetParameter(3) << "\t" << fitFunc->GetParError(3) << "\n";
    }
    for (const int &psi : psiArr)
    {
        TString histName = TString::Format("hQinvRatPsi%d",psi);
        TH1D *hSim = inpFile->Get<TH1D>(histName);
        hSim->GetYaxis()->SetRangeUser(0,1.2);
        hSim->GetXaxis()->SetRangeUser(0,500);

        hSim->Fit(fitFunc,"EMR");

        outFile <<  histName << "\n";
        outFile << fitFunc->GetParameter(0) << "\t" << fitFunc->GetParError(0) << "\n";
        outFile << fitFunc->GetParameter(1) << "\t" << fitFunc->GetParError(1) << "\n";
        outFile << fitFunc->GetParameter(2) << "\t" << fitFunc->GetParError(2) << "\n";
        outFile << fitFunc->GetParameter(3) << "\t" << fitFunc->GetParError(3) << "\n";
    }
}