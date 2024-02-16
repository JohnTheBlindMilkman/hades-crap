#include <fstream>
#include <chrono>
#include <ctime>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TDatime.h"

double TanhFunc(double *x, double *par)
{
    return par[0] * tanh(par[1] * pow(x[0],par[2])) + par[3];
}

//derivative of TanhFunc over a
double DTanhFuncDa(double x, double b, double c)
{
    return tanh(b * pow(x,c));
}

//derivative of TanhFunc over b
double DTanhFuncDb(double x, double a, double b, double c)
{
    return a * pow(x,c) / (cosh(b * pow(x,c)) * cosh(b * pow(x,c))); //sech(x) = 1/cosh(x); there is no sech function in C++
}

//derivative of TanhFunc over c
double DTanhFuncDc(double x, double a, double b, double c)
{
    return a * b * pow(x,c) * log(x) / (cosh(b * pow(x,c)) * cosh(b * pow(x,c))); //sech(x) = 1/cosh(x); there is no sech function in C++
}

// propagation of uncertainty of TanhFunc
double ErrorPropagation(double x, double a, double sa, double b, double sb, double c, double sc, double d, double sd)
{
    return sqrt(DTanhFuncDa(x,b,c)*DTanhFuncDa(x,b,c)*sa*sa + DTanhFuncDb(x,a,b,c)*DTanhFuncDb(x,a,b,c)*sb*sb + DTanhFuncDc(x,a,b,c)*DTanhFuncDc(x,a,b,c)*sc*sc + sd*sd);
}

void SetErrors(TH1D *hout, TF1 *func)
{
    const int bins = hout->GetNbinsX();
    const double *param = func->GetParameters();
    const double *errs = func->GetParErrors();

    for (int i = 1; i <= bins; ++i)
    {
        hout->SetBinContent(i,func->Eval(hout->GetBinCenter(i)));
        hout->SetBinError(i,ErrorPropagation(hout->GetBinCenter(i),param[0],errs[0],param[1],errs[1],param[2],errs[2],param[3],errs[3]));
    }
}

void fitHGeant()
{
    constexpr std::array<int,5> ktArr{1,2,3,4,5};
    constexpr std::array<int,3> yArr{1,2,3};
    constexpr std::array<int,8> psiArr{1,2,3,4,5,6,7,8};

    std::vector<TH1D*> hSimKt(ktArr.size(),nullptr), hSimY(yArr.size(),nullptr), hSimPsi(psiArr.size(),nullptr);
    std::vector<TH1D*> hErrKt(ktArr.size(),nullptr), hErrY(yArr.size(),nullptr), hErrPsi(psiArr.size(),nullptr);
    std::vector<TF1*> fFitKt(ktArr.size(),nullptr), fFitY(yArr.size(),nullptr), fFitPsi(psiArr.size(),nullptr);

    TF1 *fitFunc = new TF1("fitFunc",TanhFunc,0,500,4);
    fitFunc->SetLineColor(EColor::kRed);
    fitFunc->SetParameters(0.876870,-14.6676,-1.11005,1.01467);
    fitFunc->SetParLimits(0,0.7,1.2);
    fitFunc->SetParLimits(1,-25,-14);
    fitFunc->SetParLimits(2,-1.5,-0.7);
    fitFunc->SetParLimits(3,0.9,1.1);

    TFile *inpFile = TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_HGeant.root");

    for (const int &kt : ktArr)
    {
        TString histName = TString::Format("hQinvRatKt%d",kt);
        TString errName = TString::Format("hQinvErrKt%d",kt);
        TString fitName = TString::Format("fQinvFitKt%d",kt);
        hSimKt[kt-1] = inpFile->Get<TH1D>(histName);
        hErrKt[kt-1] = new TH1D(errName,hSimKt[kt-1]->GetTitle(),10*hSimKt[kt-1]->GetNbinsX(),hSimKt[kt-1]->GetXaxis()->GetXmin(),hSimKt[kt-1]->GetXaxis()->GetXmax());

        hSimKt[kt-1]->GetYaxis()->SetRangeUser(0,1.2);
        hSimKt[kt-1]->GetXaxis()->SetRangeUser(0,500);

        hSimKt[kt-1]->Fit(fitFunc,"EMR0");
        fFitKt[kt-1] = new TF1(*fitFunc);
        fFitKt[kt-1]->SetName(fitName);
        SetErrors(hErrKt[kt-1],fFitKt[kt-1]);
    }
    for (const int &y : yArr)
    {
        TString histName = TString::Format("hQinvRatY%d",y);
        TString errName = TString::Format("hQinvErrY%d",y);
        TString fitName = TString::Format("fQinvFitY%d",y);
        hSimY[y-1] = inpFile->Get<TH1D>(histName);
        hErrY[y-1] = new TH1D(errName,hSimY[y-1]->GetTitle(),10*hSimY[y-1]->GetNbinsX(),hSimY[y-1]->GetXaxis()->GetXmin(),hSimY[y-1]->GetXaxis()->GetXmax());

        hSimY[y-1]->GetYaxis()->SetRangeUser(0,1.2);
        hSimY[y-1]->GetXaxis()->SetRangeUser(0,500);

        hSimY[y-1]->Fit(fitFunc,"EMR0");
        fFitY[y-1] = new TF1(*fitFunc);
        fFitY[y-1]->SetName(fitName);
        SetErrors(hErrY[y-1],fFitY[y-1]);
    }
    for (const int &psi : psiArr)
    {
        TString histName = TString::Format("hQinvRatPsi%d",psi);
        TString errName = TString::Format("hQinvErrPsi%d",psi);
        TString fitName = TString::Format("fQinvFitPsi%d",psi);
        hSimPsi[psi-1] = inpFile->Get<TH1D>(histName);
        hErrPsi[psi-1] = new TH1D(errName,hSimPsi[psi-1]->GetTitle(),10*hSimPsi[psi-1]->GetNbinsX(),hSimPsi[psi-1]->GetXaxis()->GetXmin(),hSimPsi[psi-1]->GetXaxis()->GetXmax());

        hSimPsi[psi-1]->GetYaxis()->SetRangeUser(0,1.2);
        hSimPsi[psi-1]->GetXaxis()->SetRangeUser(0,500);

        hSimPsi[psi-1]->Fit(fitFunc,"EMR0");
        fFitPsi[psi-1] = new TF1(*fitFunc);
        fFitPsi[psi-1]->SetName(fitName);
        SetErrors(hErrPsi[psi-1],fFitPsi[psi-1]);
    }

    TFile *outFile = TFile::Open("/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_HGeant_fit.root","RECREATE");
    for (const int &kt : ktArr)
    {
        hSimKt[kt-1]->Write();
        hErrKt[kt-1]->Write();
        fFitKt[kt-1]->Write();
    }
    for (const int &y : yArr)
    {
        hSimY[y-1]->Write();
        hErrY[y-1]->Write();
        fFitY[y-1]->Write();
    }
    for (const int &psi : psiArr)
    {
        hSimPsi[psi-1]->Write();
        hErrPsi[psi-1]->Write();
        fFitPsi[psi-1]->Write();
    }
}