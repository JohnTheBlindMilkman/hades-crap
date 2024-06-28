#include "TFile.h"
#include "TH3D.h"
#include "TF3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "../Externals/Palettes.hxx"

enum class Projection{Out,Side,Long};

double TanhFunc(double *x, double *par)
{
    return par[0] * tanh(par[1] * (x[0] - par[4]) + par[2] * (x[1] - par[5]) + par[3] * (x[2] - par[6])) + (1 - par[0]);
}

double LogisticFunc(double *x, double *par)
{
    return 1/(1 + par[0] * exp(- (x[0] - par[6]) - (par[2] * x[1] - par[6]) - (x[2] - par[6])) + par[1] * exp(- (par[3] * x[0] - par[7]) - (par[4] * x[1] - par[7]) - (par[5] * x[2] - par[7])));
}

//derivative of TanhFunc over a
double DTanhFuncDa(double x, double y, double z, double b, double c, double d)
{
    return tanh(b*x*x + c*y*y + d*z*z);
}

//derivative of TanhFunc over b
double DTanhFuncDb(double x, double y, double z, double a, double b, double c, double d)
{
    return a*x*x/(cosh(b*x*x + c*y*y + d*z*z)*cosh(b*x*x + c*y*y + d*z*z)); //sech(x) = 1/cosh(x); there is no sech function in C++
}

//derivative of TanhFunc over c
double DTanhFuncDc(double x, double y, double z, double a, double b, double c, double d)
{
    return a*y*y/(cosh(b*x*x + c*y*y + d*z*z)*cosh(b*x*x + c*y*y + d*z*z)); //sech(x) = 1/cosh(x); there is no sech function in C++
}

//derivative of TanhFunc over d
double DTanhFuncDd(double x, double y, double z, double a, double b, double c, double d)
{
    return a*z*z/(cosh(b*x*x + c*y*y + d*z*z)*cosh(b*x*x + c*y*y + d*z*z)); //sech(x) = 1/cosh(x); there is no sech function in C++
}

// propagation of uncertainty of TanhFunc
double ErrorPropagation(double x, double y, double z, double a, double sa, double b, double sb, double c, double sc, double d, double sd, double se)
{
    return sqrt(DTanhFuncDa(x,y,z,b,c,d)*DTanhFuncDa(x,y,z,b,c,d)*sa*sa + 
    DTanhFuncDb(x,y,z,a,b,c,d)*DTanhFuncDb(x,y,z,a,b,c,d)*sb*sb + 
    DTanhFuncDc(x,y,z,a,b,c,d)*DTanhFuncDc(x,y,z,a,b,c,d)*sc*sc + 
    DTanhFuncDd(x,y,z,a,b,c,d)*DTanhFuncDd(x,y,z,a,b,c,d)*sd*sd +
    se*se);
}

void SetErrors(TH3D *hout, TF3 *func)
{
    const int binsX = hout->GetNbinsX();
    const int binsY = hout->GetNbinsY();
    const int binsZ = hout->GetNbinsZ();
    const double *param = func->GetParameters();
    const double *errs = func->GetParErrors();

    for (int i = 1; i <= binsX; ++i)
        for (int j = 1; j <= binsY; ++j)
            for (int k = 1; k <= binsZ; ++k)
            {
                hout->SetBinContent(i,j,k,func->Eval(hout->GetXaxis()->GetBinCenter(i),hout->GetYaxis()->GetBinCenter(j),hout->GetZaxis()->GetBinCenter(k)));
                hout->SetBinError(i,j,k,ErrorPropagation(hout->GetXaxis()->GetBinCenter(i),hout->GetYaxis()->GetBinCenter(j),hout->GetZaxis()->GetBinCenter(k),param[0],errs[0],param[1],errs[1],param[2],errs[2],param[3],errs[3],errs[4]));
            }
}

TH1D* GetFitProjection(TF3 *fit, Projection proj)
{
    TH1D *hist;
    const int nbins = 1000;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    fit->GetRange(Xmin,Ymin,Zmin,Xmax,Ymax,Zmax);

    switch (proj)
    {
        case Projection::Out :
            hist = new TH1D("hOut","",nbins,Xmin,Xmax);
            for (int i = 1; i <= nbins; ++i)
                hist->SetBinContent(i,fit->Eval(hist->GetBinCenter(i),0,0));
            break;

        case Projection::Side :
            hist = new TH1D("hSide","",nbins,Ymin,Ymax);
            for (int i = 1; i <= nbins; ++i)
                hist->SetBinContent(i,fit->Eval(0,hist->GetBinCenter(i),0));
            break;

        case Projection::Long :
            hist = new TH1D("hLong","",nbins,Zmin,Zmax);
            for (int i = 1; i <= nbins; ++i)
                hist->SetBinContent(i,fit->Eval(0,0,hist->GetBinCenter(i)));
            break;
    }

    return hist;
}

void fitHGeant3D()
{
    const TString fileName = "/u/kjedrzej/hades-crap/output/3Dcorr_0_10_cent_HGeant_Integ.root";
    const TString outputFile = "/u/kjedrzej/hades-crap/output/3Dcorr_0_10_cent_HGeant_Integ_fit.root";

    TF3 *fitFunc = new TF3("fitFunc",TanhFunc,0,500,0,500,0,500,7);
    fitFunc->SetLineColor(JJColor::fWutSecondaryColors[3]);
    fitFunc->SetParameters(0.4,0.2,0.025,0.2,-3,-120,-3);
    fitFunc->SetParLimits(0,0.1,1.6);
    fitFunc->SetParLimits(1,0.05,0.3);
    fitFunc->SetParLimits(2,0.01,0.05);
    fitFunc->SetParLimits(3,0.05,0.3);
    fitFunc->SetParLimits(4,-50,0);
    fitFunc->SetParLimits(5,-200,-50);
    fitFunc->SetParLimits(6,-10,0);

    TFile *inpFile = TFile::Open(fileName);

    TH3D *hSim = inpFile->Get<TH3D>("hQoslRatInteg");
    hSim->Fit(fitFunc,"EMR0");
    TH3D *hErr = new TH3D(*hSim);
    hErr->SetName("hQoslErrInteg");
    SetErrors(hErr,fitFunc);

    TH1D *hOutFit = GetFitProjection(fitFunc,Projection::Out);
    TH1D *hSideFit = GetFitProjection(fitFunc,Projection::Side);
    TH1D *hLongFit = GetFitProjection(fitFunc,Projection::Long);
    TH1D *hOutSim = inpFile->Get<TH1D>("hQoutRatInteg");
    TH1D *hSideSim = inpFile->Get<TH1D>("hQsideRatInteg");
    TH1D *hLongSim = inpFile->Get<TH1D>("hQlongRatInteg");

    TFile *otpFile = TFile::Open(outputFile,"RECREATE");
    hSim->Write();
    fitFunc->Write();
    hErr->Write();
    
    TCanvas *canv = new TCanvas("canv","",1800,600);
    canv->Divide(3,1);

    canv->cd(1);
    hOutSim->Draw();
    hOutFit->Draw("c same");

    canv->cd(2);
    hSideSim->Draw();
    hSideFit->Draw("c same");

    canv->cd(3);
    hLongSim->Draw();
    hLongFit->Draw("c same");

    canv->Write();
}