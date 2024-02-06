#include <fstream>

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

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

void CorrectHisto(TH1D *hout, TF1 *fMod)
{
    const int iterMax = hout->GetNbinsX();
    double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
    for (int i = 1; i <= iterMax; i++)
    {
        vErr = 0;
        vNum = hout->GetBinContent(i);
        eNum = hout->GetBinError(i);
        vDen = fMod->Eval(hout->GetBinCenter(i));
        eDen = ErrorPropagation(hout->GetBinCenter(i),fMod->GetParameter(0),fMod->GetParError(0),fMod->GetParameter(1),fMod->GetParError(1),fMod->GetParameter(2),fMod->GetParError(2),fMod->GetParameter(3),fMod->GetParError(3));

        // propagation of uncertainty for a function num/den with onclusion of the correlation between the constituents
        if (fabs(vDen) > 0.)
            vErr = sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));

        hout->SetBinContent(i,vNum/vDen);
        hout->SetBinError(i,vErr);
    }
}

void drawDRProton1D()
{
    gStyle->SetPalette(kPastel);

    const TString inpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent.root";
    const TString otpFile = "/home/jedkol/lxpool/hades-crap/output/1Dcorr_0_10_cent_DR.root";
    const std::string parFile = "/home/jedkol/lxpool/hades-crap/macros/DRparams.txt";
    constexpr std::array<int,5> ktArr{1,2,3,4,5};
    constexpr std::array<int,3> yArr{1,2,3};
    constexpr std::array<int,8> psiArr{1,2,3,4,5,6,7,8};
    const int rebin = 1;

    std::vector<TH1D*> hCFkt(ktArr.size(),nullptr), hCFy(yArr.size(),nullptr), hCFpsi(psiArr.size(),nullptr);
    TFile *fInpData,*fOtpFile;
    std::string name;
    std::array<double,4> parArr,errArr;
    TLine line(0,1,3000,1);
    line.SetLineColor(kGray);
    line.SetLineStyle(kDashed);

    fInpData = TFile::Open(inpFile);

    std::ifstream fParFile(parFile);
    if (!fParFile.is_open())
        throw std::runtime_error("Parameter file could not be opened");

    for (const int &i : {0,1,2})
    {
        fParFile.ignore(1000,'\n');
    }

    for(const int &kt : ktArr)
    {
        fParFile >> name;
        std::cout << name << "\n";
        fParFile >> parArr[0] >> errArr[0];
        fParFile >> parArr[1] >> errArr[1];
        fParFile >> parArr[2] >> errArr[2];
        fParFile >> parArr[3] >> errArr[3];

        TF1 *func = new TF1("func",TanhFunc,0,500,4);
        func->SetParameters(parArr.data());
        func->SetParErrors(errArr.data());

        hCFkt[kt-1] = fInpData->Get<TH1D>(name.c_str());
        CorrectHisto(hCFkt[kt-1],func);
    }
    for(const int &y : yArr)
    {
        fParFile >> name;
        std::cout << name << "\n";
        fParFile >> parArr[0] >> errArr[0];
        fParFile >> parArr[1] >> errArr[1];
        fParFile >> parArr[2] >> errArr[2];
        fParFile >> parArr[3] >> errArr[3];

        TF1 *func = new TF1("func",TanhFunc,0,500,4);
        func->SetParameters(parArr.data());
        func->SetParErrors(errArr.data());

        hCFy[y-1] = fInpData->Get<TH1D>(name.c_str());
        CorrectHisto(hCFy[y-1],func);
    }
    for(const int &psi : psiArr)
    {
        fParFile >> name;
        std::cout << name << "\n";
        fParFile >> parArr[0] >> errArr[0];
        fParFile >> parArr[1] >> errArr[1];
        fParFile >> parArr[2] >> errArr[2];
        fParFile >> parArr[3] >> errArr[3];

        TF1 *func = new TF1("func",TanhFunc,0,500,4);
        func->SetParameters(parArr.data());
        func->SetParErrors(errArr.data());

        hCFpsi[psi-1] = fInpData->Get<TH1D>(name.c_str());
        CorrectHisto(hCFpsi[psi-1],func);
    }

    fOtpFile = TFile::Open(otpFile,"RECREATE");

    TCanvas *canvKt = new TCanvas("canvKt","",1600,900);
    for (auto hist : hCFkt)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFkt.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvKt->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvKt->Write();

    TCanvas *canvY = new TCanvas("canvY","",1600,900);
    for (auto hist : hCFy)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFy.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvY->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvY->Write();

    TCanvas *canvPsi = new TCanvas("canvPsi","",1600,900);
    for (auto hist : hCFpsi)
    {
        if (hist != nullptr)
        {
            hist->Write();
            if (hist == hCFpsi.front())
            {
                hist->GetYaxis()->SetRangeUser(0,2);
                hist->Draw("hist pe pmc plc");
            }
            else
                hist->Draw("hist pe pmc plc same");
        }
    }
    canvPsi->BuildLegend(0.2,0.2,0.5,0.5,"","p");
    line.Draw("same");
    canvPsi->Write();
}