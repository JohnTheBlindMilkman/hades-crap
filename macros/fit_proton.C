#include <RtypesCore.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>

/*
 * fit_proton.C
 *
 *  Created on: 3 sie 2022
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef __CLING__
//#include "includes_hal.h"
//#include "/home/jedkol/Downloads/HAL/analysis/femto/corrfit/CorrFit.h"
// namespace Hal {};
#endif

/* #include "CorrFitKisiel.h"
#include "Femto1DCF.h"
#include "ChiSqMap2D.h"
#include "CorrFitGUI.h" */

using namespace Hal;

TH1D* MakeReducedHisto(TH1D *hInput, double maxRange)
{
  int lastBin = hInput->FindBin(maxRange);
  TH1D *hOutput = new TH1D(hInput->GetName(),hInput->GetTitle(),lastBin,0.,maxRange);
  for (int i = 1; i <= lastBin; i++) 
  {
    hOutput->SetBinContent(i,hInput->GetBinContent(i));
    hOutput->SetBinError(i,hInput->GetBinError(i));
  }

  return hOutput;
}

class MyFit : public Hal::CorrFitKisiel {
protected:
  Double_t CalculateCF(const Double_t* x, const Double_t* params) const {
    Double_t r  = params[Radius()];
    Double_t l  = params[Lambda()];
    Double_t n  = params[Norm()];
    Double_t cf = fMaps[0]->Eval(x[0], r);
    return n * ((cf - 1.0) * l + 1);
  }
  virtual void Paint(Bool_t repaint = kTRUE, Bool_t refresh = kTRUE) 
  {
      CorrFitKisiel::Paint(repaint, refresh);
      /* float val;
      for(int i = 1; i <= 40; i++)
      {
          val = (0.4/40)*i;
          cout << Eval(val) << endl;
      } */
  }

public:
  MyFit() : Hal::CorrFitKisiel(3) { SkipNumErrors(); };
  virtual ~MyFit() {};
};


void fit_proton() {
  TFile* file               = new TFile("mapa_pp.root");
    Hal::DividedHisto2D* mapH = (Hal::DividedHisto2D*) file->Get("new");
    CorrFitMapKstarRstar* map = new CorrFitMapKstarRstar(*((TH2D*) mapH->GetHist()), Hal::Femto::EKinematics::kLCMS);

    TFile *dataFile = TFile::Open("../output/1Dcorr_0_10_centforHAL.root");
    TH1D *hInp = MakeReducedHisto(dataFile->Get<TH1D>("hQinvRat4"),0.4);

  Femto1DCF* cf = new Femto1DCF("func", hInp->GetNbinsX(), hInp->GetXaxis()->GetXmin() , hInp->GetXaxis()->GetXmax(), Hal::Femto::EKinematics::kLCMS);
  // fill dummy data
  for (int i = 1; i <= hInp->GetNbinsX(); i++) {
    cf->GetNum()->SetBinContent(i, hInp->GetBinContent(i));
        cf->GetNum()->SetBinError(i, hInp->GetBinError(i));
        cf->GetDen()->SetBinContent(i, 1);
        cf->GetDen()->SetBinError(i, 0.);
  }

  MyFit* fit             = new MyFit();
  //CorrFit1DCF_Gauss* fit = new CorrFit1DCF_Gauss();
  fit->AddMap(map);
  fit->SetRLimits(1, 10);
  // fit->SetNormLimits(0.9, 1.1);
  fit->FixParameter(fit->Norm(), 1);
  fit->SetLambdaLimits(0.1, 1);
  fit->SetRange(0.01, 0.75);
  // fit->TraceFitting();
  fit->SetMinimizedFunc(Hal::CorrFit::EMinFunc::kChi2);
  fit->SetMinimizer(Hal::CorrFit::EMinAlgo::kMinuitCombined);
  TCanvas* c = new TCanvas();
  c->Divide(2, 1);
  c->cd(1);
  cf->Fit(fit);  // first fit
  ChiSqMap2D* chimap = fit->GetChiSquareMap(fit->Radius(), 20, 0.5, 10, fit->Lambda(), 20, 0, 1);
  Double_t Rpre      = chimap->GetEstX();
  Double_t Lampre    = chimap->GetEstY();
  c->cd(2);
  chimap->Draw("min");  // min -> draw minimum
  c->cd(1);


  if (Rpre < 1) Rpre = 3;           // for shitty data only!
  if (Lampre < 0.01) Lampre = 0.2;  // for shitty data only!
  fit->SetRLimits(Rpre - 1, Rpre + 1);
  fit->SetLambdaLimits(Lampre - 0.1, Lampre + 0.1);
  cf->Fit(fit);  // refit
  CorrFitDrawOptions opt(
    {
      CorrFitDrawOptions::eCommon::kDrawCF,
      CorrFitDrawOptions::eCommon::kNorm,
      CorrFitDrawOptions::eCommon::kLegend,
      CorrFitDrawOptions::eCommon::kDrawChi2  //,
                                              // CorrFitDrawOptions::eCommon::kNumError
    },
    {0, 2});
  // opt.SetLegendPos(0.5, 0.95, 0.5, 0.9);
  fit->SetDrawOption(opt);
  fit->Draw();

  Hal::CorrFitGUI* gui = new Hal::CorrFitGUI(fit);
}
