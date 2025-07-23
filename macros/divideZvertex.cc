#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"

void divideZvertex()
{
    TFile *inpFile = TFile::Open("../qaOutFile.root");
    TH1D *hZvertex = inpFile->Get<TH1D>("hZVertex");

    TF1 *fGaus1 = new TF1("fGaus1","gaus",-65.,-61.);
    TF1 *fGaus2 = new TF1("fGaus2","gaus",-60.3,-57.2);
    TF1 *fGaus3 = new TF1("fGaus3","gaus",-56,-53);
    TF1 *fGaus4 = new TF1("fGaus4","gaus",-53.,-49.9);
    TF1 *fGaus5 = new TF1("fGaus5","gaus",-48.5,-45.5);
    TF1 *fGaus6 = new TF1("fGaus6","gaus",-45.1,-43.);
    TF1 *fGaus7 = new TF1("fGaus7","gaus",-41.5,-38.5);
    TF1 *fGaus8 = new TF1("fGaus8","gaus",-38.,-36.);
    TF1 *fGaus9 = new TF1("fGaus9","gaus",-35.,-32.);
    TF1 *fGaus10 = new TF1("fGaus10","gaus",-30.2,-27.3);
    TF1 *fGaus11 = new TF1("fGaus11","gaus",-27.3,-24.);
    TF1 *fGaus12 = new TF1("fGaus12","gaus",-23.,-20.9);
    TF1 *fGaus13 = new TF1("fGaus13","gaus",-19.3,-17.);
    TF1 *fGaus14 = new TF1("fGaus14","gaus",-16.5,-13.5);
    TF1 *fGaus15 = new TF1("fGaus15","gaus",-12.7,-10.);

    hZvertex->Fit(fGaus1,"R");
    hZvertex->Fit(fGaus2,"R+");
    hZvertex->Fit(fGaus3,"R+");
    hZvertex->Fit(fGaus4,"R+");
    hZvertex->Fit(fGaus5,"R+");
    hZvertex->Fit(fGaus6,"R+");
    hZvertex->Fit(fGaus7,"R+");
    hZvertex->Fit(fGaus8,"R+");
    hZvertex->Fit(fGaus9,"R+");
    hZvertex->Fit(fGaus10,"R+");
    hZvertex->Fit(fGaus11,"R+");
    hZvertex->Fit(fGaus12,"R+");
    hZvertex->Fit(fGaus13,"R+");
    hZvertex->Fit(fGaus14,"R+");
    hZvertex->Fit(fGaus15,"R+");

    std::ofstream outFile("zVertexPeakAndSigmaFeb24.txt");
    outFile << fGaus1->GetParameter(1) << "\t" << fGaus1->GetParameter(2) << "\n";
    outFile << fGaus2->GetParameter(1) << "\t" << fGaus2->GetParameter(2) << "\n";
    outFile << fGaus3->GetParameter(1) << "\t" << fGaus3->GetParameter(2) << "\n";
    outFile << fGaus4->GetParameter(1) << "\t" << fGaus4->GetParameter(2) << "\n";
    outFile << fGaus5->GetParameter(1) << "\t" << fGaus5->GetParameter(2) << "\n";
    outFile << fGaus6->GetParameter(1) << "\t" << fGaus6->GetParameter(2) << "\n";
    outFile << fGaus7->GetParameter(1) << "\t" << fGaus7->GetParameter(2) << "\n";
    outFile << fGaus8->GetParameter(1) << "\t" << fGaus8->GetParameter(2) << "\n";
    outFile << fGaus9->GetParameter(1) << "\t" << fGaus9->GetParameter(2) << "\n";
    outFile << fGaus10->GetParameter(1) << "\t" << fGaus10->GetParameter(2) << "\n";
    outFile << fGaus11->GetParameter(1) << "\t" << fGaus11->GetParameter(2) << "\n";
    outFile << fGaus12->GetParameter(1) << "\t" << fGaus12->GetParameter(2) << "\n";
    outFile << fGaus13->GetParameter(1) << "\t" << fGaus13->GetParameter(2) << "\n";
    outFile << fGaus14->GetParameter(1) << "\t" << fGaus14->GetParameter(2) << "\n";
    outFile << fGaus15->GetParameter(1) << "\t" << fGaus15->GetParameter(2);
    outFile.close();
}