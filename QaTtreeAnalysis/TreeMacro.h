//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 15 15:45:01 2011 by ROOT version 5.28/00b
// from TTree T/QA Tree
// found on file: be1122718423901.hld_dst_aug11_Tree.root
//////////////////////////////////////////////////////////

#ifndef AnalysisMacro_h
#define AnalysisMacro_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "hqatreemaker.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

class AnalysisMacro {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   HQAFileInfoTree *TFileInfo;
   HQAStartTree    *TStart;
   HQARichTree     *TRich;
   HQAMdcTree      *TMdc;
   HQATofTree      *TTof;
   HQARpcTree      *TRpc;
   HQAShowerTree   *TShw;
   HQAPhysicsTree  *TPhy;

   // List of branches
   TBranch        *b_TFileInfo;   //!
   TBranch        *b_TStart;   //!
   TBranch        *b_TRich;   //!
   TBranch        *b_TMdc;   //!
   TBranch        *b_TTof;   //!
   TBranch        *b_TRpc;   //!
   TBranch        *b_TShw;   //!
   TBranch        *b_TPhy;   //!


   


   AnalysisMacro(TFile *file, TTree *tree=0);
   virtual ~AnalysisMacro();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     pimpHist(TH1F* hist,TString YTitle,Float_t y_low, Float_t y_up, vector<Int_t>& vec,map<Int_t,TString>& m);

   virtual void     finalize(TFile *file);

   TFile *fn;

   TDirectory *Physics;
   TDirectory *Start;
   TDirectory *RICH;
   TDirectory *MDC;
};

//#endif

//#ifdef AnalysisMacro_cxx
AnalysisMacro::AnalysisMacro(TFile *file, TTree *tree)  : fn(file)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("be1122718423901.hld_dst_aug11_Tree.root");
      if (!f) {
         f = new TFile("be1122718423901.hld_dst_aug11_Tree.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);

   Physics = fn->mkdir("Physics");
   Start   = fn->mkdir("Start");
   RICH    = fn->mkdir("RICH");
   MDC     = fn->mkdir("MDC");
}

AnalysisMacro::~AnalysisMacro()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalysisMacro::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalysisMacro::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnalysisMacro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   TFileInfo = 0;
   TStart = 0;
   TRich = 0;
   TMdc = 0;
   TTof = 0;
   TRpc = 0;
   TShw = 0;
   TPhy = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("TFileInfo", &TFileInfo, &b_TFileInfo);
   fChain->SetBranchAddress("TStart", &TStart, &b_TStart);
   fChain->SetBranchAddress("TRich", &TRich, &b_TRich);
   fChain->SetBranchAddress("TMdc", &TMdc, &b_TMdc);
   fChain->SetBranchAddress("TTof", &TTof, &b_TTof);
   fChain->SetBranchAddress("TRpc", &TRpc, &b_TRpc);
   fChain->SetBranchAddress("TShw", &TShw, &b_TShw);
   fChain->SetBranchAddress("TPhy", &TPhy, &b_TPhy);
   Notify();

}

Bool_t AnalysisMacro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnalysisMacro::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnalysisMacro::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

Int_t getIndexRunId(Int_t runid, vector<Int_t>& vec)
{

    Int_t index = -1;
    for(UInt_t i = 0; i < vec.size(); i ++){
	if(vec.at(i) == runid) {
	    index = i;
	    break;
	}
    }

    return index;
}
void AnalysisMacro::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L AnalysisMacro.C
//      Root > AnalysisMacro t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
    //by  b_branchname->GetEntry(ientry); //read only this branch


    //TFile *outputfile = new TFile("QAHistograms_1.root","RECREATE");


    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    cout <<"Entries in Tree: " <<nentries <<endl;

    vector<Int_t> vecRunId;  //Vector which stores the time of files
    map<Int_t,TString> mRunIdToName;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t i=0; i<nentries;i++) {
	Long64_t ientry = LoadTree(i);
	if (ientry < 0) break;
	nb = fChain->GetEntry(i);   nbytes += nb;
	// if (Cut(ientry) < 0) continue;

	vecRunId.push_back(TFileInfo->fTRunId + TFileInfo->fTEvB);
	mRunIdToName[TFileInfo->fTRunId + TFileInfo->fTEvB] = TFileInfo->fTFileName;

	//cout <<TFileInfo->fTFileName <<endl;

    }

    sort(vecRunId.begin(),vecRunId.end());

    Int_t nFile = vecRunId.size();

    Char_t buffer[300];
    const Int_t nStartStrips = 16;

    // List of histograms*************************************************************************************************************

    //Information from particle candidate#############################
    /* TH1F  *hProtonMassMean                 =new TH1F("hProtonMassMean","hProtonMassMean",nFile,0,nFile);
    TH1F  *hPipMassMean                    =new TH1F("hPipMassMean","hPipMassMean",nFile,0,nFile);
    TH1F  *hPimMassMean                    =new TH1F("hPimMassMean","hPimMassMean",nFile,0,nFile);
    TH1F  *hProtonMassMeanRpc              =new TH1F("hProtonMassMeanRpc","hProtonMassMeanRpc",nFile,0,nFile);
    TH1F  *hPipMassMeanRpc                 =new TH1F("hPipMassMeanRpc","hPipMassMeanRpc",nFile,0,nFile);
    TH1F  *hPimMassMeanRpc                 =new TH1F("hPimMassMeanRpc","hPimMassMeanRpc",nFile,0,nFile);
    TH1F  *hProtonMassMeanTof              =new TH1F("hProtonMassMeanTof","hProtonMassMeanTof",nFile,0,nFile);
    TH1F  *hPipMassMeanTof                 =new TH1F("hPipMassMeanTof","hPipMassMeanTof",nFile,0,nFile);
    TH1F  *hPimMassMeanTof                 =new TH1F("hPimMassMeanTof","hPimMassMeanTof",nFile,0,nFile);

    TH1F  *hTimeLepSumSys0Mean             =new TH1F("hTimeLepSumSys0Mean","hTimeLepSumSys0Mean",nFile,0,nFile);
    TH1F  *hTimeLepSumSys0Max              =new TH1F("hTimeLepSumSys0Max","hTimeLepSumSys0Max",nFile,0,nFile);
    TH1F  *hTimeLepSumSys0Sig              =new TH1F("hTimeLepSumSys0Sig","hTimeLepSumSys0Sig",nFile,0,nFile);
    TH1F  *hTimeLepSumSys1Mean             =new TH1F("hTimeLepSumSys1Mean","hTimeLepSumSys1Mean",nFile,0,nFile);
    TH1F  *hTimeLepSumSys1Max              =new TH1F("hTimeLepSumSys1Max","hTimeLepSumSys1Max",nFile,0,nFile);
    TH1F  *hTimeLepSumSys1Sig              =new TH1F("hTimeLepSumSys1Sig","hTimeLepSumSys1Sig",nFile,0,nFile);
    TH1F  *hTimeLepStaSys0Mean[nStartStrips];
    TH1F  *hTimeLepStaSys0Max[nStartStrips];
    TH1F  *hTimeLepStaSys0Sig[nStartStrips];
    TH1F  *hTimeLepStaSys1Mean[nStartStrips];
    TH1F  *hTimeLepStaSys1Max[nStartStrips];
    TH1F  *hTimeLepStaSys1Sig[nStartStrips];
    for(Int_t i=0; i<nStartStrips;i++)
    {
	sprintf(buffer,"hTimeLepStaSys0Mean[%i]",i);
	hTimeLepStaSys0Mean[i]             =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimeLepStaSys0Max[%i]",i);
	hTimeLepStaSys0Max[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimeLepStaSys0Sig[%i]",i);
	hTimeLepStaSys0Sig[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimeLepStaSys1Mean[%i]",i);
	hTimeLepStaSys1Mean[i]             =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimeLepStaSys1Max[%i]",i);
	hTimeLepStaSys1Max[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimeLepStaSys1Sig[%i]",i);
	hTimeLepStaSys1Sig[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
    }
    TH1F  *hTimePiSumSys0Mean             =new TH1F("hTimePiSumSys0Mean","hTimePiSumSys0Mean",nFile,0,nFile);
    TH1F  *hTimePiSumSys0Max              =new TH1F("hTimePiSumSys0Max","hTimePiSumSys0Max",nFile,0,nFile);
    TH1F  *hTimePiSumSys0Sig              =new TH1F("hTimePiSumSys0Sig","hTimePiSumSys0Sig",nFile,0,nFile);
    TH1F  *hTimePiSumSys1Mean             =new TH1F("hTimePiSumSys1Mean","hTimePiSumSys1Mean",nFile,0,nFile);
    TH1F  *hTimePiSumSys1Max              =new TH1F("hTimePiSumSys1Max","hTimePiSumSys1Max",nFile,0,nFile);
    TH1F  *hTimePiSumSys1Sig              =new TH1F("hTimePiSumSys1Sig","hTimePiSumSys1Sig",nFile,0,nFile);
    TH1F  *hTimePiStaSys0Mean[nStartStrips];
    TH1F  *hTimePiStaSys0Max[nStartStrips];
    TH1F  *hTimePiStaSys0Sig[nStartStrips];
    TH1F  *hTimePiStaSys1Mean[nStartStrips];
    TH1F  *hTimePiStaSys1Max[nStartStrips];
    TH1F  *hTimePiStaSys1Sig[nStartStrips];
    for(Int_t i=0; i<nStartStrips;i++)
    {
	sprintf(buffer,"hTimePiStaSys0Mean[%i]",i);
	hTimePiStaSys0Mean[i]             =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimePiStaSys0Max[%i]",i);
	hTimePiStaSys0Max[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimePiStaSys0Sig[%i]",i);
	hTimePiStaSys0Sig[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimePiStaSys1Mean[%i]",i);
	hTimePiStaSys1Mean[i]             =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimePiStaSys1Max[%i]",i);
	hTimePiStaSys1Max[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hTimePiStaSys1Sig[%i]",i);
	hTimePiStaSys1Sig[i]              =new TH1F(buffer,buffer,nFile,0,nFile);
    }


    // Information from Start detector################################
    TH1F  *hHitTimeMean                    =new TH1F("hHitTimeMean","hHitTimeMean",nFile,0,nFile);
    TH1F  *hHitTimeRms                     =new TH1F("hHitTimeRms","hHitTimeRms",nFile,0,nFile);
    TH1F  *hHitToNoveto                    =new TH1F("hHitToNoveto","hHitToNoveto",nFile,0,nFile);
    TH1F  *hVertexX                        =new TH1F("hVertexX","hVertexX",nFile,0,nFile);
    TH1F  *hVertexY                        =new TH1F("hVertexY","hVertexY",nFile,0,nFile);
    TH1F  *hVertexZ                        =new TH1F("hVertexZ","hVertexZ",nFile,0,nFile);
    TH1F  *hStCalTime[nStartStrips];
    TH1F  *hVtCalTime[nStartStrips];
    for(Int_t i=0; i<nStartStrips; i++)
    {
	sprintf(buffer,"hStCalTime[%i]",i);
	hStCalTime[i]                      =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hVtCalTime[%i]",i);
	hVtCalTime[i]                      =new TH1F(buffer,buffer,nFile,0,nFile);
    }

    // Information from RICH detector
    TH1F  *hAvChargeMax[6];
    TH1F  *hNumPadsMax[6];
    TH1F  *hNumRingsMax[6];
    TH1F  *hAvChargeMean[6];
    TH1F  *hNumPadsMean[6];
    TH1F  *hNumRingsMean[6];
    TH1F  *hHVMeanRich[6];
    TH1F  *hHVSigmaRich[6];
    TH1F  *hHVMinRich[6];
    TH1F  *hHVMaxRich[6];
    TH1F  *hHVNDataRich[6];
    for(Int_t i=0; i<6; i++)
    {
	sprintf(buffer,"hAvChargeMax[%i]",i);
	hAvChargeMax[i]                     =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hNumPadsMax[%i]",i);
	hNumPadsMax[i]                      =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hNumRingsMax[%i]",i);
	hNumRingsMax[i]                     =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hAvChargeMean[%i]",i);
	hAvChargeMean[i]                    =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hNumPadsMean[%i]",i);
	hNumPadsMean[i]                     =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hNumRingsMean[%i]",i);
	hNumRingsMean[i]                    =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hHVMeanRich[%i]",i);
	hHVMeanRich[i]                      =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hHVSigmaRich[%i]",i);
	hHVSigmaRich[i]                     =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hHVMinRich[%i]",i);
	hHVMinRich[i]                       =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hHVMaxRich[%i]",i);
	hHVMaxRich[i]                       =new TH1F(buffer,buffer,nFile,0,nFile);
	sprintf(buffer,"hHVNDataRich[%i]",i);
	hHVNDataRich[i]                     =new TH1F(buffer,buffer,nFile,0,nFile);
    }
 */
    //Information from Mdc detector
    TH1F  *hMboCounts[6][4][16];
    TH1F  *hHVMeanMdc[6][4][12];
    TH1F  *hHVSigmaMdc[6][4][12];
    TH1F  *hHVMinMdc[6][4][12];
    TH1F  *hHVMaxMdc[6][4][12];
    TH1F  *hHVNDataMdc[6][4][12];
    TH1F  *hLayerCounts[6][4][6];
    TH1F  *hToTCalMean[6][4][6];
    TH1F  *hToTCalMax[6][4][6];
    TH1F  *hTime1CalMean[6][4][6];
    TH1F  *hTime1CalMax[6][4][6];
    TH1F  *hCountsRaw[6][4];
    TH1F  *hCountsCal[6][4];
    TH1F  *hSegCounter[6][2];
    TH1F  *hSegCounterFit[6][2];
    TH1F  *hSegChi2Mean[6][2];
    TH1F  *hSegChi2Max[6][2];
    
    for(Int_t i=0; i<6; i++)
    {
	for(Int_t j=0; j<4; j++)
	{
	    for(Int_t k=0; k<16; k++)
	    {
		sprintf(buffer,"hMboCounts[%i][%i][%i]",i,j,k);
		hMboCounts[i][j][k]         =new TH1F(buffer,buffer,nFile,0,nFile);
	    }
	    for(Int_t k=0; k<12; k++)
	    {
		sprintf(buffer,"hHVMeanMdc[%i][%i][%i]",i,j,k);
		hHVMeanMdc[i][j][k]         =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hHVSigmaMdc[%i][%i][%i]",i,j,k);
		hHVSigmaMdc[i][j][k]        =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hHVMinMdc[%i][%i][%i]",i,j,k);
		hHVMinMdc[i][j][k]          =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hHVMaxMdc[%i][%i][%i]",i,j,k);
		hHVMaxMdc[i][j][k]          =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hHVNDataMdc[%i][%i][%i]",i,j,k);
		hHVNDataMdc[i][j][k]        =new TH1F(buffer,buffer,nFile,0,nFile);
	    }

	    for(Int_t k=0; k<6; k++)
	    {
		sprintf(buffer,"hLayerCounts[%i][%i][%i]",i,j,k);
		hLayerCounts[i][j][k]       =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hToTCalMean[%i][%i][%i]",i,j,k);
		hToTCalMean[i][j][k]        =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hToTCalMax[%i][%i][%i]",i,j,k);
		hToTCalMax[i][j][k]         =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hTime1CalMean[%i][%i][%i]",i,j,k);
		hTime1CalMean[i][j][k]      =new TH1F(buffer,buffer,nFile,0,nFile);
		sprintf(buffer,"hTime1CalMax[%i][%i][%i]",i,j,k);
		hTime1CalMax[i][j][k]       =new TH1F(buffer,buffer,nFile,0,nFile);
	    }
	    sprintf(buffer,"hCountsRaw[%i][%i]",i,j);
	    hCountsRaw[i][j]                =new TH1F(buffer,buffer,nFile,0,nFile);
	    sprintf(buffer,"hCountsCal[%i][%i]",i,j);
	    hCountsCal[i][j]                =new TH1F(buffer,buffer,nFile,0,nFile);
	}
	for(Int_t j=0; j<2; j++)
	{
	    sprintf(buffer,"hSegCounter[%i][%i]",i,j);
	    hSegCounter[i][j]               =new TH1F(buffer,buffer,nFile,0,nFile);
            sprintf(buffer,"hSegCounterFit[%i][%i]",i,j);
	    hSegCounterFit[i][j]            =new TH1F(buffer,buffer,nFile,0,nFile);
            sprintf(buffer,"hSegChi2Mean[%i][%i]",i,j);
	    hSegChi2Mean[i][j]              =new TH1F(buffer,buffer,nFile,0,nFile);
            sprintf(buffer,"hSegChi2Max[%i][%i]",i,j);
	    hSegChi2Max[i][j]               =new TH1F(buffer,buffer,nFile,0,nFile);
	}
    }


    // End List of histograms**********************************************************************************************************


    for (Long64_t i=0; i<nentries;i++)
    {
	Long64_t ientry = LoadTree(i);
	if (ientry < 0) break;
	nb = fChain->GetEntry(i);   nbytes += nb;
	// if (Cut(ientry) < 0) continue;

	Int_t index = getIndexRunId(TFileInfo->fTRunId + TFileInfo->fTEvB, vecRunId);

	//cout<<index <<" "<< TFileInfo->fTFileName <<" "<<TFileInfo->fTRunId <<" "<< TFileInfo->fTEvB<<endl;
	Int_t bin=index+1;


	//Start of Filling Histograms*********************************************************************
	//Information from Particle candidate
	/* hProtonMassMean->SetBinContent(bin,TPhy->fTMass_proton);
	hPipMassMean->SetBinContent(bin,TPhy->fTMass_pip);
	hPimMassMean->SetBinContent(bin,TPhy->fTMass_pim);
	hProtonMassMeanRpc->SetBinContent(bin,TPhy->fTMass_proton_RPC);
	hPipMassMeanRpc->SetBinContent(bin,TPhy->fTMass_pip_RPC);
	hPimMassMeanRpc->SetBinContent(bin,TPhy->fTMass_pim_RPC);
	hProtonMassMeanTof->SetBinContent(bin,TPhy->fTMass_proton_TOF);
	hPipMassMeanTof->SetBinContent(bin,TPhy->fTMass_pip_TOF);
	hPimMassMeanTof->SetBinContent(bin,TPhy->fTMass_pim_TOF);

	hTimeLepSumSys0Mean->SetBinContent(bin,TPhy->fTTimeLepSumSys0Mean);
	hTimeLepSumSys0Max ->SetBinContent(bin,TPhy->fTTimeLepSumSys0Max);
	hTimeLepSumSys0Sig ->SetBinContent(bin,TPhy->fTTimeLepSumSys0Sig);
	hTimeLepSumSys1Mean->SetBinContent(bin,TPhy->fTTimeLepSumSys1Mean);
	hTimeLepSumSys1Max ->SetBinContent(bin,TPhy->fTTimeLepSumSys1Max);
	hTimeLepSumSys1Sig ->SetBinContent(bin,TPhy->fTTimeLepSumSys1Sig);
	for(Int_t i=0; i<nStartStrips; i++)
	{
	    hTimeLepStaSys0Mean[i]->SetBinContent(bin,TPhy->fTTimeLepStaSys0Mean[i]);
	    hTimeLepStaSys0Max[i] ->SetBinContent(bin,TPhy->fTTimeLepStaSys0Max[i]);
	    hTimeLepStaSys0Sig[i] ->SetBinContent(bin,TPhy->fTTimeLepStaSys0Sig[i]);
	    hTimeLepStaSys1Mean[i]->SetBinContent(bin,TPhy->fTTimeLepStaSys1Mean[i]);
	    hTimeLepStaSys1Max[i] ->SetBinContent(bin,TPhy->fTTimeLepStaSys1Max[i]);
	    hTimeLepStaSys1Sig[i] ->SetBinContent(bin,TPhy->fTTimeLepStaSys1Sig[i]);
	}
	hTimePiSumSys0Mean->SetBinContent(bin,TPhy->fTTimePiSumSys0Mean);
	hTimePiSumSys0Max ->SetBinContent(bin,TPhy->fTTimePiSumSys0Max);
	hTimePiSumSys0Sig ->SetBinContent(bin,TPhy->fTTimePiSumSys0Sig);
	hTimePiSumSys1Mean->SetBinContent(bin,TPhy->fTTimePiSumSys1Mean);
	hTimePiSumSys1Max ->SetBinContent(bin,TPhy->fTTimePiSumSys1Max);
	hTimePiSumSys1Sig ->SetBinContent(bin,TPhy->fTTimePiSumSys1Sig);
	for(Int_t i=0; i<nStartStrips; i++)
	{
	    hTimePiStaSys0Mean[i]->SetBinContent(bin,TPhy->fTTimePiStaSys0Mean[i]);
	    hTimePiStaSys0Max[i] ->SetBinContent(bin,TPhy->fTTimePiStaSys0Max[i]);
	    hTimePiStaSys0Sig[i] ->SetBinContent(bin,TPhy->fTTimePiStaSys0Sig[i]);
	    hTimePiStaSys1Mean[i]->SetBinContent(bin,TPhy->fTTimePiStaSys1Mean[i]);
	    hTimePiStaSys1Max[i] ->SetBinContent(bin,TPhy->fTTimePiStaSys1Max[i]);
	    hTimePiStaSys1Sig[i] ->SetBinContent(bin,TPhy->fTTimePiStaSys1Sig[i]);
	}

	// Information from Start detector
	hHitTimeMean->SetBinContent(bin,TStart->fTHitTimeMean);
	hHitTimeRms->SetBinContent(bin,TStart->fTHitTimeRms);
	hHitToNoveto->SetBinContent(bin,TStart->fTHitToNoveto);
	hVertexX->SetBinContent(bin,TStart->fTVertexX);
	hVertexY->SetBinContent(bin,TStart->fTVertexY);
	hVertexZ->SetBinContent(bin,TStart->fTVertexZ);
	for(Int_t i=0; i<nStartStrips; i++)
	{
	    hStCalTime[i]->SetBinContent(bin,TStart->fTStCalTime[i]);
	    hVtCalTime[i]->SetBinContent(bin,TStart->fTVtCalTime[i]);
	}

	// Information from RHICH detector
	for(Int_t i=0; i<6; i++)
	{
	    hAvChargeMax[i]->SetBinContent(bin,TRich->fTAvChargeMax[i]);
	    hNumPadsMax[i]->SetBinContent(bin,TRich->fTNumPadsMax[i]);
	    hNumRingsMax[i]->SetBinContent(bin,TRich->fTNumRingsMax[i]);
	    hAvChargeMean[i]->SetBinContent(bin,TRich->fTAvChargeMean[i]);
	    hNumPadsMean[i]->SetBinContent(bin,TRich->fTNumPadsMean[i]);
	    hNumRingsMean[i]->SetBinContent(bin,TRich->fTNumRingsMean[i]);
	    hHVMeanRich[i]->SetBinContent(bin,TRich->fTHVMean[i]);
	    hHVSigmaRich[i]->SetBinContent(bin,TRich->fTHVSigma[i]);
	    hHVMinRich[i]->SetBinContent(bin,TRich->fTHVMin[i]);
	    hHVMaxRich[i]->SetBinContent(bin,TRich->fTHVMax[i]);
	    hHVNDataRich[i]->SetBinContent(bin,TRich->fTHVNData[i]);
	} */

	// Information from Mdc detectors
	for(Int_t i=0; i<6; i++)
	{
	    for(Int_t j=0; j<4; j++)
	    {
		for(Int_t k=0; k<16; k++)
		{
		    hMboCounts[i][j][k]->SetBinContent(bin,TMdc->fTMboCounts[i][j][k]);
		}
		for(Int_t k=0; k<12; k++)
		{
		    hHVMeanMdc[i][j][k]->SetBinContent(bin,TMdc->fTHVMean[i][j][k]);
		    hHVSigmaMdc[i][j][k]->SetBinContent(bin,TMdc->fTHVSigma[i][j][k]);
		    hHVMinMdc[i][j][k]->SetBinContent(bin,TMdc->fTHVMin[i][j][k]);
		    hHVMaxMdc[i][j][k]->SetBinContent(bin,TMdc->fTHVMax[i][j][k]);
		    hHVNDataMdc[i][j][k]->SetBinContent(bin,TMdc->fTHVNData[i][j][k]);
		}
		for(Int_t k=0; k<6; k++)
		{
		    hLayerCounts[i][j][k]->SetBinContent(bin,TMdc->fTLayerCounts[i][j][k]);
		    hToTCalMean[i][j][k]->SetBinContent(bin,TMdc->fTToTCalMean[i][j][k]);
		    hToTCalMax[i][j][k]->SetBinContent(bin,TMdc->fTToTCalMax[i][j][k]);
		    hTime1CalMean[i][j][k]->SetBinContent(bin,TMdc->fTTime1CalMean[i][j][k]);
		    hTime1CalMax[i][j][k]->SetBinContent(bin,TMdc->fTTime1CalMax[i][j][k]);
		}
		hCountsRaw[i][j]->SetBinContent(bin,TMdc->fTCountsRaw[i][j]);
		hCountsCal[i][j]->SetBinContent(bin,TMdc->fTCountsCal[i][j]);
	    }
	    for(Int_t j=0; j<2; j++)
	    {
		hSegCounter[i][j]->SetBinContent(bin,TMdc->fTSegCounter[i][j]);
		hSegCounterFit[i][j]->SetBinContent(bin,TMdc->fTSegCounterFit[i][j]);
		hSegChi2Mean[i][j]->SetBinContent(bin,TMdc->fTSegChi2Mean[i][j]);
		hSegChi2Max[i][j]->SetBinContent(bin,TMdc->fTSegChi2Max[i][j]);
	    }
	}


	//End of Filling Histograms*********************************************************************
    }


    //Style the Histograms##############################################################################

    //Physics
    /* pimpHist(hProtonMassMean,"Mass Proton",900,960,vecRunId,mRunIdToName);
    pimpHist(hPipMassMean,"Mass #pi^{+}",100,180,vecRunId,mRunIdToName);
    pimpHist(hPimMassMean,"Mass #pi^{-}",-180,-100,vecRunId,mRunIdToName);
    pimpHist(hProtonMassMeanRpc,"Mass Proton in Rpc",900,960,vecRunId,mRunIdToName);
    pimpHist(hPipMassMeanRpc,"Mass #pi^{+} in Rpc",100,180,vecRunId,mRunIdToName);
    pimpHist(hPimMassMeanRpc,"Mass #pi^{-} in Rpc",-180,-100,vecRunId,mRunIdToName);
    pimpHist(hProtonMassMeanTof,"Mass Proton in Tof",900,960,vecRunId,mRunIdToName);
    pimpHist(hPipMassMeanTof,"Mass #pi^{+} in Tof",100,180,vecRunId,mRunIdToName);
    pimpHist(hPimMassMeanTof,"Mass #pi^{-} in Tof",-180,-100,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys0Mean,"Lepton Tof, Mean, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys0Max,"Lepton Tof, Max, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys0Sig,"Lepton Tof, RMS, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys1Mean,"Lepton Tof, Mean, System1",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys1Max,"Lepton Tof, Max, System1",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimeLepSumSys1Sig,"Lepton Tof, RMS, System1",-10,20,vecRunId,mRunIdToName);
    for(Int_t i=0; i<nStartStrips; i++)
    {
	sprintf(buffer,"Lepton Tof, Mean, System0 StartStrip %i",i);
	pimpHist(hTimeLepStaSys0Mean[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Lepton Tof, Max, System0 StartStrip %i",i);
	pimpHist(hTimeLepStaSys0Max[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Lepton Tof, RMS, System0 StartStrip %i",i);
	pimpHist(hTimeLepStaSys0Sig[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Lepton Tof, Mean, System1 StartStrip %i",i);
	pimpHist(hTimeLepStaSys1Mean[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Lepton Tof, Max, System1 StartStrip %i",i);
	pimpHist(hTimeLepStaSys1Max[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Lepton Tof, RMS, System1 StartStrip %i",i);
	pimpHist(hTimeLepStaSys1Sig[i],buffer,-10,20,vecRunId,mRunIdToName);
    }
    pimpHist(hTimePiSumSys0Mean,"Pion Tof, Mean, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimePiSumSys0Max,"Pion Tof, Max, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimePiSumSys0Sig,"Pion Tof, RMS, System0",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimePiSumSys1Mean,"Pion Tof, Mean, System1",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimePiSumSys1Max,"Pion Tof, Max, System1",-10,20,vecRunId,mRunIdToName);
    pimpHist(hTimePiSumSys1Sig,"Pion Tof, RMS, System1",-10,20,vecRunId,mRunIdToName);
    for(Int_t i=0; i<nStartStrips; i++)
    {
	sprintf(buffer,"Pion Tof, Mean, System0 StartStrip %i",i);
	pimpHist(hTimePiStaSys0Mean[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Pion Tof, Max, System0 StartStrip %i",i);
	pimpHist(hTimePiStaSys0Max[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Pion Tof, RMS, System0 StartStrip %i",i);
	pimpHist(hTimePiStaSys0Sig[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Pion Tof, Mean, System1 StartStrip %i",i);
	pimpHist(hTimePiStaSys1Mean[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Pion Tof, Max, System1 StartStrip %i",i);
	pimpHist(hTimePiStaSys1Max[i],buffer,-10,20,vecRunId,mRunIdToName);
	sprintf(buffer,"Pion Tof, RMS, System1 StartStrip %i",i);
	pimpHist(hTimePiStaSys1Sig[i],buffer,-10,20,vecRunId,mRunIdToName);
    }

    //Start
    pimpHist(hHitTimeMean,"Start Time, Mean",-2,2,vecRunId,mRunIdToName);
    pimpHist(hHitTimeRms,"Start Time, RMS",0,2,vecRunId,mRunIdToName);
    pimpHist(hHitToNoveto,"Num Start Hits/Num Veto Hits",0,5,vecRunId,mRunIdToName);
    pimpHist(hVertexX,"Vertex X, Mean",-3,3,vecRunId,mRunIdToName);
    pimpHist(hVertexY,"Vertex Y, Mean",-10,10,vecRunId,mRunIdToName);
    pimpHist(hVertexZ,"Vertex Z, Mean",-20,-40,vecRunId,mRunIdToName);
    for(Int_t i=0; i<nStartStrips; i++)
    {
	sprintf(buffer,"Start Time for Strip %i",i);
	pimpHist(hStCalTime[i],buffer,-5,5,vecRunId,mRunIdToName);
	sprintf(buffer,"Veto Time for Strip %i",i);
	pimpHist(hVtCalTime[i],buffer,-5,5,vecRunId,mRunIdToName);
    }

    //RICH
    for(Int_t i=0; i<6; i++)
    {
        sprintf(buffer,"Average Charge, Max, Sector[%i]",i);
	pimpHist(hAvChargeMax[i],buffer,0,500,vecRunId,mRunIdToName);
	sprintf(buffer,"Num of fired pads, Max, Sector[%i]",i);
	pimpHist(hNumPadsMax[i],buffer,0,20,vecRunId,mRunIdToName);
        sprintf(buffer,"Average Num of rings per event, Max, Sector[%i]",i);
	pimpHist(hNumRingsMax[i],buffer,-2,2,vecRunId,mRunIdToName);
        sprintf(buffer,"Average Charge, Mean, Sector[%i]",i);
	pimpHist(hAvChargeMean[i],buffer,0,500,vecRunId,mRunIdToName);
	sprintf(buffer,"Num of fired pads, Mean, Sector[%i]",i);
	pimpHist(hNumPadsMean[i],buffer,0,30,vecRunId,mRunIdToName);
        sprintf(buffer,"Average Num of rings per event, Mean, Sector[%i]",i);
	pimpHist(hNumRingsMean[i],buffer,-2,10,vecRunId,mRunIdToName);
        sprintf(buffer,"High Voltage, Mean, Sector[%i]",i);
	pimpHist(hHVMeanRich[i],buffer,0,3000,vecRunId,mRunIdToName);
        sprintf(buffer,"High Voltage, Sigma, Sector[%i]",i);
	pimpHist(hHVSigmaRich[i],buffer,-10,10,vecRunId,mRunIdToName);
        sprintf(buffer,"High Voltage, Min, Sector[%i]",i);
	pimpHist(hHVMinRich[i],buffer,0,3000,vecRunId,mRunIdToName);
        sprintf(buffer,"High Voltage, Max, Sector[%i]",i);
	pimpHist(hHVMaxRich[i],buffer,0,3000,vecRunId,mRunIdToName);
        sprintf(buffer,"High Voltage, Data Points, Sector[%i]",i);
	pimpHist(hHVNDataRich[i],buffer,-2,10,vecRunId,mRunIdToName);
    } */

    //MDC
    for(Int_t i=0; i<6; i++)
    {
	for(Int_t j=0; j<4; j++)
	{
	    for(Int_t k=0; k<16; k++)
	    {
		sprintf(buffer,"av. counts per event in MBO[%i][%i][%i]",i,j,k);
		pimpHist(hMboCounts[i][j][k],buffer,-2,300,vecRunId,mRunIdToName);
	    }
	    for(Int_t k=0; k<12; k++)
	    {
		sprintf(buffer,"High Voltage, Mean, [%i][%i][%i]",i,j,k);
                pimpHist(hHVMeanMdc[i][j][k],buffer,-2,2500,vecRunId,mRunIdToName);
		sprintf(buffer,"High Voltage, Sigma, [%i][%i][%i]",i,j,k);
                pimpHist(hHVSigmaMdc[i][j][k],buffer,-2,100,vecRunId,mRunIdToName);
		sprintf(buffer,"High Voltage, Min, [%i][%i][%i]",i,j,k);
                pimpHist(hHVMinMdc[i][j][k],buffer,-2,2500,vecRunId,mRunIdToName);
		sprintf(buffer,"High Voltage, Max, [%i][%i][%i]",i,j,k);
                pimpHist(hHVMaxMdc[i][j][k],buffer,-2,2500,vecRunId,mRunIdToName);
		sprintf(buffer,"High Voltage, Data Points, [%i][%i][%i]",i,j,k);
                pimpHist(hHVNDataMdc[i][j][k],buffer,-2,500,vecRunId,mRunIdToName);
	    }

	    for(Int_t k=0; k<6; k++)
	    {
                sprintf(buffer,"av. counts per event in layer [%i][%i][%i]",i,j,k);
		pimpHist(hLayerCounts[i][j][k],buffer,-2,300,vecRunId,mRunIdToName);
                sprintf(buffer,"ToT, Mean, in layer [%i][%i][%i]",i,j,k);
		pimpHist(hToTCalMean[i][j][k],buffer,-2,200,vecRunId,mRunIdToName);
		sprintf(buffer,"ToT, Max, in layer [%i][%i][%i]",i,j,k);
		pimpHist(hToTCalMax[i][j][k],buffer,-2,200,vecRunId,mRunIdToName);
                sprintf(buffer,"Time1, Mean, in layer [%i][%i][%i]",i,j,k);
		pimpHist(hTime1CalMean[i][j][k],buffer,-2,200,vecRunId,mRunIdToName);
		sprintf(buffer,"Time1, Max, in layer [%i][%i][%i]",i,j,k);
		pimpHist(hTime1CalMax[i][j][k],buffer,-2,200,vecRunId,mRunIdToName);
	    }
	    sprintf(buffer,"av. raw counts per event in sector [%i][%i]",i,j);
	    pimpHist(hCountsRaw[i][j],buffer,-2,500,vecRunId,mRunIdToName);
            sprintf(buffer,"av. cal counts per event in sector [%i][%i]",i,j);
	    pimpHist(hCountsCal[i][j],buffer,-2,500,vecRunId,mRunIdToName);
	}
	for(Int_t j=0; j<2; j++)
	{
            sprintf(buffer,"av. counts per event in segment [%i][%i]",i,j);
	    pimpHist(hSegCounter[i][j],buffer,-2,100,vecRunId,mRunIdToName);
            sprintf(buffer,"av. counts per fitted events in segment [%i][%i]",i,j);
	    pimpHist(hSegCounterFit[i][j],buffer,-2,100,vecRunId,mRunIdToName);
            sprintf(buffer,"chi2, Mean, segment [%i][%i]",i,j);
	    pimpHist(hSegChi2Mean[i][j],buffer,-2,20,vecRunId,mRunIdToName);
	    sprintf(buffer,"chi2, Max, segment[%i][%i]",i,j);
            pimpHist(hSegChi2Max[i][j],buffer,-2,20,vecRunId,mRunIdToName);
	    
	}
    }


    //End of Styling histograms###############################################################


    



    Physics->cd();

    /* hProtonMassMean->Write();
    hPipMassMean->Write();
    hPimMassMean->Write();
    hProtonMassMeanRpc->Write();
    hPipMassMeanRpc->Write();
    hPimMassMeanRpc->Write();
    hProtonMassMeanTof->Write();
    hPipMassMeanTof->Write();
    hPimMassMeanTof->Write();

    hTimeLepSumSys0Mean->Write();
    hTimeLepSumSys0Max->Write();
    hTimeLepSumSys0Sig->Write();
    hTimeLepSumSys1Mean->Write();
    hTimeLepSumSys1Max->Write();
    hTimeLepSumSys1Sig->Write();
    for(Int_t i=0; i<nStartStrips; i++)
    {
	hTimeLepStaSys0Mean[i]->Write();
	hTimeLepStaSys0Max[i]->Write();
	hTimeLepStaSys0Sig[i]->Write();
	hTimeLepStaSys1Mean[i]->Write();
	hTimeLepStaSys1Max[i]->Write();
	hTimeLepStaSys1Sig[i]->Write();
    }
    hTimePiSumSys0Mean->Write();
    hTimePiSumSys0Max->Write();
    hTimePiSumSys0Sig->Write();
    hTimePiSumSys1Mean->Write();
    hTimePiSumSys1Max->Write();
    hTimePiSumSys1Sig->Write();
    for(Int_t i=0; i<nStartStrips; i++)
    {
	hTimePiStaSys0Mean[i]->Write();
	hTimePiStaSys0Max[i]->Write();
	hTimePiStaSys0Sig[i]->Write();
	hTimePiStaSys1Mean[i]->Write();
	hTimePiStaSys1Max[i]->Write();
	hTimePiStaSys1Sig[i]->Write();
    }


    Start->cd();

    hHitTimeMean->Write();
    hHitTimeRms->Write();
    hHitToNoveto->Write();
    hVertexX->Write();
    hVertexY->Write();
    hVertexZ->Write();
    for(Int_t i=0; i<nStartStrips; i++)
    {
	hStCalTime[i]->Write();
	hVtCalTime[i]->Write();
    }


    RICH->cd();

    for(Int_t i=0; i<6; i++)
    {
	hAvChargeMax[i]->Write();
	hNumPadsMax[i]->Write();
	hNumRingsMax[i]->Write();
	hAvChargeMean[i]->Write();
	hNumPadsMean[i]->Write();
	hNumRingsMean[i]->Write();
	hHVMeanRich[i]->Write();
	hHVSigmaRich[i]->Write();
	hHVMinRich[i]->Write();
	hHVMaxRich[i]->Write();
	hHVNDataRich[i]->Write();
    } */



    MDC->cd();

    for(Int_t i=0; i<6; i++)
    {
	for(Int_t j=0; j<4; j++)
	{
	    for(Int_t k=0; k<16; k++)
	    {
		hMboCounts[i][j][k]->Write(); 
	    }
	    for(Int_t k=0; k<12; k++)
	    {
		hHVMeanMdc[i][j][k]->Write();         
		hHVSigmaMdc[i][j][k]->Write();        
		hHVMinMdc[i][j][k]->Write();          
		hHVMaxMdc[i][j][k]->Write();          
		hHVNDataMdc[i][j][k]->Write();        
	    }

	    for(Int_t k=0; k<6; k++)
	    {
		hLayerCounts[i][j][k]->Write();       
		hToTCalMean[i][j][k]->Write();        
		hToTCalMax[i][j][k]->Write();         
		hTime1CalMean[i][j][k]->Write();      
		hTime1CalMax[i][j][k]->Write();       
	    }
	    hCountsRaw[i][j]->Write();                
	    hCountsCal[i][j]->Write();                
	}
	for(Int_t j=0; j<2; j++)
	{
	    hSegCounter[i][j]->Write();               
	    hSegCounterFit[i][j]->Write();            
	    hSegChi2Mean[i][j]->Write();              
	    hSegChi2Max[i][j]->Write();               
	}
    }








    //End of filling  Histograms #######################################################
}

void AnalysisMacro::pimpHist(TH1F* hist,TString YTitle,Float_t y_low, Float_t y_up, vector<Int_t>& vec,map<Int_t,TString>& m)
{
    //hist->SetXTitle("Filenumber");
    hist->SetYTitle(YTitle);
    hist->GetYaxis()->SetRangeUser(y_low, y_up);
    hist->SetMarkerSize(1);
    hist->SetMarkerStyle(22);
    hist->SetMarkerColor(4);
    hist->SetLineColor(2);
    hist->SetOption("P");

    for(UInt_t i=0; i<vec.size(); i++)
    {
	hist->GetXaxis()->SetBinLabel(i+1,m[vec[i]]);
    }

    hist->SetLabelSize(0.015,"X");
    hist->GetXaxis()->LabelsOption("v");

}

void AnalysisMacro::finalize(TFile *outfile)
{
	outfile->Close();
	std::cout << "File " << outfile->GetPath() << " is finalized\n";
}

#endif  //#ifdef AnalysisMacro_cxx
