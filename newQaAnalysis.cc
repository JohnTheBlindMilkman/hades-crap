#include "Includes.h"
//#include "FemtoMixer/EventCandidate.hxx"
//#include "FemtoMixer/PairCandidate.hxx"
#include "../JJFemtoMixer/JJFemtoMixer.hxx"
#include "FemtoMixer/PairUtils.hxx"
#include "FemtoMixer/EventUtils.hxx"

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>

// defining a helper function
template <typename T>
bool isSim(T *t) {return false;}
template <>
bool isSim(HParticleCandSim *t) {return true;}

struct HistogramCollection
{
	TH2D hDPhiDThtaSign,hDPhiDThtaBckg;
};

int CalcBinDistance(TH2D *hist, double recoX, double recoY, double kineX, double kineY)
{
	int binxReco = hist->GetXaxis()->FindBin(recoX);
	int binyReco = hist->GetYaxis()->FindBin(recoY);
	int binxKine = hist->GetXaxis()->FindBin(kineX);
	int binyKine = hist->GetYaxis()->FindBin(kineY);

	return std::abs(binxReco - binxKine) + std::abs(binyReco - binyKine);
}

int newQaAnalysis(TString inputlist = "", TString outfile = "qaOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = 10)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);

	constexpr bool isCustomDst{false};
	constexpr bool isSimulation{true};
	constexpr int protonPID{14};

	//--------------------------------------------------------------------------------
    // Initialization of the global ROOT object and the Hades Loop
    // The Hades Loop used as an interface to the DST data (Basically a container of a TChain).
    // kTRUE - The global HADES object is being created if not existing
    //--------------------------------------------------------------------------------
    TROOT dst_analysis("DstAnalysisMacro", "Simple DST analysis Macro");
    HLoop* loop = new HLoop(kTRUE);
	TString beamtime="apr12";
	
	Int_t mdcMods[6][4]=
		{ {1,1,1,1},
		{1,1,1,1},
		{1,1,1,1},
		{1,1,1,1},
		{1,1,1,1},
		{1,1,1,1} };
	TString asciiParFile = "";
	TString rootParFile;
	if (isSimulation)
	{
		rootParFile = "/cvmfs/hadessoft.gsi.de/param/sim/apr12/allParam_APR12_sim_run_12001_gen9_07112017.root";
	}
	else
	{
		//rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/apr12/allParam_APR12_gen9_27092017.root";
		rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/apr12/allParam_APR12_gen10_16122024.root";
		//rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/feb24/allParam_feb24_gen0_16042024.root"; // Au+Au 800 MeV
	}
	TString paramSource = "root"; // root, ascii, oracle
	TString paramrelease = "APR12_dst_gen10";
	HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,shower,wall,start,tbox");
	HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease);  // now, APR12_gen2_dst

    //--------------------------------------------------------------------------------
    // The following block finds / adds the input DST files to the HLoop
    //--------------------------------------------------------------------------------
    if (maxFiles == -1)
		loop->addMultFiles(inputlist);     //use instead of addFiles if run on batch farm
    else 
    {
		Int_t nFiles = 0;


		TString inputFolder;
		if (isSimulation) // simulation
		{
			//inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 800 MeV
			//inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV gen9
			inputFolder = "/lustre/hades/dstsim/apr12/au1230au/gen10/bmax10/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV gen10
		}
		else // data
		{
			if (isCustomDst)
			{
				inputFolder = "/lustre/hades/user/kjedrzej/customDST/apr12PlusHMdcSeg/5sec/109"; // test Robert filtered events
			}
			else
			{
				//inputFolder = "/lustre/hades/dst/feb24/gen0c/060/01/root"; // Au+Au 800 MeV
				//inputFolder = "/lustre/hades/dst/apr12/gen9/122/root"; // Au+Au 2.4 GeV
				inputFolder = "/lustre/hades/dst/apr12/gen10/122/root"; // Au+Au 2.4 GeV
			}
		}
	
		TSystemDirectory* inputDir = new TSystemDirectory("inputDir", inputFolder);
		TList* files = inputDir->GetListOfFiles();

		for (Int_t i = 0; i <= files->LastIndex() && nFiles < maxFiles; i++) 
		{
			if (((TSystemFile*) files->At(i))->IsDirectory())
				continue;
			
			loop->addFile(inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName());
			nFiles++;
		}
	}

	loop->readSectorFileList("/lustre/hades/user/sspies/SectorFileLists/Apr12AuAu1230_Gen10_Hadrons.list");
    
    //--------------------------------------------------------------------------------
    // Booking the categories to be read from the DST files.
    // By default all categories are booked therefore -* (Unbook all) first and book the ones needed
    // All required categories have to be booked except the global Event Header which is always booked
    //--------------------------------------------------------------------------------
	std::string inputString = "-*,+HParticleCand,+HParticleEvtInfo,+HWallHit";
	if (isCustomDst)
		inputString += ",+HMdcSeg";
	if (isSimulation)
		inputString += ",+HGeantKine";
    if (!loop->setInput(inputString.data())) // make sure to use +HMdcSeg if you want the wires
		exit(1);

	gHades->setBeamTimeID(HADES::kApr12); // this is needed when using the ParticleEvtChara
	
    //--------------------------------------------------------------------------------
    // Setting the cache size of the HLoop internal TChain to read data from the DSTs
    // Improves performance of the lustre storage by decreasing load on lustre META servers
    //--------------------------------------------------------------------------------
    loop->getChain()->SetCacheSize(8e6); // 8Mb
    loop->getChain()->AddBranchToCache("*", kTRUE);
    loop->getChain()->StopCacheLearningPhase();

    loop->printCategories(); // Just for informative purposes
	
    //--------------------------------------------------------------------------------
    // Creating the placeholder variables to read data from categories and getting categories (They have to be booked!)
    //--------------------------------------------------------------------------------
    HParticleCandSim*    particle_cand;	// for simulation use HParticleCandSim*; for data HParticleCand*
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
	HCategory* kine_cand_cat = (HCategory*) HCategoryManager::getCategory(catGeantKine);

	if (isSimulation && !isSim(particle_cand)) // verification if you changed particle_cand class for running simulations
	{
		throw std::runtime_error("particle candidate must be of type HParticleCandSim"); // in C++17 this can be evaluated at compile-time, c++14 doesnt support if constexpr (condition)...
	}
	else if (!isSimulation && isSim(particle_cand))
	{
		throw std::runtime_error("particle candidate must be of type HParticleCand");
	}
    
    //if (!particle_cand_cat || !kine_cand_cat) // If the category for the reconstructed tracks does not exist, the macro makes no sense
	if (!particle_cand_cat)
		exit(1);
	
    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================

	constexpr float fBeamRapidity = 0.74f; // God I hope this is correct
	constexpr float fMeVtoGeV = 1.f/1000.f;

	TH1D *hZVertex = new TH1D("hZVertex","distribution of z component of the vertex",700,-65,5);
	TH1D *hXVertex = new TH1D("hXVertex","distribution of x component of the vertex",401,-20,20);
	TH1D *hYVertex = new TH1D("hYVertex","distribution of x component of the vertex",401,-20,20);
	TH1D *hXMom = new TH1D("hXMom","p_{x} distribution of accepted protons",3000,-1500,1500);
	TH1D *hYMom = new TH1D("hYMom","p_{y} distribution of accepted protons",3000,-1500,1500);
	TH1D *hZMom = new TH1D("hZMom","p_{z} distribution of accepted protons",3000,0,3000);
	TH1D *hEne = new TH1D("hEne","Energy distribution of accepted protons",2500,900,3400);
	TH2D *hBetaMomTof = new TH2D("hBetaMomTof","#beta vs p of accepted protons (ToF);p #times c [MeV/c];#beta",125,0,2500,100,0,1.);
	TH2D *hBetaMomRpc = new TH2D("hBetaMomRpc","#beta vs p of accepted protons (RPC);p #times c [MeV/c];#beta",125,0,2500,100,0,1.);
	TH2D *hPtRap = new TH2D("hPtRap","p_{T} vs y_{c.m} of accepted protons;p_{T} [MeV/c];y_{c.m.}",200,0,2000,121,-1.15,1.25);
	TH2D *hBetaMomTofReco = new TH2D("hBetaMomTofReco","#beta vs p of accepted and well reconstructed protons (ToF);p #times c [MeV/c];#beta",125,0,2500,100,0,1.);
	TH2D *hBetaMomRpcReco = new TH2D("hBetaMomRpcReco","#beta vs p of accepted and well reconstructed protons (RPC);p #times c [MeV/c];#beta",125,0,2500,100,0,1.);
	TH2D *hPtRapReco = new TH2D("hPtRapReco","p_{T} vs y_{c.m} of accepted and well reconstructed protons;p_{T} [MeV/c];y_{c.m.}",200,0,2000,121,-1.15,1.25);
	TH2D *hM2momTof = new TH2D("hM2momTof","m^{2} vs p of accepted protons (ToF);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0.4,1.6,1250,0,2.5);
	TH2D *hM2momRpc = new TH2D("hM2momRpc","m^{2} vs p of accepted protons (RPC);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0.4,1.6,1250,0,2.5);
	TH1D *hMinvTof = new TH1D("hMinvTof","m_{inv} of accepted protons (ToF);m_{inv} [GeV/c^{2}];N",600,0.4,1.6);
	TH1D *hMinvRpc = new TH1D("hMinvRpc","m_{inv} of accepted protons (RPC);m_{inv} [GeV/c^{2}];N",600,0.4,1.6);
	TH2D *hSegNcells = new TH2D("hSegNcells","MDC segment vs number of fired cells;seg;cells",24,0.5,24.5,10,-0.5,9.5);
	TH2D *hPhiTheta = new TH2D("hPhiTheta","Angular distribution of the tracks;#phi [deg];#theta [deg]",360,0,360,90,0,90);
	TH2D *hMetaCellsToF = new TH2D("hMetaCellsToF","Hit meta cells of ToF;Sector;Meta Cell",6,-0.5,5.5,64,0,64);
	TH2D *hMetaCellsRPC = new TH2D("hMetaCellsRPC","Hit meta cells of RPC;Sector;Meta Cell",6,-0.5,5.5,186,64,250);

	TH2D *hMomResolution = new TH2D("hMomResolution","Momentum difference of protons;p_{reco} [MeV/c];#frac{1}{p_{kin}} - #frac{1}{p_{reco}} [MeV/c]",80,0,4000,1000,-0.001, 0.001);
	TH2D *hPhiResolution = new TH2D("hPhiResolution","#phi angle difference of protons;p_{reco} [MeV/c];#phi_{kin} - #phi_{reco} [deg]",80,0,4000,500,-20, 20);
	TH2D *hThetaResolution = new TH2D("hThetaResolution","#theta angle difference of protons;p_{reco}  [MeV/c];#theta_{kin} - #theta_{reco} [deg]",80,0,4000,500,-10, 10);
	TH2D *hQinvResolution = new TH2D("hQinvResolution","q_{inv} discrepancy between ideal and recunstructed proton pairs;qinv_{reco};qinv_{kine}",750,0,3000,750,0,3000);
	TH2D *hQoutResolution = new TH2D("hQoutResolution","q_{out} discrepancy between ideal and recunstructed proton pairs;q_{out}^{reco};q_{out}^{kine}",64,0,500,64,0,500);
	TH2D *hQsideResolution = new TH2D("hQsideResolution","q_{side} discrepancy between ideal and recunstructed proton pairs;q_{side}^{reco};q_{side}^{kine}",64,0,500,64,0,500);
	TH2D *hQlongResolution = new TH2D("hQlongResolution","q_{long} discrepancy between ideal and recunstructed proton pairs;q_{long}^{reco};q_{long}^{kine}",64,0,500,64,0,500);

	TH2D *hInnerChi2Phi = new TH2D("hInnerChi2Phi","#chi^{2}_{inner} vs #phi angle difference of accepted protons;#phi_{kine} - #phi_{reco} [deg]; #chi^{2}_{inner}",500,-10, 10,300,0,30);
	TH2D *hInnerChi2Theta = new TH2D("hInnerChi2Theta","#chi^{2}_{inner} vs #theta angle difference of accepted protons;#theta_{kine} - #theta_{reco} [deg]; #chi^{2}_{inner}",500,-10, 10,300,0,30);
	TH2D *hOuterChi2Phi = new TH2D("hOuterChi2Phi","#chi^{2}_{outer} vs #phi angle difference of accepted protons;#phi_{kine} - #phi_{reco} [deg]; #chi^{2}_{outer}",500,-10, 10,300,0,30);
	TH2D *hOuterChi2Theta = new TH2D("hOuterChi2Theta","#chi^{2}_{outer} vs #theta angle difference of accepted protons;#theta_{kine} - #theta_{reco} [deg]; #chi^{2}_{outer}",500,-10, 10,300,0,30);
	TH2D *hMetaQualityMom = new TH2D("hMetaQualityMom","Q_{META} vs p difference of accepted protons;p_{kine} - p_{reco} [MeV/c];Q_{META}",500,-200,200,500,0,4);
	TH2D *hChi2Mom = new TH2D("hChi2Mom","#chi^{2}_{RK} vs p difference of accepted protons;p_{kine} - p_{reco} [MeV/c];#chi2^{2}_{RK}",500,-200,200,500,0,1200);

	TH2D *hQinvSLGood = new TH2D("hQinvSLGood","q_{inv} vs Splitting Level for signal of accepted p-p CF;q_{inv} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQinvSLBad = new TH2D("hQinvSLBad","q_{inv} vs Splitting Level for signal of rejected p-p CF;q_{inv} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQinvSWGood = new TH2D("hQinvSWGood","q_{inv} vs Shared Wires for signal of accepted p-p CF;q_{inv} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQinvSWBad = new TH2D("hQinvSWBad","q_{inv} vs Shared Wires for signal of rejected p-p CF;q_{inv} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQinvBLGood = new TH2D("hQinvBLGood","q_{inv} vs Shared layers for signal of accepted p-p CF;q_{inv} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQinvBLBad = new TH2D("hQinvBLBad","q_{inv} vs Shared layers for signal of rejected p-p CF;q_{inv} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQinvMWDGood = new TH2D("hQinvMWDGood","q_{inv} vs Minimal Wire DIstance for signal of accepted p-p CF;q_{inv} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQinvMWDBad = new TH2D("hQinvMWDBad","q_{inv} vs Minimal Wire DIstance for signal of rejected p-p CF;q_{inv} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQinvSMCGood = new TH2D("hQinvSMCGood","q_{inv} vs Shared Meta Cells for signal of accepted p-p CF;q_{inv} [MeV/c];SMC",250,0,3000,4,0,4);
	TH2D *hQinvSMCBad = new TH2D("hQinvSMCBad","q_{inv} vs Shared Meta Cells for signal of rejected p-p CF;q_{inv} [MeV/c];SMC",250,0,3000,4,0,4);
	TH2D *hWiresMultiplicityGood = new TH2D("hWiresMultiplicityGood","Wire multiplicity per event of accepted protons;Sector;Layer",6,0,6,24,0,24);
	TH2D *hDPhiDThetaSignGood = new TH2D("hDPhiDThetaSignGood","Signal of angular distribution of proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaBckgGood = new TH2D("hDPhiDThetaBckgGood","Backgound of angular distribution of proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaSignBad = new TH2D("hDPhiDThetaSignBad","Signal of angular distribution of accepted proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaBckgBad = new TH2D("hDPhiDThetaBckgBad","Backgound of angular distribution of rejected proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);

	TH2D *hQoutSLGood = new TH2D("hQoutSLGood","q_{out} vs Splitting Level for signal of accepted p-p CF;q_{out} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQoutSLBad = new TH2D("hQoutSLBad","q_{out} vs Splitting Level for signal of rejected p-p CF;q_{out} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQoutSWGood = new TH2D("hQoutSWGood","q_{out} vs Shared Wires for signal of accepted p-p CF;q_{out} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQoutSWBad = new TH2D("hQoutSWBad","q_{out} vs Shared Wires for signal of rejected p-p CF;q_{out} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQoutBLGood = new TH2D("hQoutBLGood","q_{out} vs Shared layers for signal of accepted p-p CF;q_{out} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQoutBLBad = new TH2D("hQoutBLBad","q_{out} vs Shared layers for signal of rejected p-p CF;q_{out} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQoutMWDGood = new TH2D("hQoutMWDGood","q_{out} vs Minimal Wire DIstance for signal of accepted p-p CF;q_{out} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQoutMWDBad = new TH2D("hQoutMWDBad","q_{out} vs Minimal Wire DIstance for signal of rejected p-p CF;q_{out} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQoutSMCGood = new TH2D("hQoutSMCGood","q_{out} vs Shared Meta Cells for signal of accepted p-p CF;q_{out} [MeV/c];SMC",250,0,3000,4,0,4);
	TH2D *hQoutSMCBad = new TH2D("hQoutSMCBad","q_{out} vs Shared Meta Cells for signal of rejected p-p CF;q_{out} [MeV/c];SMC",250,0,3000,4,0,4);

	TH2D *hQsideSLGood = new TH2D("hQsideSLGood","q_{side} vs Splitting Level for signal of accepted p-p CF;q_{side} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQsideSLBad = new TH2D("hQsideSLBad","q_{side} vs Splitting Level for signal of rejected p-p CF;q_{side} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQsideSWGood = new TH2D("hQsideSWGood","q_{side} vs Shared Wires for signal of accepted p-p CF;q_{side} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQsideSWBad = new TH2D("hQsideSWBad","q_{side} vs Shared Wires for signal of rejected p-p CF;q_{side} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQsideBLGood = new TH2D("hQsideBLGood","q_{side} vs Shared layers for signal of accepted p-p CF;q_{side} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQsideBLBad = new TH2D("hQsideBLBad","q_{side} vs Shared layers for signal of rejected p-p CF;q_{side} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQsideMWDGood = new TH2D("hQsideMWDGood","q_{side} vs Minimal Wire DIstance for signal of accepted p-p CF;q_{side} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQsideMWDBad = new TH2D("hQsideMWDBad","q_{side} vs Minimal Wire DIstance for signal of rejected p-p CF;q_{side} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQsideSMCGood = new TH2D("hQsideSMCGood","q_{side} vs Shared Meta Cells for signal of accepted p-p CF;q_{side} [MeV/c];SMC",250,0,3000,4,0,4);
	TH2D *hQsideSMCBad = new TH2D("hQsideSMCBad","q_{side} vs Shared Meta Cells for signal of rejected p-p CF;q_{side} [MeV/c];SMC",250,0,3000,4,0,4);

	TH2D *hQlongSLGood = new TH2D("hQlongSLGood","q_{long} vs Splitting Level for signal of accepted p-p CF;q_{long} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQlongSLBad = new TH2D("hQlongSLBad","q_{long} vs Splitting Level for signal of rejected p-p CF;q_{long} [MeV/c];SL",250,0,3000,101,-2,2);
	TH2D *hQlongSWGood = new TH2D("hQlongSWGood","q_{long} vs Shared Wires for signal of accepted p-p CF;q_{long} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQlongSWBad = new TH2D("hQlongSWBad","q_{long} vs Shared Wires for signal of rejected p-p CF;q_{long} [MeV/c];SW",250,0,3000,24,0,24);
	TH2D *hQlongBLGood = new TH2D("hQlongBLGood","q_{long} vs Shared layers for signal of accepted p-p CF;q_{long} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQlongBLBad = new TH2D("hQlongBLBad","q_{long} vs Shared layers for signal of rejected p-p CF;q_{long} [MeV/c];BL",250,0,3000,24,0,24);
	TH2D *hQlongMWDGood = new TH2D("hQlongMWDGood","q_{long} vs Minimal Wire DIstance for signal of accepted p-p CF;q_{long} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQlongMWDBad = new TH2D("hQlongMWDBad","q_{long} vs Minimal Wire DIstance for signal of rejected p-p CF;q_{long} [MeV/c];MWD",250,0,3000,100,0,100);
	TH2D *hQlongSMCGood = new TH2D("hQlongSMCGood","q_{long} vs Shared Meta Cells for signal of accepted p-p CF;q_{long} [MeV/c];SMC",250,0,3000,4,0,4);
	TH2D *hQlongSMCBad = new TH2D("hQlongSMCBad","q_{long} vs Shared Meta Cells for signal of rejected p-p CF;q_{long} [MeV/c];SMC",250,0,3000,4,0,4);

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");

	std::vector<TH1D*> ktDistGood,rapDistGood,ktDistBad,rapDistBad;
	std::vector<TH2D*> hKtRapGood, hKtRapBad;
	constexpr std::array<const char *,7> centString{"overflow","0-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %"};
	for (const int i : {0,1,2,3,4})
	{
		ktDistGood.push_back(new TH1D(TString::Format("ktDist_cent%d_good",i),TString::Format("k_{T} distribution of accepted pairs for %s centrality",centString.at(i)),1000,0,2500));
		rapDistGood.push_back(new TH1D(TString::Format("rapDist_cent%d_good",i),TString::Format("pair rapidity distribution of accepted pairs for %s centrality",centString.at(i)),1000,-2,2));
		ktDistBad.push_back(new TH1D(TString::Format("ktDist_cent%d_bad",i),TString::Format("k_{T} distribution of rejected pairs for %s centrality",centString.at(i)),1000,0,2500));
		rapDistBad.push_back(new TH1D(TString::Format("rapDist_cent%d_bad",i),TString::Format("pair rapidity distribution of rejected pairs for %s centrality",centString.at(i)),1000,-2,2));
		hKtRapGood.push_back(new TH2D(TString::Format("hKtRapGoodCent%d",i),"k_{T} vs y_{c.m} of accepted proton pairs;k_{T} [MeV/c];y_{c.m.}",1000,0,2500,1000,-2,2));
		hKtRapBad.push_back(new TH2D(TString::Format("hKtRapBadCent%d",i),"k_{T} vs y_{c.m} of rejected proton pairs;k_{T} [MeV/c];y_{c.m.}",1000,0,2500,1000,-2,2));
	}

	// create objects for particle selection
	std::shared_ptr<Selection::EventCandidate> fEvent;	
	std::shared_ptr<Selection::TrackCandidate> fTrack;
	HParticleWireInfo fWireInfo;
	HGeantHeader *geantHeader;
	std::size_t tracks;
	float betaReco, betaKine, momReco, momKine, rapReco, rapKine, ptReco, ptKine;
	int binDistance;

	std::map<std::string,std::vector<std::shared_ptr<Selection::PairCandidate> > > fSignMap, fBckgMap;	

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixer;
	mixer.SetMaxBufferSize(0);
	mixer.SetEventHashingFunction(Mixing::EventGrouping{}.MakeEventGroupingFunction());
	mixer.SetPairHashingFunction(Mixing::PairGrouping{}.MakePairGroupingFunction1D());
	mixer.SetPairCuttingFunction(Mixing::PairRejection{}.MakePairRejectionFunction());
	mixer.PrintSettings();

    //--------------------------------------------------------------------------------
    // The following counter histogram is used to gather some basic information on the analysis
    //--------------------------------------------------------------------------------
    enum Counters_e {
	cNumAllEvents      = 0,
	cNumSelectedEvents = 1,
	cNumAllTracks      = 2,
	cNumSelectedTracks = 3,
	cNumAllPairs       = 4,
	cNumSelectedPairs  = 5,
	cNumCounters       = 6
    };

    TH1D* hCounter = new TH1D("hCounter", "", cNumCounters, 0, cNumCounters);

    hCounter->GetXaxis()->SetBinLabel(1, "All Events");
    hCounter->GetXaxis()->SetBinLabel(2, "Selected Events");
    hCounter->GetXaxis()->SetBinLabel(3, "All Tracks");
    hCounter->GetXaxis()->SetBinLabel(4, "Selected Tracks");
	hCounter->GetXaxis()->SetBinLabel(5, "All Pairs");
	hCounter->GetXaxis()->SetBinLabel(6, "Selected Pairs");

	//--------------------------------------------------------------------------------
	// wire information w/o HMdcSeg class access
	//--------------------------------------------------------------------------------
	HTaskSet *masterTaskSet = gHades->getTaskSet("all");
    HParticleMetaMatcher* matcher = new HParticleMetaMatcher();
    matcher->setDebug();
	matcher->setUseEMC(kFALSE);
	matcher->setRunWireManager(false);
	if (isCustomDst)
		matcher->setUseSeg(kTRUE);
    masterTaskSet->add(matcher);

	//--------------------------------------------------------------------------------
	// Momentum corrected for energy loss (look-up table)
	//--------------------------------------------------------------------------------
    HEnergyLossCorrPar enLossCorr;
    enLossCorr.setDefaultPar(beamtime);

	//--------------------------------------------------------------------------------
	// event characteristic & reaction plane
	//--------------------------------------------------------------------------------
	HParticleEvtChara evtChara;

	std::cout << "HParticleEvtChara: reading input for energy 1.23A GeV... " << std::endl;
	TString ParameterfileCVMFS;
	if (isSimulation) // Simulation
	{
		ParameterfileCVMFS = "/cvmfs/hadessoft.gsi.de/param/eventchara/centrality_epcorr_sim_au1230au_gen9vertex_UrQMD_minbias_2019_04_pass0.root";
	}
	else // Data
	{
		//ParameterfileCVMFS = "/lustre/hades/user/bkardan/param/development/centrality_epcorr_feb24_au800au_1850A_gen0c_2024_04_pass10.root";  // Au+Au 800 MeV
		ParameterfileCVMFS = "/cvmfs/hadessoft.gsi.de/param/eventchara/centrality_epcorr_apr12_gen8_2019_02_pass30.root"; // Au+Au 2.4 GeV
	}
	
	if (!evtChara.setParameterFile(ParameterfileCVMFS))
	{
		std::cout << "Parameterfile not found !!! " << std::endl;
		return kFALSE;
	}

	if (!evtChara.init())
	{
		std::cout << "HParticleEvtChara not init!!! " << std::endl;
		return kFALSE;
	}
	
	Int_t eCentEstSP  = HParticleEvtChara::kSelectedParticleCand;
	Int_t eCentEst    = HParticleEvtChara::kTOFRPC;
	Int_t eCentClass1 = HParticleEvtChara::k10;
	Int_t eEPcorr     = HParticleEvtChara::kDefault;
	std::cout << "\t selected EPcorrection method is:  "  << evtChara.getStringEventPlaneCorrection(eEPcorr) << std::endl;

	std::cout << "EVTChara for TOF+RPC hits " << std::endl;
	evtChara.printCentralityClass(eCentEst, eCentClass1);

	std::cout << "EVTChara for the selected particles " << std::endl;
	evtChara.printCentralityClass(eCentEstSP, eCentClass1);

    //--------------------------------------------------------------------------------
    // Creating and initializing the track sorter and a simple stopwatch object
    //--------------------------------------------------------------------------------
    HParticleTrackSorter sorter;
    sorter.init();
    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //--------------------------------------------------------------------------------
    // The amount of events to be processed
    //--------------------------------------------------------------------------------
    Long64_t nEvents = loop->getEntries();
    if (nDesEvents >= 0 && nEvents > nDesEvents)
		nEvents = nDesEvents;

    //--------------------------------------------------------------------------------
    // The global event loop which loops over all events in the DST files added to HLoop
    // The loop breaks if the end is reached
    //--------------------------------------------------------------------------------
    for (Long64_t event = 0; event < nEvents; event++) 
    {
		if (loop->nextEvent(event) <= 0) 
		{
			std::cout << " Last events processed " << endl;
			break;
		}

		// TString tmp; // dummy variable; required by HLoop::isNewFile
		// if (loop->isNewFile(tmp) && !isSimulation)
		// {
		// 	if (!loop->goodSector(0) || !loop->goodSector(1) || !loop->goodSector(3) || !loop->goodSector(4) || !loop->goodSector(5)) // no sector 2 in Au+Au
		// 	{
		// 		event += loop->getTree()->GetEntries() - 1;
		// 		continue;
		// 	}
		// }

		hCounter->Fill(cNumAllEvents);

		//--------------------------------------------------------------------------------
		// Just the progress of the analysis
		//--------------------------------------------------------------------------------
		HTool::printProgress(event, nEvents, 1, "Analyzed events: ");

		//--------------------------------------------------------------------------------
		// Getting the amount of tracks (Particle Candidates), the global event header, the Particle event Info object and the reconstructed global event vertex
		//--------------------------------------------------------------------------------
		Int_t nTracks           = particle_cand_cat->getEntries();
		event_header            = gHades->getCurrentEvent()->getHeader();
		particle_info           = HCategoryManager::getObject(particle_info, particle_info_cat, 0);
		HGeomVector EventVertex  = event_header->getVertexReco().getPos();
		
		Int_t centClassIndex    = evtChara.getCentralityClass(eCentEst, eCentClass1); // 0 is overflow, 1 is 0-10, etc.
		float EventPlane = -1;
		float EventPlaneA = -1;
		float EventPlaneB = -1;

		if constexpr (isSimulation)
		{
			geantHeader = loop->getGeantHeader();
			if (geantHeader == nullptr)
				continue;
			
			EventPlane = geantHeader->getEventPlane() * TMath::DegToRad();
			EventPlaneA = EventPlane;
			EventPlaneB = EventPlane;
		}
		else
		{
			EventPlane = evtChara.getEventPlane(eEPcorr);
			EventPlaneA = evtChara.getEventPlane(eEPcorr,1);
			EventPlaneB = evtChara.getEventPlane(eEPcorr,2);
		}
		
		if (EventPlane < 0)
			continue;
		if (EventPlaneA < 0 || EventPlaneB < 0)
			continue;

		fEvent = std::make_shared<Selection::EventCandidate>(event_header,particle_info,centClassIndex,EventPlane);

		//--------------------------------------------------------------------------------
		// Discarding bad events with multiple criteria and counting amount of all / good events
		//--------------------------------------------------------------------------------
        
		// remove first two to get vortex x,y,z
		if (   !particle_info->isGoodEvent(Particle::kGoodVertexClust)
			|| !particle_info->isGoodEvent(Particle::kGoodVertexCand) 
			|| !particle_info->isGoodEvent(Particle::kGoodSTART)
			|| !particle_info->isGoodEvent(Particle::kNoPileUpSTART)
			|| !particle_info->isGoodEvent(Particle::kGoodTRIGGER)
			|| !particle_info->isGoodEvent(Particle::kNoVETO)
			|| !particle_info->isGoodEvent(Particle::kGoodSTARTVETO)
			|| !particle_info->isGoodEvent(Particle::kGoodSTARTMETA))
			continue;
	
		if (particle_info->getNStartCluster() >= 5)
			continue;

		//================================================================================================================================================================
		// Put your analyses on event level here
		//================================================================================================================================================================
		
		hXVertex->Fill(fEvent->GetX());
		hYVertex->Fill(fEvent->GetY());
		hZVertex->Fill(fEvent->GetZ());

		if (!fEvent->SelectEvent<HADES::Target::Setup::Apr12>({1,2,3,4},2,2,2))
			continue;
		
		hCounter->Fill(cNumSelectedEvents);

		//--------------------------------------------------------------------------------
		// Resetting the track sorter and selecting hadrons ranked by Chi2 Runge Kutta
		//--------------------------------------------------------------------------------
		sorter.cleanUp();
		sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE);
		sorter.fill(HParticleTrackSorter::selectHadrons);
		sorter.selectBest(Particle::ESwitch::kIsBestRKSorter, Particle::ESelect::kIsHadronSorter);
	
		//--------------------------------------------------------------------------------
		// The loop over all tracks (Particle Candidates in the current event
		//--------------------------------------------------------------------------------
		for (Int_t track = 0; track < nTracks; track++) 
		{
			particle_cand = HCategoryManager::getObject(particle_cand, particle_cand_cat, track);

			// I have a vague idea about how it should be done: set momentum and then call calc4vectorproperties before using
			particle_cand->setMomentum(particle_cand->getCorrectedMomentumPID(protonPID));

			matcher->getWireInfoDirect(particle_cand,fWireInfo);
			
			//--------------------------------------------------------------------------------
			// Discarding all tracks that have been discarded by the track sorter and counting all / good tracks
			//--------------------------------------------------------------------------------
			hCounter->Fill(cNumAllTracks);
	
			if (!particle_cand->isFlagBit(Particle::kIsUsed))
				continue;
	
			//--------------------------------------------------------------------------------
			// Getting information on the current track (Not all of them necessary for all analyses)
			//--------------------------------------------------------------------------------
			if constexpr (isSimulation)
			{
				fTrack = std::make_shared<Selection::TrackCandidate>(
					particle_cand,
					static_cast<HGeantKine*>(kine_cand_cat->getObject(particle_cand->getGeantTrack() - 1)),
					HADES::MDC::CreateTrackLayers(fWireInfo),
					fEvent->GetID(),
					fEvent->GetReactionPlane(),
					track,
					protonPID);
			}
			else
			{
				fTrack = std::make_shared<Selection::TrackCandidate>(
					particle_cand,
					HADES::MDC::CreateTrackLayers(fWireInfo),
					fEvent->GetID(),
					fEvent->GetReactionPlane(),
					track,
					protonPID);
			}
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fTrack->GetSystem() == Selection::Detector::RPC)
			{
				for (const auto &hit : fTrack->GetMetaHits())
				{
					hMetaCellsRPC->Fill(fTrack->GetSector(),hit);
				}
			}
			else
			{
				for (const auto &hit : fTrack->GetMetaHits())
				{
					hMetaCellsToF->Fill(fTrack->GetSector(),hit);
				}
			}

			if (!fTrack->SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				continue;

			betaReco = fTrack->GetBeta();
			momReco = fTrack->GetP();
			rapReco = fTrack->GetRapidity();
			ptReco = fTrack->GetPt();
			if (isSimulation)
			{
				betaKine = fTrack->GetGeantKine()->GetBeta();
				momKine = fTrack->GetGeantKine()->GetP();
				rapKine = fTrack->GetGeantKine()->GetRapidity();
				ptKine = fTrack->GetGeantKine()->GetPt();
			

				hPtRap->Fill(ptReco, rapReco - fBeamRapidity);
				binDistance = CalcBinDistance(hPtRap, ptReco, rapReco - fBeamRapidity, ptKine, rapKine - fBeamRapidity);
				hPtRapReco->Fill(ptReco, rapReco - fBeamRapidity, 1./(binDistance + 1));
			
				hPhiTheta->Fill(fTrack->GetPhi(),fTrack->GetTheta());
				hInnerChi2Phi->Fill(fTrack->GetGeantKine()->GetPhi() - fTrack->GetPhi(),fTrack->GetInnerSegChi2());
				hInnerChi2Theta->Fill(fTrack->GetGeantKine()->GetTheta() - fTrack->GetTheta(),fTrack->GetInnerSegChi2());
				hOuterChi2Phi->Fill(fTrack->GetGeantKine()->GetPhi() - fTrack->GetPhi(),fTrack->GetOuterSegChi2());
				hOuterChi2Theta->Fill(fTrack->GetGeantKine()->GetTheta() - fTrack->GetTheta(),fTrack->GetOuterSegChi2());
				hMetaQualityMom->Fill(momKine - momReco,fTrack->GetMetaMatchQuality());
				hChi2Mom->Fill(momKine - momReco, fTrack->GetChi2());

				hMomResolution->Fill(momReco,1./momKine - 1./momReco);
				hPhiResolution->Fill(momReco,fTrack->GetGeantKine()->GetPhi() - fTrack->GetPhi());
				hThetaResolution->Fill(momReco,fTrack->GetGeantKine()->GetTheta() - fTrack->GetTheta());
			}
			
			for (const int &layer : HADES::MDC::WireInfo::allLayerIndexing)
			{
				hSegNcells->Fill(layer + 1,fTrack->GetWires(layer).size());
				hWiresMultiplicityGood->Fill(fTrack->GetSector(), layer, fTrack->GetWires(layer).size());
			}

			hXMom->Fill(fTrack->GetPx());
			hYMom->Fill(fTrack->GetPy());
			hZMom->Fill(fTrack->GetPz());
			hEne->Fill(fTrack->GetEnergy());

			fEvent->AddTrack(fTrack);

			if (fTrack->GetSystem() == Selection::Detector::RPC)
			{
				hBetaMomRpc->Fill(momReco * fTrack->GetCharge(),betaReco);
				binDistance = CalcBinDistance(hBetaMomRpc,momReco * fTrack->GetCharge(),betaReco,momKine * fTrack->GetCharge(),betaKine);	
				hBetaMomRpcReco->Fill(momReco * fTrack->GetCharge(), betaReco, 1./(binDistance + 1));

				hM2momRpc->Fill(fTrack->GetM2()*fMeVtoGeV*fMeVtoGeV,abs(fTrack->GetP())*fMeVtoGeV);
				hMinvRpc->Fill(fTrack->GetM()*fMeVtoGeV);
			}
			else
			{
				hBetaMomTof->Fill(momReco * fTrack->GetCharge(),betaReco);
				binDistance = CalcBinDistance(hBetaMomTof,momReco * fTrack->GetCharge(),betaReco,momKine * fTrack->GetCharge(),betaKine);	
				hBetaMomTofReco->Fill(momReco * fTrack->GetCharge(), betaReco, 1./(binDistance + 1));

				hM2momTof->Fill(fTrack->GetM2()*fMeVtoGeV*fMeVtoGeV,abs(fTrack->GetP())*fMeVtoGeV);
				hMinvTof->Fill(fTrack->GetM()*fMeVtoGeV);
			}
		} // End of track loop

		tracks = fEvent->GetTrackListSize();
		hCounter->Fill(cNumSelectedTracks,tracks);

		if (tracks > 2)
		{
			fSignMap = mixer.AddEvent(fEvent,fEvent->GetTrackList());
			fBckgMap = mixer.GetSimilarPairs(fEvent);

			for (const auto &pair : fSignMap)
			{
				for (const auto &elem : pair.second)
				{
					float qOut,qSide,qLong,qInv;
					qInv = elem->GetQinv();
					std::tie(qOut,qSide,qLong) = elem->GetOSL();

					hCounter->Fill(cNumAllPairs);
					if (pair.first != "bad") // if not rejected
					{
						if (elem->AreTracksFromTheSameSector()) // those observables only make sense when we are in the same sector
						{
							hQinvSLGood->Fill(qInv,elem->GetSplittingLevel());
							hQoutSLGood->Fill(qOut,elem->GetSplittingLevel());
							hQsideSLGood->Fill(qSide,elem->GetSplittingLevel());
							hQlongSLGood->Fill(qLong,elem->GetSplittingLevel());

							hQinvSWGood->Fill(qInv,elem->GetSharedWires());
							hQoutSWGood->Fill(qOut,elem->GetSharedWires());
							hQsideSWGood->Fill(qSide,elem->GetSharedWires());
							hQlongSWGood->Fill(qLong,elem->GetSharedWires());

							hQinvBLGood->Fill(qInv,elem->GetBothLayers());
							hQoutBLGood->Fill(qOut,elem->GetBothLayers());
							hQsideBLGood->Fill(qSide,elem->GetBothLayers());
							hQlongBLGood->Fill(qLong,elem->GetBothLayers());

							if (elem->GetMinWireDistance().has_value) 
							{
								hQinvMWDGood->Fill(qInv,elem->GetMinWireDistance().value);
								hQoutMWDGood->Fill(qOut,elem->GetMinWireDistance().value);
								hQsideMWDGood->Fill(qSide,elem->GetMinWireDistance().value);
								hQlongMWDGood->Fill(qLong,elem->GetMinWireDistance().value);
							}

							hQinvSMCGood->Fill(qInv,elem->GetSharedMetaCells());
							hQoutSMCGood->Fill(qOut,elem->GetSharedMetaCells());
							hQsideSMCGood->Fill(qSide,elem->GetSharedMetaCells());
							hQlongSMCGood->Fill(qLong,elem->GetSharedMetaCells());
						}

						if (isSimulation)
						{
							hQinvResolution->Fill(qInv,elem->GetGeantKinePair()->GetQinv());
							float qOutKine,qSideKine,qLongKine;
							std::tie(qOutKine,qSideKine,qLongKine) = elem->GetGeantKinePair()->GetOSL();

							hQoutResolution->Fill(qOut,qOutKine);
							hQsideResolution->Fill(qSide,qSideKine);
							hQlongResolution->Fill(qLong,qLongKine);
						}

						hKtRapGood.at(fEvent->GetCentrality())->Fill(elem->GetKt(),elem->GetRapidity() - fBeamRapidity);

						ktDistGood.at(fEvent->GetCentrality())->Fill(elem->GetKt());
						rapDistGood.at(fEvent->GetCentrality())->Fill(elem->GetRapidity() - fBeamRapidity);
						if (pair.first != "0") // if is not overflow in pair kinematics
						{
							hDPhiDThetaSignGood->Fill(elem->GetDPhi(),elem->GetDTheta());
							hCounter->Fill(cNumSelectedPairs);
						}
					}
					else // if rejected
					{
						if (elem->AreTracksFromTheSameSector())
						{
							hQinvSLBad->Fill(qInv,elem->GetSplittingLevel());
							hQoutSLBad->Fill(qOut,elem->GetSplittingLevel());
							hQsideSLBad->Fill(qSide,elem->GetSplittingLevel());
							hQlongSLBad->Fill(qLong,elem->GetSplittingLevel());

							hQinvSWBad->Fill(qInv,elem->GetSharedWires());
							hQoutSWBad->Fill(qOut,elem->GetSharedWires());
							hQsideSWBad->Fill(qSide,elem->GetSharedWires());
							hQlongSWBad->Fill(qLong,elem->GetSharedWires());

							hQinvBLBad->Fill(qInv,elem->GetBothLayers());
							hQoutBLBad->Fill(qOut,elem->GetBothLayers());
							hQsideBLBad->Fill(qSide,elem->GetBothLayers());
							hQlongBLBad->Fill(qLong,elem->GetBothLayers());

							if (elem->GetMinWireDistance().has_value) 
							{
								hQinvMWDBad->Fill(qInv,elem->GetMinWireDistance().value);
								hQoutMWDBad->Fill(qOut,elem->GetMinWireDistance().value);
								hQsideMWDBad->Fill(qSide,elem->GetMinWireDistance().value);
								hQlongMWDBad->Fill(qLong,elem->GetMinWireDistance().value);
							}

							hQinvSMCBad->Fill(qInv,elem->GetSharedMetaCells());
							hQoutSMCBad->Fill(qOut,elem->GetSharedMetaCells());
							hQsideSMCBad->Fill(qSide,elem->GetSharedMetaCells());
							hQlongSMCBad->Fill(qLong,elem->GetSharedMetaCells());
						}

						hKtRapBad.at(fEvent->GetCentrality())->Fill(elem->GetKt(),elem->GetRapidity() - fBeamRapidity);
						
						ktDistBad.at(fEvent->GetCentrality())->Fill(elem->GetKt());
						rapDistBad.at(fEvent->GetCentrality())->Fill(elem->GetRapidity() - fBeamRapidity);	
						hDPhiDThetaSignBad->Fill(elem->GetDPhi(),elem->GetDTheta());
					}
				}
			}
			
			for (const auto &backgroundEntry : fBckgMap)
			{
				for (const auto &entry : backgroundEntry.second)
				{
					if (backgroundEntry.first != "bad" && backgroundEntry.first != "0")
					{
						hDPhiDThetaBckgGood->Fill(entry->GetDPhi(),entry->GetDTheta());
					}
					else
					{
						hDPhiDThetaBckgBad->Fill(entry->GetDPhi(),entry->GetDTheta());
					}
				}
			}
		}
	
	} // End of event loop
	 
	hWiresMultiplicityGood->Scale(1./nEvents);

	// printing some info about the RAM I'm using to know how much memory each job should be given
	static ProcInfo_t info;
	constexpr float toGB = 1.f/1024.f/1024.f;

	gSystem->GetProcInfo(&info);
	std::cout << "\n---=== Memory Usage ===---\n";
	std::cout << "resident memory used: " << info.fMemResident*toGB << " GB\t virtual memory used: " << info.fMemVirtual*toGB << " GB\n\n";

    //--------------------------------------------------------------------------------
    // Doing some cleanup and finalization work
    //--------------------------------------------------------------------------------
    sorter.finalize();
    timer.Stop();
    std::cout << "Finished DST processing" << endl;

    //--------------------------------------------------------------------------------
    // Creating output file and storing results there
    //--------------------------------------------------------------------------------
    TFile* out = new TFile(outfile.Data(), "RECREATE");
    out->cd();

    hCounter->Write();
	
    //================================================================================================================================================================
    // Remember to write your results to the output file here
    //================================================================================================================================================================

	hXVertex->Write();
	hYVertex->Write();
	hZVertex->Write();

	hXMom->Write();
	hYMom->Write();
	hZMom->Write();
	hEne->Write();
	hPtRap->Write();
	hBetaMomRpc->Write();
	hBetaMomTof->Write();
	hPtRapReco->Write();
	hBetaMomRpcReco->Write();
	hBetaMomTofReco->Write();
	hM2momRpc->Write();
	hM2momTof->Write();
	hMinvTof->Write();
	hMinvRpc->Write();
	hSegNcells->Write();
	hPhiTheta->Write();
	hMetaCellsToF->Write();
	hMetaCellsRPC->Write();

	hMomResolution->Write();
	hPhiResolution->Write();
	hThetaResolution->Write();
	hQinvResolution->Write();
	hQoutResolution->Write();
	hQsideResolution->Write();
	hQlongResolution->Write();

	hInnerChi2Phi->Write();
	hInnerChi2Theta->Write();
	hOuterChi2Phi->Write();
	hOuterChi2Theta->Write();
	hMetaQualityMom->Write();
	hChi2Mom->Write();

	hQinvSLGood->Write();
	hQinvSWGood->Write();
	hQinvBLGood->Write();
	hQinvMWDGood->Write();
	hQinvSMCGood->Write();
	hQinvSLBad->Write();
	hQinvSWBad->Write();
	hQinvBLBad->Write();
	hQinvMWDBad->Write();
	hQinvSMCBad->Write();

	hQoutSLGood->Write();
	hQoutSWGood->Write();
	hQoutBLGood->Write();
	hQoutMWDGood->Write();
	hQoutSMCGood->Write();
	hQoutSLBad->Write();
	hQoutSWBad->Write();
	hQoutBLBad->Write();
	hQoutMWDBad->Write();
	hQoutSMCBad->Write();

	hQsideSLGood->Write();
	hQsideSWGood->Write();
	hQsideBLGood->Write();
	hQsideMWDGood->Write();
	hQsideSMCGood->Write();
	hQsideSLBad->Write();
	hQsideSWBad->Write();
	hQsideBLBad->Write();
	hQsideMWDBad->Write();
	hQsideSMCBad->Write();

	hQlongSLGood->Write();
	hQlongSWGood->Write();
	hQlongBLGood->Write();
	hQlongMWDGood->Write();
	hQlongSMCGood->Write();
	hQlongSLBad->Write();
	hQlongSWBad->Write();
	hQlongBLBad->Write();
	hQlongMWDBad->Write();
	hQlongSMCBad->Write();

	hWiresMultiplicityGood->Write();
	hDPhiDThetaSignGood->Write();
	hDPhiDThetaBckgGood->Write();
	hDPhiDThetaSignBad->Write();
	hDPhiDThetaBckgBad->Write();

	for (const int i : {0,1,2,3,4})
	{
		hKtRapGood.at(i)->Write();
		hKtRapBad.at(i)->Write();
		ktDistGood.at(i)->Write();
		rapDistGood.at(i)->Write();
		ktDistBad.at(i)->Write();
		rapDistBad.at(i)->Write();
	}

    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    std::cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

