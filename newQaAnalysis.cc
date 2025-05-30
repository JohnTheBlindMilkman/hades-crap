#include "Includes.h"
#include "FemtoMixer/EventCandidate.hxx"
#include "FemtoMixer/PairCandidate.hxx"
#include "../JJFemtoMixer/JJFemtoMixer.hxx"
#include "../HFiredWires/HFiredWires.hxx"

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

std::string EventHashing(const std::shared_ptr<Selection::EventCandidate> &evt)
{
    //return JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetCentrality()),2) + JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
	return JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetNCharged()/10),2) + 
	JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetReactionPlane()/10),2) + 
	JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
}

std::string TrackHashing(const std::shared_ptr<Selection::TrackCandidate> &trck)
{
	//return std::to_string(trck->GetPt()/1000) + std::to_string(trck->GetRapidity());
	return std::to_string(static_cast<std::size_t>(std::abs(trck->GetPx()/100))) + 
		std::to_string(static_cast<std::size_t>(std::abs(trck->GetPy()/100))) + 
		std::to_string(static_cast<std::size_t>(std::abs(trck->GetPz()/100)));
}

std::string PairHashing(const std::shared_ptr<Selection::PairCandidate> &pair)
{
	// collection of slices: (array[n],array[n+1]>
    // make sure this is ordered!
    constexpr std::array<float,11> ktArr{0,200,400,600,800,1000,1200,1400,1600,1800,2000};
    constexpr std::array<float,14> yArr{0.09,0.19,0.29,0.39,0.49,0.59,0.69,0.79,0.89,0.99,1.09,1.19,1.29,1.39};
    constexpr std::array<float,9> EpArr{-202.5,-157.5,-112.5,-67.5,-22.5,22.5,67.5,112.5,157.5};

    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair->GetKt()) - ktArr.begin();
    std::size_t yCut = std::lower_bound(yArr.begin(),yArr.end(),pair->GetRapidity()) - yArr.begin();
    std::size_t EpCut = std::lower_bound(EpArr.begin(),EpArr.end(),pair->GetPhi()) - EpArr.begin();

	// reject if value is below first slice or above the last
	if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > yArr.size()-1 || EpCut == 0 || EpCut > EpArr.size()-1)
		return "0";
	else
    	return std::to_string((ktCut)) + std::to_string(yCut) + std::to_string(EpCut);
}

bool PairCut(const std::shared_ptr<Selection::PairCandidate> &pair)
{
	using Behaviour = Selection::PairCandidate::Behaviour;

	if (pair->AreTracksFromTheSameSector())
	{
		return pair->RejectPairByCloseHits<Behaviour::OneUnder>(0.7,2) ||
			pair->GetBothLayers() < 20 ||
			pair->GetSharedMetaCells() > 0;
	}
	else
	{
		return false;
	}
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
	HGeantKine*       geant_kine;
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
	HCategory* kine_cand_cat = (HCategory*) HCategoryManager::getCategory(catGeantKine);
	//HCategory *mdc_seg_cat = (HCategory*) gHades->getCurrentEvent()->getCategory(catMdcSeg); // so... this is my stupid fix for now, without it HFiredWires doesn't work

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

	//TH1D *hFemtoMixerTest = new TH1D("hFemtoMixerTest","",999999,0,999999);

	TH1D *hZVertex = new TH1D("hZVertex","distribution of z component of the vertex",700,-65,5);
	TH1D *hXVertex = new TH1D("hXVertex","distribution of x component of the vertex",401,-20,20);
	TH1D *hYVertex = new TH1D("hYVertex","distribution of x component of the vertex",401,-20,20);
	TH2D *hBetaMomTof = new TH2D("hBetaMomTof","#beta vs p of accepted protons (ToF);p #times c [MeV/c];#beta",1250,0,2500,200,0,1.);
	TH2D *hBetaMomRpc = new TH2D("hBetaMomRpc","#beta vs p of accepted protons (RPC);p #times c [MeV/c];#beta",1250,0,2500,200,0,1.);
	TH2D *hPtRap = new TH2D("hPtRap","p_{T} vs y_{c.m} of accepted protons;p_{T} [MeV/c];y_{c.m.}",2000,0,2000,121,-1.15,1.25);
	TH2D *hM2momTof = new TH2D("hM2momTof","m^{2} vs p of accepted protons (ToF);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0.4,1.6,1250,0,2.5);
	TH2D *hM2momRpc = new TH2D("hM2momRpc","m^{2} vs p of accepted protons (RPC);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0.4,1.6,1250,0,2.5);
	TH1D *hMinvTof = new TH1D("hMinvTof","m_{inv} of accepted protons (ToF);m_{inv} [GeV/c^{2}];N",600,0.4,1.6);
	TH1D *hMinvRpc = new TH1D("hMinvRpc","m_{inv} of accepted protons (RPC);m_{inv} [GeV/c^{2}];N",600,0.4,1.6);
	TH2D *hSegNcells = new TH2D("hSegNcells","MDC segment vs number of fired cells;seg;cells",24,0.5,24.5,10,-0.5,9.5);
	TH2D *hPhiTheta = new TH2D("hPhiTheta","Angular distribution of the tracks;#phi [deg];#theta [deg]",360,0,360,90,0,90);
	TH2D *hMetaCellsToF = new TH2D("hMetaCellsToF","Hit meta cells of ToF;Sector;Meta Cell",6,0,6,32,0,32);
	TH2D *hMetaCellsRPC = new TH2D("hMetaCellsRPC","Hit meta cells of RPC;Sector;Meta Cell",6,0,6,32,64,96);

	TH2D *hMomResolution = new TH2D("hMomResolution","Momentum difference of protons;p_{reco} [MeV/c];p_{kin} - p_{reco}  [MeV/c]",80,0,4000,500,-500, 500);
	TH2D *hPhiResolution = new TH2D("hPhiResolution","#phi angle difference of protons;p_{reco} [MeV/c];#phi_{kin} - #phi_{reco} [deg]",80,0,4000,500,-20, 20);
	TH2D *hThetaResolution = new TH2D("hThetaResolution","#theta angle difference of protons;p_{reco}  [MeV/c];#theta_{kin} - #theta_{reco} [deg]",80,0,4000,500,-10, 10);
	TH2D *hQinvResolution = new TH2D("hQinvResolution","q_{inv} discrepancy between ideal and recunstructed proton pairs;qinv_{reco};qinv_{kine}",750,0,3000,750,0,3000);

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
	TH2D *hKtRapGood = new TH2D("hKtRapGood","k_{T} vs y_{c.m} of accepted proton pairs;k_{T} [MeV/c];y_{c.m.}",2000,0,2000,201,-2,2);
	TH2D *hKtRapBad = new TH2D("hKtRapBad","k_{T} vs y_{c.m} of rejected proton pairs;k_{T} [MeV/c];y_{c.m.}",2000,0,2000,201,-2,2);
	TH2D *hWiresMultiplicityGood = new TH2D("hWiresMultiplicityGood","Wire multiplicity per event of accepted protons;Sector;Layer",6,0,6,24,0,24);
	TH2D *hDPhiDThetaSignGood = new TH2D("hDPhiDThetaSignGood","Signal of angular distribution of proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaBckgGood = new TH2D("hDPhiDThetaBckgGood","Backgound of angular distribution of proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaSignBad = new TH2D("hDPhiDThetaSignBad","Signal of angular distribution of accepted proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);
	TH2D *hDPhiDThetaBckgBad = new TH2D("hDPhiDThetaBckgBad","Backgound of angular distribution of rejected proton pairs;#Delta #phi [deg];#Delta #theta [deg]",721,-360,360,181,-90,90);

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	std::vector<TH1D*> ktDistGood,rapDistGood,ktDistBad,rapDistBad;
	constexpr std::array<const char *,7> centString{"overflow","0-10 %","10-20 %","20-30 %","30-40 %","40-50 %","50-60 %"};
	for (const int i : {0,1,2,3,4,5,6})
	{
		ktDistGood.push_back(new TH1D(TString::Format("ktDist_cent%d_good",i),TString::Format("k_{T} distribution of accepted pairs for %s centrality",centString.at(i)),1000,0,2000));
		rapDistGood.push_back(new TH1D(TString::Format("rapDist_cent%d_good",i),TString::Format("pair rapidity distribution of accepted pairs for %s centrality",centString.at(i)),1000,-2,2));
		ktDistBad.push_back(new TH1D(TString::Format("ktDist_cent%d_bad",i),TString::Format("k_{T} distribution of rejected pairs for %s centrality",centString.at(i)),1000,0,2000));
		rapDistBad.push_back(new TH1D(TString::Format("rapDist_cent%d_bad",i),TString::Format("pair rapidity distribution of rejected pairs for %s centrality",centString.at(i)),1000,-2,2));
	}

	// create objects for particle selection
	std::shared_ptr<Selection::EventCandidate> fEvent;	
	std::shared_ptr<Selection::TrackCandidate> fTrack;
	HParticleWireInfo fWireInfo;
	HGeantHeader *geantHeader;
	std::size_t tracks;

	std::map<std::string,std::vector<std::shared_ptr<Selection::PairCandidate> > > fSignMap, fBckgMap;	

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixer;
	mixer.SetMaxBufferSize(200);
	mixer.SetEventHashingFunction(EventHashing);
	//mixer.SetTrackHashingFunction(TrackHashing);
	mixer.SetPairHashingFunction(PairHashing);
	mixer.SetPairCuttingFunction(PairCut);
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

		TString tmp; // dummy variable; required by HLoop::isNewFile
		if (loop->isNewFile(tmp) && !isSimulation)
		{
			if (!loop->goodSector(0) || !loop->goodSector(1) || !loop->goodSector(3) || !loop->goodSector(4) || !loop->goodSector(5)) // no sector 2 in Au+Au
			{
				event += loop->getTree()->GetEntries() - 1;
				continue;
			}
		}

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
			matcher->getWireInfoDirect(particle_cand,fWireInfo);

			// for Feb24 there is no ene loss correction (yet?)
			particle_cand->setMomentum(enLossCorr.getCorrMom(protonPID,particle_cand->getMomentum(),particle_cand->getTheta()));
			
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
					Selection::TrackCandidate::CreateWireArray(fWireInfo),
					fEvent->GetID(),
					fEvent->GetReactionPlane(),
					track,
					protonPID);
			}
			else
			{
				fTrack = std::make_shared<Selection::TrackCandidate>(
					particle_cand,
					Selection::TrackCandidate::CreateWireArray(fWireInfo),
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
				// for (const auto &hit : fTrack->GetMetaHits())
				// {
				// 	hMetaCellsRPC->Fill(fTrack->GetSector(),hit);
				// }
			}
			else
			{
				// for (const auto &hit : fTrack->GetMetaHits())
				// {
				// 	hMetaCellsToF->Fill(fTrack->GetSector(),hit);
				// }
			}

			if (!fTrack->SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				continue;

			hPtRap->Fill(fTrack->GetPt(), fTrack->GetRapidity() - fBeamRapidity);
			hPhiTheta->Fill(fTrack->GetPhi(),fTrack->GetTheta());
			//hFemtoMixerTest->Fill(mixer.GetTrackHash(fTrack));

			geant_kine = static_cast<HGeantKine*>(kine_cand_cat->getObject(particle_cand->getGeantTrack() - 1));
			double kineMom = geant_kine->getTotalMomentum();
			double kinePhi = geant_kine->getPhiDeg();
			double kineTheta = geant_kine->getThetaDeg();

			hMomResolution->Fill(fTrack->GetP(),kineMom - fTrack->GetP());
			hPhiResolution->Fill(fTrack->GetP(),kinePhi - fTrack->GetPhi());
			hThetaResolution->Fill(fTrack->GetP(),kineTheta - fTrack->GetTheta());
			
			for (const int &layer : HADES::MDC::WireInfo::allLayerIndexing)
			{
				hSegNcells->Fill(layer + 1,fTrack->GetWires(layer).size());
				hWiresMultiplicityGood->Fill(fTrack->GetSector(), layer, fTrack->GetWires(layer).size());
			}

			fEvent->AddTrack(fTrack);

			if (fTrack->GetSystem() == Selection::Detector::RPC)
			{
				hBetaMomRpc->Fill(fTrack->GetP()*fTrack->GetCharge(),fTrack->GetBeta());
				hM2momRpc->Fill(fTrack->GetM2()*fMeVtoGeV*fMeVtoGeV,abs(fTrack->GetP())*fMeVtoGeV);
				hMinvRpc->Fill(fTrack->GetM()*fMeVtoGeV);
			}
			else
			{
				hBetaMomTof->Fill(fTrack->GetP()*fTrack->GetCharge(),fTrack->GetBeta());
				hM2momTof->Fill(fTrack->GetM2()*fMeVtoGeV*fMeVtoGeV,abs(fTrack->GetP())*fMeVtoGeV);
				hMinvTof->Fill(fTrack->GetM()*fMeVtoGeV);
			}
		} // End of track loop

		tracks = fEvent->GetTrackListSize();
		hCounter->Fill(cNumSelectedTracks,tracks);

		if (tracks > 2)
		{
			hCounter->Fill(cNumAllPairs,tracks*(tracks-1)/2); // all combinations w/o repetitions

			fSignMap = mixer.AddEvent(fEvent,fEvent->GetTrackList());
			fBckgMap = mixer.GetSimilarPairs(fEvent);

			for (const auto &pair : fSignMap)
			{
				for (const auto &elem : pair.second)
				{
					if (pair.first != "bad") // if not rejected
					{
						hQinvSLGood->Fill(elem->GetQinv(),elem->GetSplittingLevel());
						hQinvSWGood->Fill(elem->GetQinv(),elem->GetSharedWires());
						hQinvBLGood->Fill(elem->GetQinv(),elem->GetBothLayers());
						hQinvMWDGood->Fill(elem->GetQinv(),elem->GetMinWireDistance());
						hQinvSMCGood->Fill(elem->GetQinv(),elem->GetSharedMetaCells());

						if (isSimulation)
							hQinvResolution->Fill(elem->GetQinv(),elem->GetGeantKinePair()->GetQinv());

						hKtRapGood->Fill(elem->GetKt(),elem->GetRapidity() - fBeamRapidity);

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
						hQinvSLBad->Fill(elem->GetQinv(),elem->GetSplittingLevel());
						hQinvSWBad->Fill(elem->GetQinv(),elem->GetSharedWires());
						hQinvBLBad->Fill(elem->GetQinv(),elem->GetBothLayers());
						hQinvMWDBad->Fill(elem->GetQinv(),elem->GetMinWireDistance());
						hQinvSMCBad->Fill(elem->GetQinv(),elem->GetSharedMetaCells());

						hKtRapBad->Fill(elem->GetKt(),elem->GetRapidity() - fBeamRapidity);
						
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
					//std::cout << "Bckg pair: " << entry->GetID() << "\t SW: " << entry->GetSharedWires() << "\t BL: " << entry->GetBothLayers() << "\t SMC: " << entry->GetSharedMetaCells() << "\n";
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

	//hFemtoMixerTest->Write();
	hPtRap->Write();
	hBetaMomRpc->Write();
	hBetaMomTof->Write();
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

	hQinvSLGood->Write();
	hQinvSWGood->Write();
	hQinvBLGood->Write();
	hQinvMWDGood->Write();
	hQinvSMCGood->Write();
	hKtRapGood->Write();
	hQinvSLBad->Write();
	hQinvSWBad->Write();
	hQinvBLBad->Write();
	hQinvMWDBad->Write();
	hQinvSMCBad->Write();
	hKtRapBad->Write();
	hWiresMultiplicityGood->Write();
	hDPhiDThetaSignGood->Write();
	hDPhiDThetaBckgGood->Write();
	hDPhiDThetaSignBad->Write();
	hDPhiDThetaBckgBad->Write();

	for (const int i : {0,1,2,3,4,5,6})
	{
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

