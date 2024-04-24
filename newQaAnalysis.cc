#include "Includes.h"
#include "FemtoMixer/EventCandidate.hxx"
#include "FemtoMixer/PairCandidate.hxx"
#include "FemtoMixer/JJFemtoMixer.hxx"
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
	TH1D hQinvSign,hQinvBckg;
	TH2D hDphiDthetaSign,hDphiDthetaBckg;
	TH3D hQoslSign,hQoslBckg;
};

std::size_t EventHashing(const Selection::EventCandidate &evt)
{
    return evt.Centrality*1e1 + evt.TargetPlate;
}

std::size_t PairHashing(const Selection::PairCandidate &pair)
{
	// collection of slices: (array[n],array[n+1]>
    // make sure this is ordered!
    constexpr std::array<float,6> ktArr{150,450,750,1050,1350,1650};
    constexpr std::array<float,4> yArr{0,0.49,0.99,1.49};
    constexpr std::array<float,9> EpArr{-202.5,-157.5,-112.5,-67.5,-22.5,22.5,67.5,112.5,157.5};

    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair.Kt) - ktArr.begin();
    std::size_t yCut = std::lower_bound(yArr.begin(),yArr.end(),pair.Rapidity) - yArr.begin();
    std::size_t EpCut = std::lower_bound(EpArr.begin(),EpArr.end(),pair.AzimuthalAngle) - EpArr.begin();

	// reject if value is below first slice or above the last
	if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > yArr.size()-1 || EpCut == 0 || EpCut > EpArr.size()-1 /* || (abs(pair.DeltaPhi) < 8 && abs(pair.DeltaTheta) < 2) */)
		return 0;
	else
    	return ktCut*1e2 + yCut*1e1 + EpCut;
}

int newQaAnalysis(TString inputlist = "", TString outfile = "qaOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = 10)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);

	constexpr bool isSimulation = false; // for now this could be easly just const

	//--------------------------------------------------------------------------------
    // Initialization of the global ROOT object and the Hades Loop
    // The Hades Loop used as an interface to the DST data (Basically a container of a TChain).
    // kTRUE - The global HADES object is being created if not existing
    //--------------------------------------------------------------------------------
    TROOT dst_analysis("DstAnalysisMacro", "Simple DST analysis Macro");
    HLoop* loop = new HLoop(kTRUE);

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
			inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 800 MeV
			//inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV
		}
		else // data
		{
			inputFolder = "/lustre/hades/dst/feb24/gen0c/039/01/root"; // Au+Au 800 MeV
			//inputFolder = "/lustre/hades/dst/apr12/gen9/122/root"; // Au+Au 2.4 GeV
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
    
    //--------------------------------------------------------------------------------
    // Booking the categories to be read from the DST files.
    // By default all categories are booked therefore -* (Unbook all) first and book the ones needed
    // All required categories have to be booked except the global Event Header which is always booked
    //--------------------------------------------------------------------------------
    if (!loop->setInput("-*,+HGeantKine,+HParticleCand,+HParticleEvtInfo,+HWallHit"))
		exit(1);

	gHades->setBeamTimeID(HADES::kFeb24); // this is needed when using the ParticleEvtChara
	
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
    HParticleCand*    particle_cand;	// for simulation use HParticleCandSim*; for data HParticleCand*
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);

	if (isSimulation && !isSim(particle_cand)) // verification if you changed particle_cand class for running simulations
	{
		throw std::runtime_error("particle candidate must be of type HParticleCandSim"); // in C++17 this can be evaluated at compile-time, c++14 doesnt support if constexpr (condition)...
	}
	else if (!isSimulation && isSim(particle_cand))
	{
		throw std::runtime_error("particle candidate must be of type HParticleCand");
	}
    
    if (!particle_cand_cat) // If the category for the reconstructed trackes does not exist the macro makes no sense
		exit(1);
	
    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================

	constexpr float fBeamRapidity = 0.74; // God I hope this is correct
	constexpr float fMeVtoGeV = 1.f/1000.f;

	TH2D *hBetaMomTof = new TH2D("hBetaMomTof","#beta vs p of accepted protons (ToF);p #times c [MeV/c];#beta",1250,0,2500,200,0,1.);
	TH2D *hBetaMomRpc = new TH2D("hBetaMomRpc","#beta vs p of accepted protons (RPC);p #times c [MeV/c];#beta",1250,0,2500,200,0,1.);
	TH2D *hPtRap = new TH2D("hPtRap","p_{T} vs y_{c.m} of accepted protons;p_{T} [MeV/c];y_{c.m.}",2000,0,2000,121,-1.15,1.25);
	TH2D *hM2momTof = new TH2D("hM2momTof","m^{2} vs p of accepted protons (ToF);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0,1.2,1250,0,2.5);
	TH2D *hM2momRpc = new TH2D("hM2momRpc","m^{2} vs p of accepted protons (RPC);m^{2} [GeV^{2}/c^{4}];p [GeV/c]",600,0,1.2,1250,0,2.5);

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	// create objects for particle selection
	Selection::EventCandidate fEvent;	
	Selection::TrackCandidate fTrack;

    //--------------------------------------------------------------------------------
    // The following counter histogram is used to gather some basic information on the analysis
    //--------------------------------------------------------------------------------
    enum Counters_e {
	cNumAllEvents      = 0,
	cNumSelectedEvents = 1,
	cNumAllTracks      = 2,
	cNumSelectedTracks = 3,
	cNumCounters       = 4
    };

    TH1D* hCounter = new TH1D("hCounter", "", cNumCounters, 0, cNumCounters);

    hCounter->GetXaxis()->SetBinLabel(1, "All Events");
    hCounter->GetXaxis()->SetBinLabel(2, "Selected Events");
    hCounter->GetXaxis()->SetBinLabel(3, "All Tracks");
    hCounter->GetXaxis()->SetBinLabel(4, "Selected Tracks");

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
		ParameterfileCVMFS = "/lustre/hades/user/bkardan/param/development/centrality_epcorr_feb24_au800au_1850A_gen0c_2024_04_pass10.root";  // Au+Au 800 MeV
		//ParameterfileCVMFS = "/cvmfs/hadessoft.gsi.de/param/eventchara/centrality_epcorr_apr12_gen8_2019_02_pass30.root"; // Au+Au 2.4 GeV
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
			cout << " Last events processed " << endl;
			break;
		}

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
		
		Float_t vertZ = EventVertex.Z();
		Float_t vertX = EventVertex.X();
		Float_t vertY = EventVertex.Y();
		Int_t centClassIndex    = evtChara.getCentralityClass(eCentEst, eCentClass1); // 0 is overflow, 1 is 0-10, etc.
		Float_t EventPlane = 0;
		/* Float_t EventPlane = evtChara.getEventPlane(eEPcorr);
		Float_t EventPlaneA = evtChara.getEventPlane(eEPcorr,1);
		Float_t EventPlaneB = evtChara.getEventPlane(eEPcorr,2); */

		/* if (!isSimulation) // EP is reconstructed only for real data
		{
			if (EventPlane < 0)
				continue;
			if (EventPlaneA < 0 || EventPlaneB < 0)
				continue;
		} */

		Selection::CreateEvent(fEvent,vertX,vertY,vertZ,centClassIndex,EventPlane);

		//--------------------------------------------------------------------------------
		// Discarding bad events with multiple criteria and counting amount of all / good events
		//--------------------------------------------------------------------------------
        hCounter->Fill(cNumAllEvents);
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
	
		//================================================================================================================================================================
		// Put your analyses on event level here
		//================================================================================================================================================================
		
		if (!Selection::SelectEvent(fEvent)) 
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
			
			//--------------------------------------------------------------------------------
			// Discarding all tracks that have been discarded by the track sorter and counting all / good tracks
			//--------------------------------------------------------------------------------
			hCounter->Fill(cNumAllTracks);
	
			if (!particle_cand->isFlagBit(Particle::kIsUsed))
				continue;
	
			//--------------------------------------------------------------------------------
			// Getting information on the current track (Not all of them necessary for all analyses)
			//--------------------------------------------------------------------------------
			Selection::CreateTrack(fTrack,particle_cand);
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fTrack.System == Selection::Detector::RPC)
			{
				// fill RPC monitors for all tracks
			}
			else
			{
				// fill ToF monitors for all tracks
			}

			if (!Selection::SelectTrack(fTrack,betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				continue;

			hCounter->Fill(cNumSelectedTracks);
			hPtRap->Fill(fTrack.TransverseMomentum, fTrack.Rapidity - fBeamRapidity);

			fEvent.trackList.push_back(fTrack);

			if (fTrack.System == Selection::Detector::RPC)
			{
				hBetaMomRpc->Fill(fTrack.TotalMomentum*fTrack.Charge,fTrack.Beta);
				hM2momRpc->Fill(fTrack.Mass2*fMeVtoGeV*fMeVtoGeV,abs(fTrack.TotalMomentum)*fMeVtoGeV);
			}
			else
			{
				hBetaMomTof->Fill(fTrack.TotalMomentum*fTrack.Charge,fTrack.Beta);
				hM2momTof->Fill(fTrack.Mass2*fMeVtoGeV*fMeVtoGeV,abs(fTrack.TotalMomentum)*fMeVtoGeV);
			}

		} // End of track loop
	} // End of event loop
	
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
    cout << "Finished DST processing" << endl;

    //--------------------------------------------------------------------------------
    // Creating output file and storing results there
    //--------------------------------------------------------------------------------
    TFile* out = new TFile(outfile.Data(), "RECREATE");
    out->cd();

    hCounter->Write();
	
    //================================================================================================================================================================
    // Remember to write your results to the output file here
    //================================================================================================================================================================

	hPtRap->Write();
	hBetaMomRpc->Write();
	hBetaMomTof->Write();
	hM2momRpc->Write();
	hM2momTof->Write();

    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

