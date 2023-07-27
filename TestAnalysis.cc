#include "Includes.h"
#include <iostream>
#include <string>
#include <vector>

int TestAnalysis(TString inputlist = "", TString outfile = "testOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = 1)	//for simulation set approx 100 files and output name testOutFileSim.root
{
	//gStyle->SetOptStat(0); // JS. Mozesz wlaczyc jak nie chcesz boxa ze statystykami na wykresie
	gROOT->SetBatch(kTRUE);
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

		//simulation
		//TString inputFolder = "/lustre/nyx/hades/dstsim/mar19/ag1580ag/gen3/bmax14/no_enhancement_gcalor/root";
		
		//data
		TString inputFolder = "/lustre/hades/dst/apr12/gen9/122/root";
	
		TSystemDirectory* inputDir = new TSystemDirectory("inputDir", inputFolder);
		TList* files = inputDir->GetListOfFiles();

		for (Int_t i = 0; i <= files->LastIndex() && nFiles < maxFiles; i++) 
		{
			if (((TSystemFile*) files->At(i))->IsDirectory())
				continue;
			
			//if(TString(((TSystemFile*)files->At(i))->GetName()).Contains("accepted"))		//for data
			//{
				loop->addFile(inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName());
				nFiles++;
			//}
			
			//loop->addFile(inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName());		//for simulation
			//	nFiles++;
		}
	}
    
    //--------------------------------------------------------------------------------
    // Booking the categories to be read from the DST files.
    // By default all categories are booked therefore -* (Unbook all) first and book the ones needed
    // All required categories have to be booked except the global Event Header which is always booked
    //--------------------------------------------------------------------------------
    if (!loop->setInput("-*,+HGeantKine,+HParticleCand,+HParticleEvtInfo"))
		exit(1);

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
    HParticleCand*    particle_cand;	//dla symulacji jest TParticleCandSim*, bo inny obiekt (posiada inne informacje); dla danych jest TParticleCand*
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
    
    if (!particle_cand_cat) // If the category for the reconstructed trackes does not exist the macro makes no sense
		exit(1);

    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================
	
	


	//----------------------TCutG by Simon Spies -----------------------------
    
	TFile *file = new TFile("/lustre/nyx/hades/user/jszewczy/Mar19AgAg1580.root");
	TCutG *piPRPCCut = (TCutG*)file -> Get("tcgPBetaPiPRPC2Sig");
	TCutG *piNRPCCut = (TCutG*)file -> Get("tcgPBetaPiMRPC2Sig");
	TCutG *piPToFCut = (TCutG*)file -> Get("tcgPBetaPiPToF2Sig");
	TCutG *piNToFCut = (TCutG*)file -> Get("tcgPBetaPiMToF2Sig");
	
	//------------------------------------------------------------------------
    
	
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
		
		Float_t vertZ 			= EventVertex.Z();

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
	
		hCounter->Fill(cNumSelectedEvents);
	
		//================================================================================================================================================================
		// Put your analyses on event level here
		//================================================================================================================================================================
		
		int nTOFRPC = particle_info->getSumRpcMultHitCut() + particle_info->getSumTofMultCut();
		// remove comment to see the vortex
		//hXZ->Fill(EventVertex.Z(), EventVertex.X());
		//hYZ->Fill(EventVertex.Z(), EventVertex.Y());
		
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
	
			hCounter->Fill(cNumSelectedTracks);
	
			//--------------------------------------------------------------------------------
			// Getting information on the current track (Not all of them necessary for all analyses)
			//--------------------------------------------------------------------------------
			Short_t fCharge   = particle_cand->getCharge();
			Float_t fMomentum = particle_cand->getMomentum();
			Float_t fBeta	  = particle_cand->getBeta();
			Short_t fSystem   = particle_cand->getSystemUsed();
			//Short_t PID	    = particle_cand->getGeantPID();
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			


		} // End of track loop


	} // End of event loop
	
	

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
	
	
	
    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

