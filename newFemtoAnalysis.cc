#include "Includes.h"
#include "FemtoMixer/EventCandidate.hxx"
#include "FemtoMixer/PairCandidate.hxx"
#include "FemtoMixer/JJFemtoMixer.hxx"
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>

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
    // make sure this is ordered!
    constexpr std::array<float,7> ktArr{150,450,750,1050,1350,1650,1651};
    constexpr std::array<float,3> yArr{-0.5,0.5,0.6};
    constexpr std::array<float,6> EpArr{-180,-90,0,90,180,181};

    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair.Kt) - ktArr.begin();
    std::size_t yCut = std::lower_bound(yArr.begin(),yArr.end(),pair.Rapidity) - yArr.begin();
    std::size_t EpCut = std::lower_bound(EpArr.begin(),EpArr.end(),pair.AzimuthalAngle) - EpArr.begin();

	if (ktCut == 0 || ktCut == ktArr.size()-1 || yCut == 0 || yCut == yArr.size()-1 || EpCut == 0 || EpCut == EpArr.size()-1)
		return 0;
	else
    	return ktCut*1e2 + yCut*1e1 + EpCut;
}

int newFemtoAnalysis(TString inputlist = "", TString outfile = "femtoOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)	//for simulation set approx 100 files and output name testOutFileSim.root
{
	const int fTargetPlates = 15; // number of target plates

	const std::vector<std::pair<float,float> > fPairAzimuthBins = {{-180,-90},{-90,0},{0,90},{90,180}};
	const std::vector<std::pair<int,int> > fCentralityBins = {{1,2}};
	const std::vector<std::pair<float,float> > fKtBins = {{150,450},{450,750},{750,1050},{1050,1350},{1350,1650}};
	const std::vector<std::pair<float,float> > fPairRapidityBins = {{-0.5,0.5}};
	const int fEventsToMix = 50;

	gStyle->SetOptStat(0);
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
		//TString inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root";
		
		//data
		TString inputFolder = "/lustre/hades/dst/apr12/gen9/122/root";
	
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
    if (!loop->setInput("-*,+HGeantKine,+HParticleCand,+HParticleEvtInfo"))
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
    HParticleCand*    particle_cand;	//dla symulacji jest HParticleCandSim*, bo inny obiekt (posiada inne informacje); dla danych jest HParticleCand*
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
    
    if (!particle_cand_cat) // If the category for the reconstructed trackes does not exist the macro makes no sense
		exit(1);
	
    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================

	TH2D *hPhiTheta = new TH2D("hPhiTheta","#phi vs #theta distribution of tracks;#phi [deg];#theta [deg]",361,0,360,91,0,90);

	std::unordered_map<std::size_t, HistogramCollection> fMapFoHistograms;

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	// create objects for particle selection and mixing
	Selection::EventCandidate fEvent;	
	Selection::TrackCandidate fTrack;

	std::unordered_map<std::size_t,std::vector<Selection::PairCandidate> > fSignMap, fBckgMap;

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixer;
	mixer.SetMaxBufferSize(50);
	mixer.SetEventHashingFunction(EventHashing);
	mixer.SetPairHashingFunction(PairHashing);
	
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
	TString ParameterfileCVMFS = "/cvmfs/hadessoft.gsi.de/param/eventchara/centrality_epcorr_apr12_gen8_2019_02_pass30.root";

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
		Float_t EventPlane = evtChara.getEventPlane(eEPcorr);

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
	
		hCounter->Fill(cNumSelectedEvents);
	
		//================================================================================================================================================================
		// Put your analyses on event level here
		//================================================================================================================================================================
		
		if (!Selection::SelectEvent(fEvent)) 
			continue;
		
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

			fEvent.trackList.push_back(fTrack);
			hPhiTheta->Fill(fTrack.AzimimuthalAngle,fTrack.PolarAngle);

			if (fTrack.System == Selection::Detector::RPC)
			{
				// fill RPC monitors for accepted tracks
			}
			else
			{
				// fill ToF monitors for accepted tracks
			}

		} // End of track loop

		if (fEvent.trackList.size()) // if track vector has entries
		{
            fSignMap = mixer.AddEvent(fEvent,fEvent.trackList);
			fBckgMap = mixer.GetSimilarPairs(fEvent);

			for (const auto &signalEntry : fSignMap)
			{
				for (const auto &entry : signalEntry.second)
				{
					if (fMapFoHistograms.find(signalEntry.first) == fMapFoHistograms.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSign_%lu",signalEntry.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckg_%lu",signalEntry.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000),
						TH2D(TString::Format("hDphiDthetaSign_%lu",signalEntry.first), "#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90),
						TH2D(TString::Format("hDphiDthetaBckg_%lu",signalEntry.first), "#Delta#phi vs #Delta#theta distribution of backgound 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90),
						TH3D(),
						TH3D()};
						fMapFoHistograms.emplace(signalEntry.first,histos);
					}
					fMapFoHistograms.at(signalEntry.first).hQinvSign.Fill(entry.QInv);
					fMapFoHistograms.at(signalEntry.first).hDphiDthetaSign.Fill(entry.DeltaPhi,entry.DeltaTheta);
				}
			}

			for (const auto &backgroundEntry : fBckgMap)
			{
				for (const auto &entry : backgroundEntry.second)
				{
					if (fMapFoHistograms.find(backgroundEntry.first) == fMapFoHistograms.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSign_%lu",backgroundEntry.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckg_%lu",backgroundEntry.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000),
						TH2D(TString::Format("hDphiDthetaSign_%lu",backgroundEntry.first), "#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90),
						TH2D(TString::Format("hDphiDthetaBckg_%lu",backgroundEntry.first), "#Delta#phi vs #Delta#theta distribution of backgound 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90),
						TH3D(),
						TH3D()};
						fMapFoHistograms.emplace(backgroundEntry.first,histos);
					}
					fMapFoHistograms.at(backgroundEntry.first).hQinvBckg.Fill(entry.QInv);
					fMapFoHistograms.at(backgroundEntry.first).hDphiDthetaBckg.Fill(entry.DeltaPhi,entry.DeltaTheta);
				}
			}
			
		} // End of femto mixing
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

	for (auto &histos : fMapFoHistograms)
	{
		histos.second.hQinvSign.Write();
		histos.second.hQinvBckg.Write();
		histos.second.hDphiDthetaSign.Write();
		histos.second.hDphiDthetaBckg.Write();
		histos.second.hQoslSign.Write();
		histos.second.hQoslBckg.Write();
	}

	hPhiTheta->Write();
	
    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

