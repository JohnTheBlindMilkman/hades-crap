#include "Includes.h"
#include "FemtoMixer/EventCandidate.hxx"
#include "FemtoMixer/PairCandidate.hxx"
#include "FemtoMixer/JJFemtoMixer.hxx"
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
	TH1D hQinvSign,hQinvBckg;
	TH2D hDphiDthetaSign,hDphiDthetaBckg;
	TH3D hQoslSign,hQoslBckg;
};

std::size_t EventHashing(const Selection::EventCandidate &evt)
{
    return evt.GetCentrality()*1e1 + evt.GetPlate();
}

std::size_t PairHashing(const Selection::PairCandidate &pair)
{
	// collection of slices: (array[n],array[n+1]>
    // make sure this is ordered!
    constexpr std::array<float,6> ktArr{150,450,750,1050,1350,1650};
    constexpr std::array<float,4> yArr{0,0.49,0.99,1.49};
    constexpr std::array<float,9> EpArr{-202.5,-157.5,-112.5,-67.5,-22.5,22.5,67.5,112.5,157.5};

    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair.GetKt()) - ktArr.begin();
    std::size_t yCut = std::lower_bound(yArr.begin(),yArr.end(),pair.GetRapidity()) - yArr.begin();
    std::size_t EpCut = std::lower_bound(EpArr.begin(),EpArr.end(),pair.GetPhi()) - EpArr.begin();

	// reject if value is below first slice or above the last
	if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > yArr.size()-1 || EpCut == 0 || EpCut > EpArr.size()-1)
		return 0;
	else
    	return ktCut*1e2 + yCut*1e1 + EpCut;
}

bool PairRejection(const Selection::PairCandidate &pair)
{
	using Behaviour = Selection::PairCandidate::Behaviour;
	return pair.RejectPairByCloseHits<Behaviour::Uniform>(0,1) || pair.GetSharedMetaCells() > 0;
}

int newFemtoAnalysis(TString inputlist = "", TString outfile = "femtoOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)	//for simulation set approx 100 files and output name testOutFileSim.root
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
	TString beamtime="apr12";
	
	Int_t mdcMods[6][4]=
	{ {1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1} };
	TString asciiParFile     = "";
	TString rootParFile      = "/cvmfs/hadessoft.gsi.de/param/real/apr12/allParam_APR12_gen9_27092017.root";
	TString paramSource      = "root"; // root, ascii, oracle
	TString paramrelease     = "APR12_dst_gen9"; 
	HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,shower,wall,start,tbox");
	HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease); 

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
			//inputFolder = "/lustre/hades/user/kjedrzej/customDST/apr12PlusHMdcSeg/5sec/109"; // test Robert filtered events
			//inputFolder = "/lustre/hades/dst/feb24/gen0c/039/01/root"; // Au+Au 800 MeV
			inputFolder = "/lustre/hades/dst/apr12/gen9/122/root"; // Au+Au 2.4 GeV
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

	TH2D *hPhiTheta = new TH2D("hPhiTheta","#phi vs #theta distribution of tracks;#phi [deg];#theta [deg]",360,0,360,90,0,90);

	std::unordered_map<std::size_t, HistogramCollection> fMapFoHistograms;

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	// create objects for particle selection and mixing
	Selection::EventCandidate fEvent;	
	Selection::TrackCandidate fTrack;

	// create object for getting MDC wires
	HParticleWireInfo fWireInfo;

	std::unordered_map<std::size_t,std::vector<Selection::PairCandidate> > fSignMap, fBckgMap;

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixer;
	mixer.SetMaxBufferSize(50);
	mixer.SetEventHashingFunction(EventHashing);
	mixer.SetPairHashingFunction(PairHashing);
	mixer.SetPairCuttingFunction(PairRejection);
	
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
    masterTaskSet->add(matcher);

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
			cout << " Last events processed " << endl;
			break;
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
		
		Float_t vertZ = EventVertex.Z();
		Float_t vertX = EventVertex.X();
		Float_t vertY = EventVertex.Y();
		Int_t centClassIndex    = evtChara.getCentralityClass(eCentEst, eCentClass1); // 0 is overflow, 1 is 0-10, etc.
		Float_t EventPlane = evtChara.getEventPlane(eEPcorr);
		Float_t EventPlaneA = evtChara.getEventPlane(eEPcorr,1);
		Float_t EventPlaneB = evtChara.getEventPlane(eEPcorr,2);

		if (EventPlane < 0)
			continue;
		if (EventPlaneA < 0 || EventPlaneB < 0)
			continue;

		fEvent = Selection::EventCandidate(
			std::to_string(event_header->getEventRunNumber())+std::to_string(event_header->getEventSeqNumber()),
			vertX,
			vertY,
			vertZ,
			centClassIndex,
			EventPlane);

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
	
		//================================================================================================================================================================
		// Put your analyses on event level here
		//================================================================================================================================================================
		
		if (! fEvent.SelectEvent()) 
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
			HParticleWireManager &wire_manager = matcher->getWireManager();
			wire_manager.setWireRange(0);
			wire_manager.getWireInfo(track,fWireInfo,particle_cand);
			
			//--------------------------------------------------------------------------------
			// Discarding all tracks that have been discarded by the track sorter and counting all / good tracks
			//--------------------------------------------------------------------------------
			hCounter->Fill(cNumAllTracks);
	
			if (!particle_cand->isFlagBit(Particle::kIsUsed))
				continue;
	
			//--------------------------------------------------------------------------------
			// Getting information on the current track (Not all of them necessary for all analyses)
			//--------------------------------------------------------------------------------
			fTrack = Selection::TrackCandidate(
				particle_cand,Selection::TrackCandidate::CreateWireArray(fWireInfo),
				fEvent.GetID(),
				fEvent.GetReactionPlane(),
				track,
				14);
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fTrack.GetSystem() == Selection::Detector::RPC)
			{
				// fill RPC monitors for all tracks
			}
			else
			{
				// fill ToF monitors for all tracks
			}

			if (!fTrack.SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				continue;

			fEvent.AddTrack(fTrack);
			hPhiTheta->Fill(fTrack.GetPhi(),fTrack.GetTheta());
			hCounter->Fill(cNumSelectedTracks);

			if (fTrack.GetSystem() == Selection::Detector::RPC)
			{
				// fill RPC monitors for accepted tracks
			}
			else
			{
				// fill ToF monitors for accepted tracks
			}

		} // End of track loop

		if (fEvent.GetTrackListSize()) // if track vector has entries
		{
			hCounter->Fill(cNumAllPairs,fEvent.GetTrackListSize()*(fEvent.GetTrackListSize()-1)/2); // all combinations w/o repetitions

            fSignMap = mixer.AddEvent(fEvent,fEvent.GetTrackList());
			fBckgMap = mixer.GetSimilarPairs(fEvent);

			for (const auto &signalEntry : fSignMap)
			{
				for (const auto &entry : signalEntry.second)
				{
					if (fMapFoHistograms.find(signalEntry.first) == fMapFoHistograms.end())
					{
						HistogramCollection histos{
						TH1D(/* TString::Format("hQinvSign_%lu",signalEntry.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000 */),
						TH1D(/* TString::Format("hQinvBckg_%lu",signalEntry.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000 */),
						TH2D(/* TString::Format("hDphiDthetaSign_%lu",signalEntry.first), "#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45 */),
						TH2D(/* TString::Format("hDphiDthetaBckg_%lu",signalEntry.first), "#Delta#phi vs #Delta#theta distribution of backgound 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45 */),
						TH3D(TString::Format("hQoslSign_%lu",signalEntry.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",126,-500,500,126,-500,500,126,-500,500),
						TH3D(TString::Format("hQoslBckg_%lu",signalEntry.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",126,-500,500,126,-500,500,126,-500,500)
						};
						fMapFoHistograms.emplace(signalEntry.first,std::move(histos));
					}
					//fMapFoHistograms.at(signalEntry.first).hQinvSign.Fill(entry.GetQinv());
					//fMapFoHistograms.at(signalEntry.first).hDphiDthetaSign.Fill(entry.DeltaPhi,entry.DeltaTheta);
					float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(signalEntry.first).hQoslSign.Fill(qout,qside,qlong);

					if (signalEntry.first != 0)
						hCounter->Fill(cNumSelectedPairs);
				}
			}

			for (const auto &backgroundEntry : fBckgMap)
			{
				for (const auto &entry : backgroundEntry.second)
				{
					if (fMapFoHistograms.find(backgroundEntry.first) == fMapFoHistograms.end())
					{
						HistogramCollection histos{
						TH1D(/* TString::Format("hQinvSign_%lu",backgroundEntry.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000 */),
						TH1D(/* TString::Format("hQinvBckg_%lu",backgroundEntry.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000 */),
						TH2D(/* TString::Format("hDphiDthetaSign_%lu",backgroundEntry.first), "#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45 */),
						TH2D(/* TString::Format("hDphiDthetaBckg_%lu",backgroundEntry.first), "#Delta#phi vs #Delta#theta distribution of backgound 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",180,-45,45,180,-45,45 */),
						TH3D(TString::Format("hQoslSign_%lu",backgroundEntry.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",126,-500,500,126,-500,500,126,-500,500),
						TH3D(TString::Format("hQoslBckg_%lu",backgroundEntry.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",126,-500,500,126,-500,500,126,-500,500)
						};
						fMapFoHistograms.emplace(backgroundEntry.first,std::move(histos));
					}
					//fMapFoHistograms.at(backgroundEntry.first).hQinvBckg.Fill(entry.GetQinv());
					//fMapFoHistograms.at(backgroundEntry.first).hDphiDthetaBckg.Fill(entry.DeltaPhi,entry.DeltaTheta);
					float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(backgroundEntry.first).hQoslBckg.Fill(qout,qside,qlong);
				}
			}
			
		} // End of femto mixing
	} // End of event loop
	
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

	for (auto &histos : fMapFoHistograms)
	{
		//histos.second.hQinvSign.Write();
		//histos.second.hQinvBckg.Write();
		//histos.second.hDphiDthetaSign.Write();
		//histos.second.hDphiDthetaBckg.Write();
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

