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

bool PairCut(const Selection::PairCandidate &pair)
{
	using Behaviour = Selection::PairCandidate::Behaviour;

	return pair.RejectPairByCloseHits<Behaviour::OneUnder>(0.7,2) ||
		pair.GetBothLayers() < 20 ||
		pair.GetSharedMetaCells() > 0; 
}

int newMomentumResolutionAnalysis(TString inputlist = "", TString outfile = "momResOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);

	constexpr int protonPID{14};
	constexpr std::size_t mixerBuffer{200};

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
	TString rootParFile = "/cvmfs/hadessoft.gsi.de/param/sim/apr12/allParam_APR12_sim_run_12001_gen9_07112017.root";

	TString paramSource = "root"; // root, ascii, oracle
	TString paramrelease = "APR12_dst_gen9";  // now, APR12_gen2_dst APR12_gen5_dst
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


		TString inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV
	
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
    if (!loop->setInput("-*,+HParticleCand,+HParticleEvtInfo,+HWallHit,+HGeantKine")) // make sure to use +HMdcSeg if you want the wires
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
	//HCategory *mdc_seg_cat = (HCategory*) gHades->getCurrentEvent()->getCategory(catMdcSeg); // so... this is my stupid fix for now, without it HFiredWires doesn't work
    
    if (!particle_cand_cat || !kine_cand_cat) // If the category for the reconstructed tracks does not exist, the macro makes no sense
		exit(1);
	
    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================

	std::map<std::size_t, HistogramCollection> fHistogramsGeantKine,fHistogramsPartCand;

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");

	// create objects for particle selection
	Selection::EventCandidate fEventGeantKine, fEventPartCand;
	Selection::TrackCandidate fTrackGeantKine, fTrackPartCand;
	HParticleWireInfo fWireInfo;
	std::size_t tracks;

	std::map<std::size_t,std::vector<Selection::PairCandidate> > fSignMapGeantKine, fBckgMapGeantKine, fSignMapPartCand, fBckgMapPartCand;

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixerGeantKine, mixerPartCand;
	mixerGeantKine.SetMaxBufferSize(mixerBuffer);
	mixerGeantKine.SetEventHashingFunction(EventHashing);
	mixerGeantKine.SetPairHashingFunction(PairHashing);
	mixerGeantKine.SetPairCuttingFunction(PairCut);
	mixerGeantKine.PrintSettings();

	mixerPartCand.SetMaxBufferSize(mixerBuffer);
	mixerPartCand.SetEventHashingFunction(EventHashing);
	mixerPartCand.SetPairHashingFunction(PairHashing);
	mixerPartCand.SetPairCuttingFunction(PairCut);
	mixerPartCand.PrintSettings();

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

	//--------------------------------------------------------------------------------
	// wire information w/o HMdcSeg class access
	//--------------------------------------------------------------------------------
	HTaskSet *masterTaskSet = gHades->getTaskSet("all");
    HParticleMetaMatcher* matcher = new HParticleMetaMatcher();
    matcher->setDebug();
	matcher->setUseEMC(kFALSE);
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
	TString ParameterfileCVMFS = "/cvmfs/hadessoft.gsi.de/param/eventchara/centrality_epcorr_sim_au1230au_gen9vertex_UrQMD_minbias_2019_04_pass0.root";
	
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

		fEventGeantKine = Selection::EventCandidate(
			std::to_string(event_header->getEventRunNumber())+std::to_string(event_header->getEventSeqNumber()),
			vertX,
			vertY,
			vertZ,
			centClassIndex,
			EventPlane);

		fEventPartCand = Selection::EventCandidate(
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
		
		if (!fEventGeantKine.SelectEvent({1,2,3,4,5,6},2,2,2) || !fEventPartCand.SelectEvent({1,2,3,4,5,6},2,2,2)) 
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
			//innerSeg = HCategoryManager::getObject(innerSeg,mdc_seg_cat,particle_cand->getInnerSegInd());
			//outerSeg = HCategoryManager::getObject(outerSeg,mdc_seg_cat,particle_cand->getOuterSegInd());
			//HParticleWireManager &wire_manager = matcher->getWireManager();
			//wire_manager.setDebug();
			//wire_manager.setWireRange(0);
			//wire_manager.getWireInfo(track,fWireInfo,particle_cand);
			matcher->getWireInfoDirect(particle_cand,fWireInfo);

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
			fTrackGeantKine = Selection::TrackCandidate(
				static_cast<HGeantKine*>(kine_cand_cat->getObject(particle_cand->getGeantTrack() - 1)),
				{},
				fEventGeantKine.GetID(),
				fEventGeantKine.GetReactionPlane(),
				track,
				14);
			fTrackPartCand = Selection::TrackCandidate(
				particle_cand,
				Selection::TrackCandidate::CreateWireArray(fWireInfo),
				fEventPartCand.GetID(),
				fEventPartCand.GetReactionPlane(),
				track,
				14);
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fTrackGeantKine.SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				fEventGeantKine.AddTrack(fTrackGeantKine);
			if (fTrackPartCand.SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				fEventPartCand.AddTrack(fTrackPartCand);

		} // End of track loop

		/* tracks = fEvent.GetTrackListSize();
		hCounter->Fill(cNumSelectedTracks,tracks); */

		if (fEventGeantKine.GetTrackListSize())
		{
			hCounter->Fill(cNumAllPairs,tracks*(tracks-1)/2); // all combinations w/o repetitions

			fSignMapGeantKine = mixerGeantKine.AddEvent(fEventGeantKine,fEventGeantKine.GetTrackList());
			fBckgMapGeantKine = mixerGeantKine.GetSimilarPairs(fEventGeantKine);

			for (const auto &pair : fSignMapGeantKine)
				for (const auto &elem : pair.second)
				{
					if (fHistogramsGeantKine.find(pair.first) == fHistogramsGeantKine.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSignGK_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckgGK_%lu",pair.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",pair.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fHistogramsGeantKine.emplace(pair.first,std::move(histos));
					}
					fHistogramsGeantKine.at(pair.first).hQinvSign.Fill(elem.GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(pair.first).hQoslSign.Fill(qout,qside,qlong); */

					if (pair.first != 0)
						hCounter->Fill(cNumSelectedPairs);
				}

			for (const auto &pair : fBckgMapGeantKine)
				for (const auto &elem : pair.second)
				{
					if (fHistogramsGeantKine.find(pair.first) == fHistogramsGeantKine.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSignGK_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckgGK_%lu",pair.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",pair.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fHistogramsGeantKine.emplace(pair.first,std::move(histos));
					}
					fHistogramsGeantKine.at(pair.first).hQinvBckg.Fill(elem.GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(pair.first).hQoslSign.Fill(qout,qside,qlong); */
				}
		}

		if (fEventPartCand.GetTrackListSize())
		{
			hCounter->Fill(cNumAllPairs,tracks*(tracks-1)/2); // all combinations w/o repetitions

			fSignMapPartCand = mixerPartCand.AddEvent(fEventPartCand,fEventPartCand.GetTrackList());
			fBckgMapPartCand = mixerPartCand.GetSimilarPairs(fEventPartCand);

			for (const auto &pair : fSignMapPartCand)
				for (const auto &elem : pair.second)
				{
					if (fHistogramsPartCand.find(pair.first) == fHistogramsPartCand.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSignPC_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckgPC_%lu",pair.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",pair.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fHistogramsPartCand.emplace(pair.first,std::move(histos));
					}
					fHistogramsPartCand.at(pair.first).hQinvSign.Fill(elem.GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(pair.first).hQoslSign.Fill(qout,qside,qlong); */

					if (pair.first != 0)
						hCounter->Fill(cNumSelectedPairs);
				}

			for (const auto &pair : fBckgMapPartCand)
				for (const auto &elem : pair.second)
				{
					if (fHistogramsPartCand.find(pair.first) == fHistogramsPartCand.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSignPC_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckgPC_%lu",pair.first),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",pair.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",pair.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fHistogramsPartCand.emplace(pair.first,std::move(histos));
					}
					fHistogramsPartCand.at(pair.first).hQinvBckg.Fill(elem.GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(pair.first).hQoslSign.Fill(qout,qside,qlong); */
				}
		}

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

	for(const auto elem : fHistogramsGeantKine)
	{
		elem.second.hQinvSign.Write();
		elem.second.hQinvBckg.Write();
	}
	for(const auto elem : fHistogramsPartCand)
	{
		elem.second.hQinvSign.Write();
		elem.second.hQinvBckg.Write();
	}
		
    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

