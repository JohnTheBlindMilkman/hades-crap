#include "Includes.h"
#include "FemtoMixer/EventCandidate.hxx"
#include "FemtoMixer/PairCandidate.hxx"
#include "../JJFemtoMixer/JJFemtoMixer.hxx"
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>
#include <type_traits>

// defining a helper function
template <typename T>
bool isSim(T *t) {return false;}
template <>
bool isSim(HParticleCandSim *t) {return true;}

struct HistogramCollection
{
	TH1D hQinvSign;
	TH3D hQoslSign;
};

std::string EventHashing(const std::shared_ptr<Selection::EventCandidate> &evt)
{
    return std::to_string(static_cast<std::size_t>(evt->GetCentrality())) + JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
	// return JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetNCharged()/10),2) + 
	// 	JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetReactionPlane()/10),2) + 
	// 	JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
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

int newPurityAnalysis(TString inputlist = "", TString outfile = "purityOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = 10)	//for simulation set approx 100 files and output name testOutFileSim.root
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);
	
	constexpr int protonPID{14};
	constexpr std::size_t mixerBuffer{0};

	//--------------------------------------------------------------------------------
    // Initialization of the global ROOT object and the Hades Loop
    // The Hades Loop used as an interface to the DST data (Basically a container of a TChain).
    // kTRUE - The global HADES object is being created if not existing
    //--------------------------------------------------------------------------------
    TROOT dst_analysis("DstAnalysisMacro", "Simple DST analysis Macro");
    HLoop* loop = new HLoop(kTRUE);
	const TString beamtime="apr12";
	
	Int_t mdcMods[6][4]=
	{ {1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1},
	{1,1,1,1} };
	TString asciiParFile     = "";
	TString rootParFile = "/cvmfs/hadessoft.gsi.de/param/sim/apr12/allParam_APR12_sim_run_12001_gen9_07112017.root";

	TString paramSource      = "root"; // root, ascii, oracle
	TString paramrelease     = "APR12_dst_gen10"; 
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

		//inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 800 MeV
		//inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV gen9
		inputFolder = "/lustre/hades/dstsim/apr12/au1230au/gen10/bmax10/no_enhancement_gcalor/root"; // Au+Au 2.4 GeV gen10
	
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
    HParticleCandSim*    particle_cand;	//dla symulacji jest HParticleCandSim*, bo inny obiekt (posiada inne informacje); dla danych jest HParticleCand*
    HEventHeader*     event_header;
    HParticleEvtInfo* particle_info;

    HCategory* particle_info_cat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particle_cand_cat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
	HCategory* kine_cand_cat = (HCategory*) HCategoryManager::getCategory(catGeantKine);

    if (!particle_cand_cat || !kine_cand_cat) // If the category for the reconstructed trackes does not exist the macro makes no sense
		exit(1);
	
    //================================================================================================================================================================
    // Put your object declarations here
    //================================================================================================================================================================

	std::map<std::string, HistogramCollection> fMapFoHistogramsNum,fMapFoHistogramsDen;

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	// create objects for particle selection and mixing
	std::shared_ptr<Selection::EventCandidate> fEventNum,fEventDen;
	std::shared_ptr<Selection::TrackCandidate> fTrackNum,fTrackDen;
	HGeantHeader *geantHeader;

	// create object for getting MDC wires
	HParticleWireInfo fWireInfo;

	std::map<std::string,std::vector<std::shared_ptr<Selection::PairCandidate> > > fSignMapNum, fSignMapDen;

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixerNum,mixerDen;
	mixerNum.SetMaxBufferSize(mixerBuffer);
	mixerNum.SetEventHashingFunction(EventHashing);
	mixerNum.SetPairHashingFunction(PairHashing);
	mixerNum.SetPairCuttingFunction(PairCut);
	mixerNum.PrintSettings();

	mixerDen.SetMaxBufferSize(mixerBuffer);
	mixerDen.SetEventHashingFunction(EventHashing);
	mixerDen.SetPairHashingFunction(PairHashing);
	mixerDen.SetPairCuttingFunction(PairCut);
	mixerDen.PrintSettings();
	
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
			std::cout << " Last events processed " << endl;
			break;
		}

		TString tmp; // dummy variable; required by HLoop::isNewFile
		if (loop->isNewFile(tmp))
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
		
		Float_t vertZ = EventVertex.Z();
		Float_t vertX = EventVertex.X();
		Float_t vertY = EventVertex.Y();
		Int_t centClassIndex    = evtChara.getCentralityClass(eCentEst, eCentClass1); // 0 is overflow, 1 is 0-10, etc.

		geantHeader = loop->getGeantHeader();
		if (geantHeader == nullptr)
			continue;
		float EventPlane = geantHeader->getEventPlane() * TMath::DegToRad();
		float EventPlaneA = EventPlane;
		float EventPlaneB = EventPlane;

		if (EventPlane < 0)
			continue;
		if (EventPlaneA < 0 || EventPlaneB < 0)
			continue;
		
		fEventNum = std::make_shared<Selection::EventCandidate>(event_header,particle_info,centClassIndex,EventPlane);
		fEventDen = std::make_shared<Selection::EventCandidate>(event_header,particle_info,centClassIndex,EventPlane);

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
		
		if (! fEventNum->SelectEvent<HADES::Target::Setup::Apr12>({1},2,2,2) || ! fEventDen->SelectEvent<HADES::Target::Setup::Apr12>({1},2,2,2))
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
			//fWireManager = matcher->getWireManager();
			matcher->getWireInfoDirect(particle_cand,fWireInfo);
			//fWireManager.setWireRange(0);
			//fWireManager.getWireInfo(track,fWireInfo,particle_cand);
			
			// I have no freakin idea if this is how it should be done
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
			fTrackNum = std::make_shared<Selection::TrackCandidate>(
					particle_cand,
					static_cast<HGeantKine*>(kine_cand_cat->getObject(particle_cand->getGeantTrack() - 1)),
					HADES::MDC::CreateTrackLayers(fWireInfo),
					fEventNum->GetID(),
					fEventNum->GetReactionPlane(),
					track,
					protonPID);

			fTrackDen = std::make_shared<Selection::TrackCandidate>(
					particle_cand,
					static_cast<HGeantKine*>(kine_cand_cat->getObject(particle_cand->getGeantTrack() - 1)),
					HADES::MDC::CreateTrackLayers(fWireInfo),
					fEventDen->GetID(),
					fEventDen->GetReactionPlane(),
					track,
					protonPID);
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================


			if (fTrackNum->SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				fEventNum->AddTrack(fTrackNum);

			if (fTrackDen->SelectTrack(betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom,false))
				fEventDen->AddTrack(fTrackDen);

			hCounter->Fill(cNumSelectedTracks);

		} // End of track loop

		if (fEventNum->GetTrackListSize() > 2) // if track vector has entries
		{
            fSignMapNum = mixerNum.AddEvent(fEventNum,fEventNum->GetTrackList());

			for (const auto &signalEntry : fSignMapNum)
			{
				for (const auto &entry : signalEntry.second)
				{
					if (fMapFoHistogramsNum.find(signalEntry.first) == fMapFoHistogramsNum.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvNum_%s",signalEntry.first.data()),"Numerator of Proton Purity 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslNum_%lu",signalEntry.first),"Purity of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						};
						fMapFoHistogramsNum.emplace(signalEntry.first,std::move(histos));
					}
					fMapFoHistogramsNum.at(signalEntry.first).hQinvSign.Fill(entry->GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistogramsNum.at(signalEntry.first).hQoslSign.Fill(qout,qside,qlong); */
				}
			}
		}

		if (fEventDen->GetTrackListSize() > 2) // if track vector has entries
		{
            fSignMapDen = mixerDen.AddEvent(fEventDen,fEventDen->GetTrackList());

			for (const auto &signalEntry : fSignMapDen)
			{
				for (const auto &entry : signalEntry.second)
				{
					if (fMapFoHistogramsDen.find(signalEntry.first) == fMapFoHistogramsDen.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvDen_%s",signalEntry.first.data()),"Denominator of Proton Purity 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslDen_%lu",signalEntry.first),"Purity of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						};
						fMapFoHistogramsDen.emplace(signalEntry.first,std::move(histos));
					}
					fMapFoHistogramsDen.at(signalEntry.first).hQinvSign.Fill(entry->GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistogramsDen.at(signalEntry.first).hQoslSign.Fill(qout,qside,qlong); */
				}
			}
		}

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
    std::cout << "Finished DST processing" << endl;

	//--------------------------------------------------------------------------------
    // Showing how much of the buffer was used for each event hash
    //--------------------------------------------------------------------------------
	mixerNum.PrintStatus();
	mixerDen.PrintStatus();

    //--------------------------------------------------------------------------------
    // Creating output file and storing results there
    //--------------------------------------------------------------------------------
    TFile* out = new TFile(outfile.Data(), "RECREATE");
    out->cd();

    hCounter->Write();
	
    //================================================================================================================================================================
    // Remember to write your results to the output file here
    //================================================================================================================================================================

	for (const auto &histos : fMapFoHistogramsNum)
	{
		histos.second.hQinvSign.Write();
		//histos.second.hQoslSign.Write();
	}

	for (const auto &histos : fMapFoHistogramsDen)
	{
		histos.second.hQinvSign.Write();
		//histos.second.hQoslSign.Write();
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

