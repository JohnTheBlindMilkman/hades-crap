#include "Includes.h"
#include "../../../../JJFemtoMixer/JJFemtoMixer.hxx"
#include "PairUtils.hxx"
#include "EventUtils.hxx"
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
	TH1D hQinvSign,hQinvBckg;
	TH3D hQoslSign,hQoslBckg;
};

int newFemtoAnalysis(TString inputlist = "", TString outfile = "femtoOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)
{
	gStyle->SetOptStat(0);
	gROOT->SetBatch(kTRUE);
	
	constexpr bool isCustomDst{false};
	constexpr bool isSimulation{false}; // for now this could be easly just const
	constexpr int protonPID{14};

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
	TString rootParFile;
	if (isSimulation)
	{
		rootParFile = "/cvmfs/hadessoft.gsi.de/param/sim/apr12/allParam_APR12_sim_run_12001_gen9_07112017.root";
	}
	else
	{
		//rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/apr12/allParam_APR12_gen9_27092017.root"; //gen9
		rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/apr12/allParam_APR12_gen10_16122024.root"; //gen10
		//rootParFile = "/cvmfs/hadessoft.gsi.de/param/real/feb24/allParam_feb24_gen0_16042024.root"; // Au+Au 800 MeV
	}
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
				//inputFolder = "/lustre/hades/dst/apr12/gen9/122/root"; // Au+Au 2.4 GeV gen9
				inputFolder = "/lustre/hades/dst/apr12/gen10/122/root"; // Au+Au 2.4 GeV gen10
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
    if (!loop->setInput(inputString.data()))
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
    HParticleCand*    particle_cand;	
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

	std::map<std::string, HistogramCollection> fMapFoHistograms;

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_3sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_3.0");
	TCutG* betamom_3sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_3.0");
	
	// create objects for particle selection and mixing
	std::shared_ptr<Selection::EventCandidate> fEvent;	
	std::shared_ptr<Selection::TrackCandidate> fTrack;

	// create object for getting MDC wires
	HParticleWireInfo fWireInfo;
	HGeantHeader *geantHeader;

	std::map<std::string,std::vector<std::shared_ptr<Selection::PairCandidate> > > fSignMap, fBckgMap;

    Mixing::JJFemtoMixer<Selection::EventCandidate,Selection::TrackCandidate,Selection::PairCandidate> mixer;
	mixer.SetMaxBufferSize((isSimulation) ? 200 : 50); // ana=50, sim=200
	mixer.SetEventHashingFunction(Mixing::EventGrouping{}.MakeEventGroupingFunction());
	mixer.SetPairHashingFunction(Mixing::PairGrouping{}.MakePairGroupingFunction1D());
	mixer.SetPairCuttingFunction(Mixing::PairRejection{}.MakePairRejectionFunction());
	
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
		Float_t EventPlane = -1;
		Float_t EventPlaneA = -1;
		Float_t EventPlaneB = -1;

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
		
		if (! fEvent->SelectEvent<HADES::Target::Setup::Apr12>({1},2,2,2))
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
			
			//fWireManager = matcher->getWireManager();
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
					nullptr,
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
				// fill RPC monitors for all tracks
			}
			else
			{
				// fill ToF monitors for all tracks
			}

			if (!fTrack->SelectTrack(betamom_3sig_p_rpc_pionCmom,betamom_3sig_p_tof_pionCmom))
				continue;

			//fSmearer.SmearMomenta(fTrack); // this will smear your momenta

			fEvent->AddTrack(fTrack);
			hPhiTheta->Fill(fTrack->GetPhi(),fTrack->GetTheta());
			hCounter->Fill(cNumSelectedTracks);

			if (fTrack->GetSystem() == Selection::Detector::RPC)
			{
				// fill RPC monitors for accepted tracks
			}
			else
			{
				// fill ToF monitors for accepted tracks
			}

		} // End of track loop

		if (fEvent->GetTrackListSize() > 2) // if track vector has entries
		{
            fSignMap = mixer.AddEvent(fEvent,fEvent->GetTrackList());
			fBckgMap = mixer.GetSimilarPairs(fEvent);

			for (const auto &signalEntry : fSignMap)
			{
				for (const auto &entry : signalEntry.second)
				{
					hCounter->Fill(cNumAllPairs);
					if (fMapFoHistograms.find(signalEntry.first) == fMapFoHistograms.end())
					{
						HistogramCollection histos{
						TH1D(TString::Format("hQinvSign_%s",signalEntry.first.data()),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckg_%s",signalEntry.first.data()),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",signalEntry.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",signalEntry.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fMapFoHistograms.emplace(signalEntry.first,std::move(histos));
					}
					fMapFoHistograms.at(signalEntry.first).hQinvSign.Fill(entry->GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(signalEntry.first).hQoslSign.Fill(qout,qside,qlong); */

					if (signalEntry.first != "bad" && signalEntry.first != "0")
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
						TH1D(TString::Format("hQinvSign_%s",backgroundEntry.first.data()),"Signal of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH1D(TString::Format("hQinvBckg_%s",backgroundEntry.first.data()),"Backgound of Protons 0-10%% centrality;q_{inv} [MeV/c];CF(q_{inv})",750,0,3000),
						TH3D(/* TString::Format("hQoslSign_%lu",backgroundEntry.first),"Signal of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */),
						TH3D(/* TString::Format("hQoslBckg_%lu",backgroundEntry.first),"Background of Protons 0-10%% centrality;q_{out} [MeV/c];q_{side} [MeV/c];q_{long} [MeV/c];CF(q_{inv})",64,0,500,64,0,500,64,0,500 */)
						};
						fMapFoHistograms.emplace(backgroundEntry.first,std::move(histos));
					}
					fMapFoHistograms.at(backgroundEntry.first).hQinvBckg.Fill(entry->GetQinv());
					/* float qout,qside,qlong;
					std::tie(qout,qside,qlong) = entry.GetOSL();
					fMapFoHistograms.at(backgroundEntry.first).hQoslBckg.Fill(qout,qside,qlong); */
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
    std::cout << "Finished DST processing" << endl;

	//--------------------------------------------------------------------------------
    // Showing how much of the buffer was used for each event hash
    //--------------------------------------------------------------------------------
	mixer.PrintStatus();

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
		//histos.second.hQoslSign.Write();
		//histos.second.hQoslBckg.Write();
	}

	hPhiTheta->Write();
	
    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    std::cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

