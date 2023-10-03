#include "Includes.h"
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>

enum class Detector {RPC, ToF};

struct Proton
{
	TLorentzVector momVec;
	int centrality;
	float zVertex;
	Detector sys;
};

struct ProtonPair
{
	Proton particle1;
	Proton particle2;
	float qInv;
	float qOut;
	float qSide;
	float qLong;
	float kT;
};

ProtonPair CFKinematics(Proton part1, Proton part2)
{
	ProtonPair pair;
	pair.particle1 = part1;
	pair.particle2 = part2;

	// adapted from https://github.com/DanielWielanek/HAL/blob/main/analysis/femto/base/FemtoPairKinematics.cxx
	Double_t tPx = part1.momVec.Px() + part2.momVec.Px();
    Double_t tPy = part1.momVec.Py() + part2.momVec.Py();
    Double_t tPz = part1.momVec.Pz() + part2.momVec.Pz();
    Double_t tE = part1.momVec.E() + part2.momVec.E();
    Double_t tPt = tPx * tPx + tPy * tPy;
    Double_t tMt = tE * tE - tPz * tPz;  // mCVK;
    tMt = TMath::Sqrt(tMt);
    pair.kT = TMath::Sqrt(tPt);
    Double_t tBeta  = tPz / tE;
    Double_t tGamma = tE / tMt;

    // Transform to LCMS

    Double_t particle1lcms_pz = tGamma * (part1.momVec.Pz() - tBeta * part1.momVec.E());
    Double_t particle1lcms_e  = tGamma * (part1.momVec.E() - tBeta * part1.momVec.Pz());
    Double_t particle2lcms_pz = tGamma * (part2.momVec.Pz() - tBeta * part2.momVec.E());
    Double_t particle2lcms_e  = tGamma * (part2.momVec.E() - tBeta * part2.momVec.Pz());

    // Rotate in transverse plane

    Double_t particle1lcms_px = (part1.momVec.Px() * tPx + part1.momVec.Py() * tPy) / pair.kT;
    Double_t particle1lcms_py = (-part1.momVec.Px() * tPy + part1.momVec.Py() * tPx) / pair.kT;

    Double_t particle2lcms_px = (part2.momVec.Px() * tPx + part2.momVec.Py() * tPy) / pair.kT;
    Double_t particle2lcms_py = (-part2.momVec.Px() * tPy + part2.momVec.Py() * tPx) / pair.kT;

    pair.qOut = abs(particle1lcms_px - particle2lcms_px);
    pair.qSide = abs(particle1lcms_py - particle2lcms_py);
    pair.qLong = abs(particle1lcms_pz - particle2lcms_pz);
    Double_t mDE = particle1lcms_e - particle2lcms_e;
    pair.qInv = TMath::Sqrt(TMath::Abs(pair.qOut * pair.qOut + pair.qSide * pair.qSide + pair.qLong * pair.qLong - mDE * mDE));

	return pair;
}

void CalcSignal(ROOT::VecOps::RVec<ProtonPair> &pairVec, const ROOT::VecOps::RVec<Proton> &trackVec)
{ 
	for (size_t iter1 = 0; iter1 < trackVec.size(); iter1++)
		for (size_t iter2 = iter1+1; iter2 < trackVec.size(); iter2++)
			pairVec.push_back(CFKinematics(trackVec.at(iter1),trackVec.at(iter2)));
}

void CalcBackground(ROOT::VecOps::RVec<ProtonPair> &pairVec, const ROOT::VecOps::RVec<Proton> &trackVec, const std::deque<Proton> &bckgVec)
{
	for (auto &track : trackVec)
		for (auto &bckg : bckgVec)
				pairVec.push_back(CFKinematics(track,bckg));
}

bool MixBckg(const ROOT::VecOps::RVec<Proton> &trackVec, std::deque<Proton> &bckgVec, const int &maxSize, std::mt19937 &generator)
{
	if (trackVec.size())
	{
		std::uniform_int_distribution<int> dist(0,trackVec.size()-1);
		bckgVec.push_back(trackVec.at(dist(generator)));
		if (bckgVec.size() > maxSize)
			bckgVec.pop_front();

		return true;
	}
	else
		return false;
}

void FillAndClear(ROOT::VecOps::RVec<ProtonPair> &pairVec, TH1D *hInp1, TH3D *hInp3)
{

	for (auto &pair : pairVec)
	{
		hInp1->Fill(pair.qInv);
		hInp3->Fill(pair.qOut,pair.qSide,pair.qLong);
	}

	pairVec.clear();
	pairVec.resize(0);
}

int testAnalysis(TString inputlist = "", TString outfile = "testOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = 10)	//for simulation set approx 100 files and output name testOutFileSim.root
{
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
		//TString inputFolder = "/lustre/nyx/hades/dstsim/mar19/ag1580ag/gen3/bmax14/no_enhancement_gcalor/root";
		
		//data
		TString inputFolder = "/u/kjedrzej/lustre/hades/dst/apr12/gen9/122/root";
	
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
	
	int nBinsX = 600, nMinX = -3000, nMaxX = 3000, nBinsY = 200;
	float nMinY = 0, nMaxY = 1.2;

	TH1D* hZvertex = new TH1D("hZvertex","",700,-65,5);
	TH1D *hAccTracks = new TH1D("hAccTracks","Number of accepted tracks",100,0,100);
	TH1D *hM2AccTof = new TH1D("hM2AccTof","m^{2} distribution fo accepted protons",1000,5e5,1.5e6);
	TH1D *hM2AccRpc = new TH1D("hM2AccRpc","m^{2} distribution fo accepted protons",1000,5e5,1.5e6);

    TH2D *hBetaMomTofAll = new TH2D("hBetaMomTofAll",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomTofProt = new TH2D("hBetaMomTofProt",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomRpcAll = new TH2D("hBetaMomRpcAll",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomRpcProt = new TH2D("hBetaMomRpcProt",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);

	TH1D *hQinvSign = new TH1D("hQinvSign","Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",250,0,1000);
	TH3D *hQoslSign = new TH3D("hQoslSign","Signal of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);
	TH1D *hQinvBckg = new TH1D("hQinvBckg","Background of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",250,0,1000);
	TH3D *hQoslBckg = new TH3D("hQoslBckg","Background of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);

	TFile *cutfile_betamom_pionCmom = new TFile("/u/kjedrzej/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	//create RNG
	std::random_device rd;
	std::mt19937 gen(rd());
	
	// create objects for mixing
	Proton fProton;
	Detector fDetector;
	bool fIsAccepted = false;
	ROOT::VecOps::RVec<Proton> fTrackVec;
	std::deque<Proton> fBckgCocktail;
	ROOT::VecOps::RVec<ProtonPair> fSignVec, fBckgVec;
	const int partToMix = 50;
	
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
		
		Float_t vertZ 			= EventVertex.Z();
		Int_t centClassIndex    = evtChara.getCentralityClass(eCentEst, eCentClass1); // 0 is overflow, 1 is 0-10, etc.

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
		
		if (centClassIndex != 1) // if centrality is not 0-10% 
			continue;

		hZvertex->Fill(vertZ);

		// remove comment to see the vertex
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
			Float_t fMass2 = fMomentum*fMomentum*(1-fBeta*fBeta)/(fBeta*fBeta);
			Short_t fSystem   = particle_cand->getSystemUsed();
			Bool_t fIsAtMdcEdge = particle_cand->isAtAnyMdcEdge();
			//Short_t PID	    = particle_cand->getGeantPID();
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fIsAtMdcEdge || fSystem == -1)
				continue;

			float fMomCharg = fCharge*fMomentum;

			switch (fSystem)
			{
				case 0:
					hBetaMomRpcAll->Fill(fMomCharg,fBeta);
					if (betamom_2sig_p_rpc_pionCmom->IsInside(fMomCharg,fBeta))
					{
						hM2AccRpc->Fill(fMass2);
						hBetaMomRpcProt->Fill(fMomCharg,fBeta);
						fIsAccepted = true;
						fDetector = Detector::RPC;
					}
					else
					{
						fIsAccepted = false;
					}
					break;

				case 1:
					hBetaMomTofAll->Fill(fMomCharg,fBeta);
					if (betamom_2sig_p_tof_pionCmom->IsInside(fMomCharg,fBeta))
					{
						hM2AccTof->Fill(fMass2);
						hBetaMomTofProt->Fill(fMomCharg,fBeta);
						fIsAccepted = true;
						fDetector = Detector::ToF;
					}
					else
					{
						fIsAccepted = false;
					}
					break;
				
				default:
					fIsAccepted = false;
					break;
			}

			if(fIsAccepted)
			{
				particle_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
				fProton.momVec = *particle_cand;
				fProton.centrality = centClassIndex;
				fProton.zVertex = 0;
				fProton.sys = fDetector;
				fTrackVec.push_back(fProton);
			}

		} // End of track loop

		if (MixBckg(fTrackVec,fBckgCocktail,partToMix,gen)) // if track vector has entries
		{
			hAccTracks->Fill(fTrackVec.size());

			CalcSignal(fSignVec,fTrackVec);
			CalcBackground(fBckgVec,fTrackVec,fBckgCocktail);
			FillAndClear(fSignVec,hQinvSign,hQoslSign);
			FillAndClear(fBckgVec,hQinvBckg,hQoslBckg);

			fTrackVec.clear();
			fTrackVec.resize(0);
		}
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
	
	hZvertex->Write();
	hAccTracks->Write();
	hM2AccRpc->Write();
	hM2AccTof->Write();

	hBetaMomRpcAll->Write();
	hBetaMomRpcProt->Write();
	hBetaMomTofAll->Write();
	hBetaMomTofProt->Write();

	hQinvSign->Write();
	hQoslSign->Write();
	hQinvBckg->Write();
	hQoslBckg->Write();
	
    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------
    out->Save();
    out->Close();

    cout << "####################################################" << endl;
	gROOT->SetBatch(kFALSE);
	return 0;
	}

