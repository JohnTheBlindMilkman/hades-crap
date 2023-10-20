#include "Includes.h"
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>

enum class Detector {RPC, ToF};

struct Observables
{
	float charge, momentum, beta, mass2, rapidity, pt;
	short int system, PID;
	bool isAtMdcEdge;
};

struct Proton
{
	TLorentzVector momVec;
	int centrality;
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

std::tuple<bool,int> SelectEvent(float xVertex, float yVertex, float zVertex, int nSigma = 1)
{
	// copied from zVertexPeakAndSigma.txt file
	static std::vector<std::pair<float,float>> plateVector = {{-54.7598,0.755846},{-51.6971,0.783591},{-47.7996,0.763387},{-44.5473,0.769386},{-40.569,0.781312},{-37.2151,0.762538},{-33.2948,0.76901},{-30.3726,0.742618},{-26.648,0.748409},{-22.5492,0.738462},{-18.9649,0.747727},{-15.5259,0.749724},{-11.8726,0.740386},{-8.45083,0.742672},{-4.58076,0.712394}};
	static std::pair<float,float> xPosition = {0.1951,0.619};
	static std::pair<float,float> yPosition = {0.7045,0.6187};

	if ((xVertex < (xPosition.first - nSigma*xPosition.second)) || (xVertex > (xPosition.first + nSigma*xPosition.second)))
		return std::make_tuple<bool,int>(false,-999);
	if ((yVertex < (yPosition.first - nSigma*yPosition.second)) || (yVertex > (yPosition.first + nSigma*yPosition.second)))
		return std::make_tuple<bool,int>(false,-999);

	int plateNo = 0;
	for(const auto &plate : plateVector)
	{
		if ((zVertex < (plate.first - nSigma*plate.second)) || (zVertex > (plate.first + nSigma*plate.second)))
			return std::make_tuple<bool,int>(true,std::move(plateNo));
		plateNo++;
	}

	return std::make_tuple<bool,int>(false,-999);
}

std::tuple<bool,Detector> SelectProton(Observables obs, TCutG *rpcCut, TCutG *tofCut)
{
	bool isAccepted = false;
	Detector det;
	switch (obs.system)
	{
		case 0:
			det = Detector::RPC;
			if (rpcCut->IsInside(obs.momentum*obs.charge,obs.beta))
				isAccepted = true;

			break;

		case 1:
			det = Detector::ToF;
			if (tofCut->IsInside(obs.momentum*obs.charge,obs.beta))
				isAccepted = true;

			break;
	}
	
	return std::make_tuple<bool,Detector>(std::move(isAccepted),std::move(det));
}

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

void CalcSignal(std::vector<ProtonPair> &pairVec, const std::vector<Proton> &trackVec, float openAngle = 0.0)
{ 
	for (size_t iter1 = 0; iter1 < trackVec.size(); iter1++)
		for (size_t iter2 = iter1+1; iter2 < trackVec.size(); iter2++)
			if (trackVec.at(iter1).momVec.Angle(trackVec.at(iter2).momVec.Vect()) > openAngle)
				pairVec.push_back(CFKinematics(trackVec.at(iter1),trackVec.at(iter2)));
}

void CalcBackground(std::vector<ProtonPair> &pairVec, const std::vector<Proton> &trackVec, const std::deque<Proton> &bckgVec, float openAngle = 0.0)
{
	for (auto &track : trackVec)
		for (auto &bckg : bckgVec)
			if (track.momVec.Angle(bckg.momVec.Vect()) > openAngle)
				pairVec.push_back(CFKinematics(track,bckg));
}

void MixBckg(const std::vector<Proton> &trackVec, std::deque<Proton> &bckgVec, const int &maxSize, std::mt19937 &generator)
{
	std::uniform_int_distribution<int> dist(0,trackVec.size()-1);
	bckgVec.push_back(trackVec.at(dist(generator)));
	if (bckgVec.size() > maxSize)
		bckgVec.pop_front();
}

int AssignKt(float kt)
{
	static std::vector<std::pair<float,float> > ktBinVec = {{150,450},{450,700},{700,1050},{1050,1350},{1350,1650}};
	for (std::size_t i = 0; i < ktBinVec.size(); i++)
		if (kt >= ktBinVec.at(i).first && kt <= ktBinVec.at(i).second)
			return i;

	return -999;
}

void FillAndClear(std::vector<ProtonPair> &pairVec, std::array<TH1D*,5> hInp1, std::array<TH3D*,5> hInp3)
{
	int ktBin = 0;
	for (auto &pair : pairVec)
	{
		ktBin = AssignKt(pair.kT);
		if (ktBin == -999)
			continue;

		hInp1.at(ktBin)->Fill(pair.qInv);
		hInp3.at(ktBin)->Fill(pair.qOut,pair.qSide,pair.qLong);
	}

	pairVec.clear();
	pairVec.resize(0);
}

int testAnalysis(TString inputlist = "", TString outfile = "testOutFile2.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)	//for simulation set approx 100 files and output name testOutFileSim.root
{
	const int fTargetPlates = 15; // number of target plates
	const int fKtBins = 5; // number of kT bins
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

	TH1D* hXvertex = new TH1D("hXvertex","",401,-20,20);
	TH1D* hYvertex = new TH1D("hYvertex","",401,-20,20);
	TH1D* hZvertex = new TH1D("hZvertex","",700,-65,5);
	TH1D *hAccTracks = new TH1D("hAccTracks","Number of accepted tracks",100,0,100);
	TH1D *hM2AccTof = new TH1D("hM2AccTof","m^{2} distribution fo accepted protons",1000,5e5,1.5e6);
	TH1D *hM2AccRpc = new TH1D("hM2AccRpc","m^{2} distribution fo accepted protons",1000,5e5,1.5e6);

    TH2D *hBetaMomTofAll = new TH2D("hBetaMomTofAll",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomTofProt = new TH2D("hBetaMomTofProt",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomRpcAll = new TH2D("hBetaMomRpcAll",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);
	TH2D *hBetaMomRpcProt = new TH2D("hBetaMomRpcProt",";p #times c;#beta",nBinsX,-nMaxX,nMaxX,nBinsY,nMinY,nMaxY);

	TH2D *hPtYTofProt = new TH2D("hPtYTofProt",";y_{c.m.};p_{T} [MeV/c]",121,-0.4,2,2000,0,2000);
	TH2D *hPtYRpcProt = new TH2D("hPtYRpcProt",";y_{c.m.};p_{T} [MeV/c]",121,-0.4,2,2000,0,2000);

	TH1D *hKtSign = new TH1D("hKtSign","Distribution of k_{T} in signal;k_{T} [MeV]",1000,0,2500);
	TH1D *hKtBckg = new TH1D("hKtBckg","Distribution of k_{T} in background;k_{T} [MeV]",1000,0,2500);
	std::array<TH1D*,fKtBins> hQinvSign,hQinvBckg;
	std::array<TH3D*,fKtBins> hQoslSign,hQoslBckg;
	for (int i = 0; i < fKtBins; i++)
	{
		hQinvSign.at(i) = new TH1D(TString::Format("hQinvSign%d",i),"Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",250,0,1000);
		hQoslSign.at(i) = new TH3D(TString::Format("hQoslSign%d",i),"Signal of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);
		hQinvBckg.at(i) = new TH1D(TString::Format("hQinvBckg%d",i),"Background of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",250,0,1000);
		hQoslBckg.at(i) = new TH3D(TString::Format("hQoslBckg%d",i),"Background of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);
	}

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	//create RNG
	std::random_device rd;
	std::mt19937 gen(rd());

	// create temporary obj for basic data storage
	std::tuple<bool,Detector> fIsTrackAccepted;
	std::tuple<bool,int> fIsEventAccepted;
	Observables fHADES;
	
	// create objects for mixing
	Proton fProton;
	std::vector<Proton> fTrackVec;
	std::array<std::deque<Proton>,fTargetPlates> fBckgCocktail;
	std::vector<ProtonPair> fSignVec, fBckgVec;
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
		
		Float_t vertZ = EventVertex.Z();
		Float_t vertX = EventVertex.X();
		Float_t vertY = EventVertex.Y();
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

		fIsEventAccepted = SelectEvent(vertX,vertY,vertZ);
		if (!std::get<0>(fIsEventAccepted))
			continue;
		//hZvertex->Fill(vertZ);
		//hXvertex->Fill(vertX);
		//hYvertex->Fill(vertY);

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
			fHADES.charge   = particle_cand->getCharge();
			fHADES.momentum = particle_cand->getMomentum();
			fHADES.beta	  = particle_cand->getBeta();
			fHADES.mass2 = fHADES.momentum*fHADES.momentum*(1-fHADES.beta*fHADES.beta)/(fHADES.beta*fHADES.beta);
			fHADES.system   = particle_cand->getSystemUsed();
			fHADES.isAtMdcEdge = particle_cand->isAtAnyMdcEdge();
			//fPID	    = particle_cand->getGeantPID();
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fHADES.isAtMdcEdge || fHADES.system == -1)
				continue;

			fIsTrackAccepted = SelectProton(fHADES,betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom);

			if (std::get<1>(fIsTrackAccepted) == Detector::RPC)
				hBetaMomRpcAll->Fill(fHADES.momentum*fHADES.charge,fHADES.beta);
			else
				hBetaMomTofAll->Fill(fHADES.momentum*fHADES.charge,fHADES.beta);

			if(std::get<0>(fIsTrackAccepted))
			{
				particle_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
				fProton.momVec = *particle_cand;
				fProton.centrality = centClassIndex;
				fProton.sys = std::get<1>(fIsTrackAccepted);
				fTrackVec.push_back(fProton);

				fHADES.rapidity = fProton.momVec.Rapidity();
				fHADES.pt = fProton.momVec.Perp();

				if (std::get<1>(fIsTrackAccepted) == Detector::RPC)
				{
					hM2AccRpc->Fill(fHADES.mass2);
					hBetaMomRpcProt->Fill(fHADES.momentum*fHADES.charge,fHADES.beta);
					hPtYRpcProt->Fill(fHADES.rapidity,fHADES.pt);
				}
				else
				{
					hM2AccTof->Fill(fHADES.mass2);
					hBetaMomTofProt->Fill(fHADES.momentum*fHADES.charge,fHADES.beta);
					hPtYTofProt->Fill(fHADES.rapidity,fHADES.pt);
				}
			}

		} // End of track loop

		if (fTrackVec.size()) // if track vector has entries
		{
			hAccTracks->Fill(fTrackVec.size());

			if (fBckgCocktail.at(std::get<1>(fIsEventAccepted)).size())
			{
				CalcSignal(fSignVec,fTrackVec);
				CalcBackground(fBckgVec,fTrackVec,fBckgCocktail.at(std::get<1>(fIsEventAccepted)));
			}
			MixBckg(fTrackVec,fBckgCocktail.at(std::get<1>(fIsEventAccepted)),partToMix,gen);

			for (const auto &val : fSignVec)
				hKtSign->Fill(val.kT);
			for (const auto &val : fBckgVec)
				hKtBckg->Fill(val.kT);

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
	
	hAccTracks->Write();
	hM2AccRpc->Write();
	hM2AccTof->Write();
	hPtYRpcProt->Write();
	hPtYTofProt->Write();

	hBetaMomRpcAll->Write();
	hBetaMomRpcProt->Write();
	hBetaMomTofAll->Write();
	hBetaMomTofProt->Write();

	hKtSign->Write();
	hKtBckg->Write();
	for (int i = 0; i < fKtBins; i++)
	{
		hQinvSign.at(i)->Write();
		hQoslSign.at(i)->Write();
		hQinvBckg.at(i)->Write();
		hQoslBckg.at(i)->Write();
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

