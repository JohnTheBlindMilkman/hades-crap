#include "Includes.h"
#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <random>

enum class Detector {RPC, ToF};

struct TrackCandidate
{
	Detector System;
	bool IsAtMdcEdge;
	short int PID, Charge;
	float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass2, Beta, PolarAngle, AzimimuthalAngle;
};

struct EventCandidate
{
	int Centrality;
	short int TargetPlate;
	float ReactionPlaneAngle;
	float X, Y, Z;
	std::vector<TrackCandidate> trackList;
};

struct PairCandidate
{
	TrackCandidate Particle1, Particle2;
	float QInv, QOut, QSide, QLong, Kt, AzimuthalAngle, OpeningAngle, DeltaPhi, DeltaTheta;
};

void CreateEvent(EventCandidate &eventCand, const float &vertx, const float &verty, const float &vertz, const int &cent)
{
	eventCand.Centrality = cent;
	eventCand.trackList.clear();
	eventCand.trackList.resize(0);
	eventCand.X = vertx;
	eventCand.Y = verty;
	eventCand.Z = vertz;
}

bool SelectEvent(EventCandidate &eventCand, const int centIndex = 1, const int nSigma = 1)
{
	// copied from zVertexPeakAndSigma.txt file
	// first: mean, second: stdev
	static const std::vector<std::pair<float,float>> plateVector = {{-54.7598,0.755846},{-51.6971,0.783591},{-47.7996,0.763387},{-44.5473,0.769386},{-40.569,0.781312},{-37.2151,0.762538},{-33.2948,0.76901},{-30.3726,0.742618},{-26.648,0.748409},{-22.5492,0.738462},{-18.9649,0.747727},{-15.5259,0.749724},{-11.8726,0.740386},{-8.45083,0.742672},{-4.58076,0.712394}};
	static const std::pair<float,float> xPosition = {0.1951,0.619};
	static const std::pair<float,float> yPosition = {0.7045,0.6187};

	if (eventCand.Centrality != centIndex)
		return false;

	if ((eventCand.X < (xPosition.first - nSigma*xPosition.second)) || (eventCand.X > (xPosition.first + nSigma*xPosition.second)))
		return false;
	if ((eventCand.Y < (yPosition.first - nSigma*yPosition.second)) || (eventCand.Y > (yPosition.first + nSigma*yPosition.second)))
		return false;

	short int plateNo = 0;
	for(const auto &plate : plateVector)
	{
		if ((eventCand.Z < (plate.first - nSigma*plate.second)) || (eventCand.Z > (plate.first + nSigma*plate.second)))
		{
			eventCand.TargetPlate = plateNo;
			return true;
		}
		plateNo++;
	}

	return false;
}

void CreateTrack(TrackCandidate &trackCand, HParticleCand* particle_cand)
{
	particle_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
	TLorentzVector vecTmp = *particle_cand;

	trackCand.AzimimuthalAngle = particle_cand->getPhi();
	trackCand.Beta = particle_cand->getBeta();
	trackCand.Charge = particle_cand->getCharge();
	trackCand.Energy = vecTmp.E();
	trackCand.IsAtMdcEdge =particle_cand->isAtAnyMdcEdge();
	trackCand.Mass2 = particle_cand->getMass2();
	//trackCand.PID = particle_cand->getGeantPID(); // only for simulations
	trackCand.PolarAngle = particle_cand->getTheta();
	trackCand.Px = vecTmp.Px();
	trackCand.Py = vecTmp.Py();
	trackCand.Pz = vecTmp.Pz();
	trackCand.Rapidity = vecTmp.Rapidity();
	if (particle_cand->getSystem() == 0)
		trackCand.System = Detector::RPC;
	else if (particle_cand->getSystem() == 1)
		trackCand.System = Detector::ToF;
	trackCand.TotalMomentum = vecTmp.P();
	trackCand.TransverseMomentum = vecTmp.Pt();
}

bool SelectTrack(TrackCandidate &trackCand, TCutG *rpcCut, TCutG *tofCut)
{
	//if (trackCand.PID != 14)
		//return false;
	if (trackCand.IsAtMdcEdge)
		return false;
	if (trackCand.TotalMomentum > 2000.)
		return false;
	if (trackCand.TransverseMomentum > 1400.)
		return false;
	if (trackCand.Beta < 0.2)
		return false;

	switch (trackCand.System)
	{
		case Detector::RPC:
			if (rpcCut->IsInside(trackCand.TotalMomentum*trackCand.Charge,trackCand.Beta))
				return true;
			break;

		case Detector::ToF:
			if (tofCut->IsInside(trackCand.TotalMomentum*trackCand.Charge,trackCand.Beta))
				return true;
			break;
	}
	
	return false;
}

PairCandidate CFKinematics(const TrackCandidate &part1, const TrackCandidate &part2)
{
	PairCandidate pair;
	pair.Particle1 = part1;
	pair.Particle2 = part2;

	// adapted from https://github.com/DanielWielanek/HAL/blob/main/analysis/femto/base/FemtoPairKinematics.cxx
	Double_t tPx = part1.Px + part2.Px;
    Double_t tPy = part1.Py + part2.Py;
    Double_t tPz = part1.Pz + part2.Pz;
    Double_t tE = part1.Energy + part2.Energy;
    Double_t tPt = tPx * tPx + tPy * tPy;
    Double_t tMt = tE * tE - tPz * tPz;  // mCVK;
    tMt = TMath::Sqrt(tMt);
    pair.Kt = TMath::Sqrt(tPt);
    Double_t tBeta  = tPz / tE;
    Double_t tGamma = tE / tMt;

    // Transform to LCMS

    Double_t particle1lcms_pz = tGamma * (part1.Pz - tBeta * part1.Energy);
    Double_t particle1lcms_e  = tGamma * (part1.Energy - tBeta * part1.Pz);
    Double_t particle2lcms_pz = tGamma * (part2.Pz - tBeta * part2.Energy);
    Double_t particle2lcms_e  = tGamma * (part2.Energy - tBeta * part2.Pz);

    // Rotate in transverse plane

    Double_t particle1lcms_px = (part1.Px * tPx + part1.Py * tPy) / pair.Kt;
    Double_t particle1lcms_py = (-part1.Px * tPy + part1.Py * tPx) / pair.Kt;

    Double_t particle2lcms_px = (part2.Px * tPx + part2.Py * tPy) / pair.Kt;
    Double_t particle2lcms_py = (-part2.Px * tPy + part2.Py * tPx) / pair.Kt;

    pair.QOut = abs(particle1lcms_px - particle2lcms_px);
    pair.QSide = abs(particle1lcms_py - particle2lcms_py);
    pair.QLong = abs(particle1lcms_pz - particle2lcms_pz);
    Double_t mDE = particle1lcms_e - particle2lcms_e;
    pair.QInv = TMath::Sqrt(TMath::Abs(pair.QOut * pair.QOut + pair.QSide * pair.QSide + pair.QLong * pair.QLong - mDE * mDE));

	return pair;
}

float CalcOpeningAngle(const TrackCandidate &part1, const TrackCandidate &part2)
{
	float ptot = sqrt((part1.Px*part1.Px + part1.Py*part1.Py + part1.Pz*part1.Pz) * (part2.Px*part2.Px + part2.Py*part2.Py + part2.Pz*part2.Pz));
	if (ptot <= 0.)
	{	
		return 0.;
	}
	else
	{
		float arg = (part1.Px*part2.Px + part1.Py*part2.Py + part1.Pz*part2.Pz) / ptot;
		if (arg > 1.) return 1.;
		if (arg < -1.) return -1.;
		return acos(arg);
	}
}

void CreatePair(PairCandidate &pairCand, const TrackCandidate &part1, const TrackCandidate &part2)
{
	pairCand = CFKinematics(part1,part2);
	pairCand.OpeningAngle = CalcOpeningAngle(part1,part2);
	pairCand.AzimuthalAngle = (part1.AzimimuthalAngle + part2.AzimimuthalAngle) / 2.;
	pairCand.DeltaPhi = part1.AzimimuthalAngle - part2.AzimimuthalAngle;
	pairCand.DeltaTheta = part1.PolarAngle - part2.PolarAngle;
}

bool SelectPair(PairCandidate &pairCand, const float openAngle = 0.0, const float minKt = 150., const float maxKt = 1650.)
{
	std::vector<PairCandidate> acceptedPairs;

	if (pairCand.OpeningAngle <= openAngle)
		return false;
	if (pairCand.Kt < minKt || pairCand.Kt > maxKt)
		return false;

	return true;
}

void CalcSignal(std::vector<PairCandidate> &pairVec, const EventCandidate &eventCand, bool reorder, float openAngle = 0.0)
{
	PairCandidate pairCand;
	for (size_t iter1 = 0; iter1 < eventCand.trackList.size(); iter1++)
		for (size_t iter2 = iter1+1; iter2 < eventCand.trackList.size(); iter2++)
		{
			if (reorder)
				CreatePair(pairCand,eventCand.trackList.at(iter2),eventCand.trackList.at(iter1));
			else
				CreatePair(pairCand,eventCand.trackList.at(iter1),eventCand.trackList.at(iter2));
			if (SelectPair(pairCand,openAngle))
				pairVec.push_back(pairCand);
		}		
}

void CalcBackground(std::vector<PairCandidate> &pairVec, const EventCandidate &eventCand, const std::deque<TrackCandidate> &bckgVec, float openAngle = 0.0)
{
	PairCandidate pairCand;
	for (auto &track : eventCand.trackList)
		for (auto &bckg : bckgVec)
		{
			CreatePair(pairCand,track,bckg);
			if (SelectPair(pairCand,openAngle))
				pairVec.push_back(pairCand);
		}
				
}

void MixBckg(const EventCandidate &eventCand, std::deque<TrackCandidate> &bckgVec, const int &maxSize, std::mt19937 &generator)
{
	std::uniform_int_distribution<int> dist(0,eventCand.trackList.size()-1);
	bckgVec.push_back(eventCand.trackList.at(dist(generator)));
	if (bckgVec.size() > maxSize)
		bckgVec.pop_front();
}

int AssignKt(float kt)
{
	static std::vector<std::pair<float,float> > ktBinVec = {{150,450},{450,750},{750,1050},{1050,1350},{1350,1650}};
	for (std::size_t i = 0; i < ktBinVec.size(); i++)
		if (kt >= ktBinVec.at(i).first && kt <= ktBinVec.at(i).second)
			return i;

	return -999;
}

float DegToRad(const int &angle)
{
	return (TMath::Pi()/180)*angle;
}

bool IntToBool(const int &integ)
{
	return integ;
}

void FillAndClear(std::vector<PairCandidate> &pairVec, std::array<TH1D*,5> hInp1, std::array<TH2D*,5> hInp2, std::array<TH3D*,5> hInp3)
{
	int ktBin = 0;
	for (auto &pair : pairVec)
	{
		ktBin = AssignKt(pair.Kt);
		if (ktBin == -999)
			continue;

		hInp1.at(ktBin)->Fill(pair.QInv);
		hInp2.at(ktBin)->Fill(pair.DeltaPhi,pair.DeltaTheta);
		hInp3.at(ktBin)->Fill(pair.QOut,pair.QSide,pair.QLong);
	}

	pairVec.clear();
	pairVec.resize(0);
}

int baseAnalysis(TString inputlist = "", TString outfile = "baseOutFile.root", Long64_t nDesEvents = -1, Int_t maxFiles = -1)	//for simulation set approx 100 files and output name testOutFileSim.root
{
	const int fTargetPlates = 15; // number of target plates
	std::vector<std::pair<float,float> > fKtBinVec = {{150,450},{450,750},{750,1050},{1050,1350},{1350,1650}};
	const int fKtBins = 5; // number of kT bins
	const int fAlphas = 6; // number of opening angle cuts
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
		//TString inputFolder = "/lustre/hades/dstsim/apr12/gen9vertex/no_enhancement_gcalor/root/";
		
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

	std::array<TH1D*,fKtBins> hQinvSign,hQinvBckg;
	std::array<TH2D*,fKtBins> hDphiDthetaSign,hDphiDthetaBckg;
	std::array<TH3D*,fKtBins> hQoslSign,hQoslBckg;
	for (int i = 0; i < fKtBins; i++)
	{
		hQinvSign.at(i) = new TH1D(TString::Format("hQinvSign%d",i),"Signal of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
		hDphiDthetaSign.at(i) = new TH2D(TString::Format("hDphiDthetaSign%d",i), "#Delta#phi vs #Delta#theta distribution of signal 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90);
		hQoslSign.at(i) = new TH3D(TString::Format("hQoslSign%d",i),"Signal of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);
		hQinvBckg.at(i) = new TH1D(TString::Format("hQinvBckg%d",i),"Background of Protons 0-10%% centrality;q_{inv} [MeV];CF(q_{inv})",750,0,3000);
		hDphiDthetaBckg.at(i) = new TH2D(TString::Format("hDphiDthetaBckg%d",i), "#Delta#phi vs #Delta#theta distribution of background 0-10%%;#Delta#phi [deg]; #Delta#theta [deg]",360,-360,360,90,-90,90);
		hQoslBckg.at(i) = new TH3D(TString::Format("hQoslBckg%d",i),"Background of Protons 0-10%% centrality;q_{out} [MeV];q_{side} [MeV];q_{long} [MeV];CF(q_{inv})",250,0,1000,250,0,1000,250,0,1000);
	}

	TFile *cutfile_betamom_pionCmom = new TFile("/lustre/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
	TCutG* betamom_2sig_p_tof_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_TOF_2.0");
	TCutG* betamom_2sig_p_rpc_pionCmom = cutfile_betamom_pionCmom->Get<TCutG>("BetaCutProton_RPC_2.0");
	
	//create RNG
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> uniInt(0,1);
	
	// create objects for mixing
	EventCandidate fEvent;
	TrackCandidate fTrack;
	std::array<std::deque<TrackCandidate>,fTargetPlates> fBckgCocktail;
	std::vector<PairCandidate> fSignVec, fBckgVec;
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
		CreateEvent(fEvent,vertX,vertY,vertZ,centClassIndex);

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
		
		if (!SelectEvent(fEvent)) 
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
			CreateTrack(fTrack,particle_cand);
			
			//================================================================================================================================================================
			// Put your analyses on track level here
			//================================================================================================================================================================
			
			if (fTrack.System == Detector::RPC)
			{
				// fill RPC monitors for all tracks
			}
			else
			{
				// fill ToF monitors for all tracks
			}

			if (!SelectTrack(fTrack,betamom_2sig_p_rpc_pionCmom,betamom_2sig_p_tof_pionCmom))
				continue;
			fEvent.trackList.push_back(fTrack);
			hPhiTheta->Fill(fTrack.AzimimuthalAngle,fTrack.PolarAngle);

			if (fTrack.System == Detector::RPC)
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
			if (fBckgCocktail.at(fEvent.TargetPlate).size())
			{
				CalcSignal(fSignVec,fEvent,IntToBool(uniInt(gen)),DegToRad(0.0));
				CalcBackground(fBckgVec,fEvent,fBckgCocktail.at(fEvent.TargetPlate),DegToRad(0.0));
			}
			MixBckg(fEvent,fBckgCocktail.at(fEvent.TargetPlate),partToMix,gen);

			FillAndClear(fSignVec,hQinvSign,hDphiDthetaSign,hQoslSign);
			FillAndClear(fBckgVec,hQinvBckg,hDphiDthetaBckg,hQoslBckg);
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

	for (int i = 0; i < fKtBins; i++)
	{
		hQinvSign.at(i)->Write();
		hDphiDthetaSign.at(i)->Write();
		hQoslSign.at(i)->Write();
		hQinvBckg.at(i)->Write();
		hDphiDthetaBckg.at(i)->Write();
		hQoslBckg.at(i)->Write();
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

