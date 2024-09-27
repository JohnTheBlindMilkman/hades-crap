#include "Includes.h" // A lot of header files from Root and Hydra (general ones)
#include "GeomFunct.h" // DCA stuff
#include "emcdef.h"

#include <vector>
#include <algorithm>
#include <cmath>

#include "henergylosscorrpar.h" //energy loss correction

using namespace std;

Int_t getTargetPlateNumber(Double_t eventZvertex){ //for mar19 target
	if(eventZvertex >= -63 && eventZvertex < -60) return 1;
    else if(eventZvertex >= -60 && eventZvertex < -56.5) return 2;
    else if(eventZvertex >= -56.5 && eventZvertex < -53) return 3;
    else if(eventZvertex >= -53 && eventZvertex < -49) return 4;
    else if(eventZvertex >= -49 && eventZvertex < -45.5) return 5;
    else if(eventZvertex >= -45.5 && eventZvertex < -42) return 6;
    else if(eventZvertex >= -42 && eventZvertex < -39) return 7;
    else if(eventZvertex >= -39 && eventZvertex < -35.5) return 8;
    else if(eventZvertex >= -35.5 && eventZvertex < -31.5) return 9;
    else if(eventZvertex >= -31.5 && eventZvertex < -28) return 10;
    else if(eventZvertex >= -28 && eventZvertex < -24.5) return 11;
    else if(eventZvertex >= -24.5 && eventZvertex < -21) return 12;
    else if(eventZvertex >= -21 && eventZvertex < -16.5) return 13;
    else if(eventZvertex >= -16.5 && eventZvertex < -13.5) return 14;
    else if(eventZvertex >= -13.5 && eventZvertex < -11) return 15;
    else return 1000;
}

struct DParticle {
  TLorentzVector vec;
  Int_t geantTrackNum;
};
 
//================================================================================
// Main Function of the Macro with 4 parameters
// 1. Comma separated DST file list to run the macro on the cluster
// 2. Name / Path of the outputfile to be created
// 3. Amount of events to process (-1 means all events)
// 4. If not -1 amount of random files to read instead of files from file list parameter
//================================================================================
Int_t MomentumResolution(TString inputlist = "", TString outfile = "test.root", Long64_t nDesEvents = -1, Int_t maxFiles = 1){
    //--------------------------------------------------------------------------------
    // Initialization of the global ROOT object and the Hades Loop
    // The Hades Loop used as an interface to the DST data (Basically a container of a TChain).
    // kTRUE - The global HADES object is being created if not existing
    //--------------------------------------------------------------------------------
	
    gROOT -> SetBatch(kTRUE);
    gStyle -> SetOptStat(111111);
    TROOT dst_analysis("DstAnalysisMacro", "Simple DST analysis Macro");
    HLoop* loop = new HLoop(kTRUE);
    TH1::SetDefaultSumw2();

    //--------------------------------------------------------------------------------
    // The following block finds / adds the input DST files to the HLoop
    //--------------------------------------------------------------------------------
    if (maxFiles == -1)
        loop->addMultFiles(inputlist);//use instead of addFiles if run on batch farm
    else {
        Int_t nFiles = 0;
        TString inputFolder = "/lustre/nyx/hades/dstsim/mar19/ag1580ag/gen5/bmax8/no_enhancement_gcalor/root"; //sim data
        TSystemDirectory* inputDir = new TSystemDirectory("inputDir", inputFolder);
        TList* files = inputDir->GetListOfFiles();

        for (Int_t i = 0; i <= files->LastIndex() && nFiles < maxFiles; i++){
            
            if (((TSystemFile*) files->At(i))->IsDirectory()) continue;
                loop->addFile(inputFolder + "/" + ((TSystemFile*) files->At(i))->GetName());
                nFiles++;
        }
    }
    //--------------------------------------------------------------------------------
    // Booking the categories to be read from the DST files.
    // By default all categories are booked therefore -* (Unbook all) first and book the ones needed
    // All required categories have to be booked except the global Event Header which is always booked
    //--------------------------------------------------------------------------------
    if (!loop->setInput("-*,+HParticleEvtInfo,+HParticleCandSim,+HGeantKine"))//,+HEmcClusterSim, +HEmcCal, +HWallHit"))
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
    // Creating the placeholder variables to read data from categories and getting categories (ThEy have to be booked!)
    //--------------------------------------------------------------------------------
    HParticleCandSim* particleCand;
    HEventHeader* eventHeader;
    HParticleEvtInfo* particleInfo;

    HCategory* particleInfoCat = (HCategory*) HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* particleCandCat = (HCategory*) HCategoryManager::getCategory(catParticleCand);
    HCategory* kineCandCat = (HCategory*) HCategoryManager::getCategory(catGeantKine);
    // HCategory* emcClusterCat = (HCategory *)  HCategoryManager::getCategory(catEmcCluster);
    
    if (!particleCandCat || !kineCandCat)//|| !emcClusterCat) // If the category for the reconstructed trackes does not exist the macro makes no sense
    exit(1);
	
    //================================================================================================================================================================
    // Object declatarions (histograms, trees ect.)
    //================================================================================================================================================================
    TFile* outFile = new TFile(outfile.Data(), "RECREATE");

    TH2D* hMomResolution_Lambda = new TH2D("hMomResolution","Mom difference; Mom(reco); #DeltaMom(kine-reco)",20,0,1000,200,-100, 100);
    TH2D* hPhiResolution_Lambda = new TH2D("hPhiResolution","#phi difference; E(reco); #Delta#phi(kine-reco)",20,0,1000,200,-0.2, 0.2);
    TH2D* hThetaResolution_Lambda = new TH2D("hThetaResolution","#theta difference; E(reco); #Delta#theta(kine-reco)",20,0,1000,200,-0.053, 0.053);

    TH2D* hMomResolution_Deuteron = new TH2D("hMomResolution_deuteron","Mom difference deuteron; Mom(reco); #DeltaMom(kine-reco)",20,0,1000,200,-100, 100);
    TH2D* hPhiResolution_Deuteron = new TH2D("hPhiResolution_deuteron","#phi difference deuteron; E(reco); #Delta#phi(kine-reco)",20,0,1000,200,-0.2, 0.2);
    TH2D* hThetaResolution_Deuteron = new TH2D("hThetaResolution_deuteron","#theta difference deuteron; E(reco); #Delta#theta(kine-reco)",20,0,1000,200,-0.053, 0.053);
 
    //------------------------------------------------------------------------------------------------------
    // Loading & setting Cut files
    //------------------------------------------------------------------------------------------------------

    //energy loss correction
    HEnergyLossCorrPar enLossCorr;
    enLossCorr.setDefaultPar("mar19");

    TFile *file = new TFile("/lustre/hades/user/mstefan/sub/loopDST/Mar19AgAg1580.root");
  //----------------2-Sigma----------------------------------
  TCutG *protonRPCCut = (TCutG*)file->Get("tcgPBetaProtonRPC2Sig");
  TCutG *protonToFCut = (TCutG*)file->Get("tcgPBetaProtonToF2Sig");
  TCutG *piNRPCCut = (TCutG*)file->Get("tcgPBetaPiMRPC2Sig");
  TCutG *piNToFCut = (TCutG*)file->Get("tcgPBetaPiMToF2Sig");

  TFile *file2 = new TFile("/lustre/hades/user/mstefan/sub/loopDST/Mar19AgAg1580_Gen5.root");
  TCutG *dRPCCut = (TCutG*)file2 -> Get("tcgPdEdxMDCDeuteronRPC2Sig");
  TCutG *dToFCut = (TCutG*)file2 -> Get("tcgPdEdxMDCDeuteronTOF2Sig");

  vector<DParticle> LambdaVector;

  vector<TLorentzVector> usedProtonRPC;
  vector<TLorentzVector> usedpionNRPC;
  vector<TLorentzVector> HBTLambdaRPC;
  vector<TLorentzVector> HBTProtonRPC;

    vector<HParticleCandSim*> protons;
    vector<HParticleCandSim*> pions;
    vector<HParticleCandSim*> deuterons; 

	 
	//-------------------------------------------------------------------------------
    // Event & track counter
    //-------------------------------------------------------------------------------
    enum Counters {
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
    // track sorter and a simple stopwatch object
    //--------------------------------------------------------------------------------
    HParticleTrackSorter sorter;
    sorter.setIgnoreInnerMDC(); //do not reject double inner MDC hits
    sorter.init();

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    //--------------------------------------------------------------------------------
    //  Event chara setup
    //--------------------------------------------------------------------------------
    HParticleEvtChara evtChara;
    TString ParameterfileCVMFS = "/lustre/nyx/hades/user/bkardan/param/centrality_epcorr_mar19_ag158ag_3200A_glauber_gen5_pass3_2021_08.root";
    if(!evtChara.setParameterFile(ParameterfileCVMFS)){
      cout << "Parameterfile not found !!! " << endl;return kFALSE;
    }
    if(!evtChara.init()) {
      cout << "HParticleEvtChara not init!!! " << endl;return kFALSE;
    }

    Int_t eCentEst   = HParticleEvtChara::kTOFRPC;
    Int_t eCentClass = HParticleEvtChara::k10; //k:10 - 10% intervals, k:5 - 5%

    evtChara.printCentralityClass(eCentEst, eCentClass);
    //--------------------------------------------------------------------------------
    // The amount of events to be processed
    //--------------------------------------------------------------------------------
    Long64_t nEvents = loop->getEntries();
    if (nDesEvents >= 0 && nEvents > nDesEvents)
    nEvents = nDesEvents;

    //============================================================================================================================
    // ---EVENT LOOP---
    //============================================================================================================================
    for (Long64_t event = 0; event < nEvents; event++) {
        if (loop->nextEvent(event) <= 0) {
            cout << " Last events processed " << endl;
            break;
        }
    
        //--------------------------------------------------------------------------------
        // Just the progress of the analysis
        //--------------------------------------------------------------------------------
        HTool::printProgress(event, nEvents, 1, "Analyzed events: ");

        //-------------------------------------------------------------------------------------
        // Getting the amount of tracks, event header, Particle Event Info & global event vertex
        //------------------------------------------------------------------------------------
        Int_t nTracks = particleCandCat->getEntries();
        eventHeader = gHades->getCurrentEvent()->getHeader();
        particleInfo = HCategoryManager::getObject(particleInfo, particleInfoCat, 0);
        Float_t EventVertexZ = eventHeader->getVertexReco().getZ(); 

        HGeomVector EventVertex = eventHeader->getVertexReco().getPos();
    
        hCounter->Fill(cNumAllEvents);
        //====================================================================================
        // Centrality & plate selection
        //====================================================================================
        Int_t centrality = evtChara.getCentralityClass(eCentEst, eCentClass);
        Int_t plate = getTargetPlateNumber(EventVertexZ);
        
        if(centrality > 4 || centrality < 1) continue; 
        if(plate==1000) continue; //skip wrong plates

        hCounter->Fill(cNumSelectedEvents); 

        //--------------------------------------------------------------------------------
        // Sorter setup for event
        //--------------------------------------------------------------------------------
        sorter.cleanUp();
        sorter.resetFlags(kTRUE, kTRUE, kTRUE, kTRUE);
        sorter.fill(HParticleTrackSorter::selectHadrons);

        sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);  
        // sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE); // reset all flags for flags (0-28) ,reject,used,lepton
        // sorter.fill(HParticleTrackSorter::selectLeptons); // fill only good leptons

        // sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsLeptonSorter);
        // sorter.fill(HParticleTrackSorter::selectHadrons); // fill only good hadrons (already marked good leptons will be skipped)

        // sorter.selectBest(Particle::kIsBestRKSorter, Particle::kIsHadronSorter);

        //----------------------------------------
        // Event variables
        //----------------------------------------
        // Int_t firedEcalModules[6][163]= {0};

        protons.clear();
        pions.clear();
        deuterons.clear();
        LambdaVector.clear();

        //============================================================================================================================
        //  TRACK LOOP (particleCands) 
        //============================================================================================================================
        for (Int_t track = 0; track<nTracks; track++) {
            particleCand = HCategoryManager::getObject(particleCand, particleCandCat, track);
            
            Float_t dEdx = particleCand -> getMdcdEdx();
            Short_t chargeCand = particleCand -> getCharge();
            Float_t momCand = particleCand -> getMomentum();
            Float_t betaCand = particleCand ->getBeta();
            Float_t theta = particleCand->getTheta();
            Float_t mass     = particleCand->getMass();
            Short_t sysCand    = particleCand -> getSystemUsed(); // detector: 0 RPC, 1 ToF

            //--------------------------------------------------------------------------------
            // Discarding all tracks that have been discarded by the track sorter and counting all / good tracks
            //--------------------------------------------------------------------------------
            hCounter->Fill(cNumAllTracks);
            if (!particleCand->isFlagBit(Particle::kIsUsed))
            continue;
    
            hCounter->Fill(cNumSelectedTracks);

            //gather pi- & proton, check if same parent (Lambda), get parent 4vector

            if (protonRPCCut->IsInside(chargeCand*momCand, betaCand) || protonToFCut->IsInside(chargeCand*momCand, betaCand)){ // RPC+ToF

                particleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
                Double_t momCandproton = enLossCorr.getCorrMom(14,momCand,theta);
                particleCand->setMomentum(momCandproton);
                protons.push_back(particleCand);
            }
            else if (piNRPCCut->IsInside(chargeCand*momCand, betaCand) || piNToFCut->IsInside(chargeCand*momCand, betaCand)) { // RPC+ToF

                particleCand->calc4vectorProperties(HPhysicsConstants::mass(9));
                Double_t momCandpion = enLossCorr.getCorrMom(9,momCand,theta);
                particleCand->setMomentum(momCandpion);
                pions.push_back(particleCand);
            }

            if( ( mass > 1600 && mass < 2150) && (
            (sysCand == 0 && dRPCCut->IsInside(chargeCand*momCand, dEdx)) ||
            (sysCand == 1 && dToFCut->IsInside(chargeCand*momCand, dEdx)))){

                particleCand->setMomentum(particleCand->getCorrectedMomentumPID(45));
                particleCand -> calc4vectorProperties(HPhysicsConstants::mass(45)); // d

                if (particleCand->getGeantPID() == 45)
                    deuterons.push_back(particleCand);
            }

            //calcSegVector() for pi- & proton, check their ID and parent geant track number, get their TLorentz's, add them
		
            // //Geant addition
            // geantTrack=particleCand->getGeantTrack(); 
            // kineCand = (HGeantKine*)kineCandCat->getObject(geantTrack -1); //because geant numbering scheme it need -1 here...
            // idKine = kineCand->getGeantPID();

            // HGeantKine* kineParent = kineCand->getgeantParentTrackNum(geantTrack);
            // if(kineParent!=NULL){
            //     idParent = kineParent ->getGeantPID();   
            //     trackParent = kineParent -> getGeantTrack();
            // }


            // if (idParent!=18) continue; //18 = lambda ID?

            // Double_t energyReco, energyGen;
            // Double_t phiReco, phiGen;
            // Double_t thetaReco, thetaGen;


            // energyReco = particleCand->getEnergy();
            // phiReco = particleCand->getPhi() * TMath::DegToRad(); //0-360 degrees
            // thetaReco = particleCand->getTheta() * TMath::DegToRad();
            // energyGen = kineCand -> getE();
            // phiGen = kineCand -> getPhiDeg() * TMath::DegToRad();
            // thetaGen = kineCand -> getThetaDeg() * TMath::DegToRad();

            // if(phiGen>TMath::Pi()) phiGen -= 2 * TMath::Pi();

            // hEnergyResolution_Lambda -> Fill(energyReco, energyGen - energyReco);
            // hPhiResolution_Lambda -> Fill(energyReco, phiGen - phiReco);
            // hThetaResolution_Lambda -> Fill(energyReco, thetaGen - thetaReco);

                        
        } // End of particleCand track loop

        for (vector<HParticleCandSim*>::iterator itDau1 = protons.begin(); itDau1 != protons.end(); itDau1++) {
            HParticleCandSim* hpcDau1 = *itDau1;

            HGeomVector hgvBaseDau1, hgvDirDau1;
            HParticleTool::calcSegVector(hpcDau1->getZ(), hpcDau1->getR(), (TMath::DegToRad()*hpcDau1->getPhi()), (TMath::DegToRad()*hpcDau1->getTheta()), hgvBaseDau1, hgvDirDau1);

            Double_t VerDistA = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvBaseDau1, hgvDirDau1, EventVertex);

            for (vector<HParticleCandSim*>::iterator itDau2 = pions.begin(); itDau2 != pions.end(); itDau2++) {
                HParticleCandSim* hpcDau2 = *itDau2;

                HGeomVector hgvBaseDau2, hgvDirDau2;
                HParticleTool::calcSegVector(hpcDau2->getZ(), hpcDau2->getR(), (TMath::DegToRad()*hpcDau2->getPhi()), (TMath::DegToRad()*hpcDau2->getTheta()), hgvBaseDau2, hgvDirDau2);

                Double_t VerDistB = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvBaseDau2, hgvDirDau2, EventVertex);

                HGeomVector hgvDecayVertex = HParticleTool::calcVertexAnalytical(hgvBaseDau1, hgvDirDau1, hgvBaseDau2, hgvDirDau2);

                Double_t VerDistX = (hgvDecayVertex - EventVertex).length();

                Double_t Alpha = TMath::RadToDeg()*TMath::ACos(hgvDirDau1.scalarProduct(hgvDirDau2) / (hgvDirDau1.length() * hgvDirDau2.length()));

                TLorentzVector tlvMot = *hpcDau1 + *hpcDau2;
                TLorentzVector tlvMot1 = *hpcDau1;
                TLorentzVector tlvMot2 = *hpcDau2;
                HGeomVector hgvDirMot(tlvMot.Px(), tlvMot.Py(), tlvMot.Pz());

                Double_t VerDistC = HParticleTool::calculateMinimumDistanceStraightToPoint(hgvDecayVertex, hgvDirMot, EventVertex);

                Double_t MinTrackDist = HParticleTool::calculateMinimumDistance(hgvBaseDau1, hgvDirDau1, hgvBaseDau2, hgvDirDau2);

                if(VerDistX > 65 && VerDistA > 8 && VerDistB > 24 && VerDistC < 5 && MinTrackDist < 6 && Alpha > 15){
                    
                    DParticle Lambda;
                    double LambdaMass = tlvMot.M();          

                    if(hpcDau1->getGeantPID()==14 && hpcDau2->getGeantPID()==9){ //checking if we have true proton & pion

                        if(hpcDau1->getGeantParentPID()==18 && hpcDau2->getGeantParentPID()==18 && (hpcDau1->getGeantParentTrackNum()==hpcDau2->getGeantParentTrackNum())){
                            
                            if(LambdaMass > 1110 && LambdaMass < 1120){
                                usedProtonRPC.push_back(tlvMot1);
                                usedpionNRPC.push_back(tlvMot2);

                                Lambda.vec = tlvMot;
                                Lambda.geantTrackNum = hpcDau1->getGeantParentTrackNum();
                                LambdaVector.push_back(Lambda);

                            } //Mass cut
                        }
                    }
                } // DCA cuts
            } // End of second dauther track loop
        } // first daugther track loop

        for (std::size_t i = 0; i< LambdaVector.size(); i++)
        {
            DParticle particle = LambdaVector[i];
            HGeantKine* kine = (HGeantKine*)kineCandCat->getObject(particle.geantTrackNum -1);

            Double_t momReco = particle.vec.P();
            Double_t phiReco = particle.vec.Phi();
            Double_t thetaReco = particle.vec.Theta();

            Double_t momKine = kine->getTotalMomentum();
            Double_t phiKine = kine->getPhiDeg()*TMath::DegToRad(); //0-2pi
            Double_t thetaKine = kine->getThetaDeg()*TMath::DegToRad();

            if(phiKine>TMath::Pi()) phiKine-=2*TMath::Pi();

            hMomResolution_Lambda -> Fill(momReco, momKine - momReco);
            hPhiResolution_Lambda -> Fill(phiReco, phiKine - phiReco);
            hThetaResolution_Lambda -> Fill(thetaReco, thetaKine - thetaReco);

        }

        for(vector<HParticleCandSim*>::iterator itDeut = deuterons.begin(); itDeut != deuterons.end(); itDeut++) {
          HParticleCandSim* hpcDeut = *itDeut;

          Int_t geantTrack=particleCand->getGeantTrack(); 
          HGeantKine* kine = (HGeantKine*)kineCandCat->getObject(geantTrack -1); //because geant numbering scheme it need -1 here...

          Double_t momReco = hpcDeut->getMomentum();
          Double_t phiReco = hpcDeut->getPhi()*TMath::DegToRad();
          Double_t thetaReco = hpcDeut->getTheta()*TMath::DegToRad();

          Double_t momKine = kine->getTotalMomentum();
          Double_t phiKine = kine->getPhiDeg()*TMath::DegToRad(); //0-2pi
          Double_t thetaKine = kine->getThetaDeg()*TMath::DegToRad();

          hMomResolution_Deuteron -> Fill(momReco, momKine - momReco);
          hPhiResolution_Deuteron -> Fill(phiReco, phiKine - phiReco);
          hThetaResolution_Deuteron -> Fill(thetaReco, thetaKine - thetaReco);
          
        } // end for for Deuterons

    }// End of event loop

    //--------------------------------------------------------------------------------
    // Doing some cleanup and finalization work
    //--------------------------------------------------------------------------------
    sorter.finalize();
    timer.Stop();

    LambdaVector.clear();
    protons.clear();
    pions.clear();
    deuterons.clear();

    //--------------------------------------------------------------------------------
    // storing results
    //--------------------------------------------------------------------------------
    
    outFile->cd();

    hMomResolution_Lambda -> Write();
    hPhiResolution_Lambda -> Write();
    hThetaResolution_Lambda -> Write();

    hMomResolution_Deuteron -> Write();
    hPhiResolution_Deuteron -> Write();
    hThetaResolution_Deuteron -> Write();

	outFile->Save();
	outFile->Close();

    //--------------------------------------------------------------------------------
    // Closing file and finalization
    //--------------------------------------------------------------------------------

    cout << "####################################################" << endl;
    cout << "  Analyse finished succesfully! Go back to work!" << endl;
    cout << "####################################################" << endl;
    return 0;
}

