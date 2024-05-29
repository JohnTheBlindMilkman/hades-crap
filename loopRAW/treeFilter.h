

#include "hparticlecand.h"
#include "hcategorymanager.h"
#include "hcategory.h"

#include "hparticletree.h"
#include "htime.h"

#include "hparticledef.h"
#include "hstartdef.h"
#include "walldef.h"
#include "richdef.h"
#include "showerdef.h"
#include "hstartdef.h"
#include "tofdef.h"
#include "rpcdef.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"


#include "TObjArray.h"
#include "TSystem.h"

////////////////////////////////////////////////////////////////////////////////
//
//
// HParticleTree
//
// Reconstructor to create an filtered root output.
// The output is generated independend for each HParticleTree
// reconstructor and das not affect Hades. The purpose
// is to create different output files with different event
// structure in parallel and allow to reduce the data volume
// for special analysis tasks (for example rare dilepton decays
// in Au+Au collisions).
// The filtering starts from HParticleCand. Only objects
// which are flagged kIsUsed are copied to the output.
// The user can specify selectLeptons() and selectHadrons()
// functions. By default the functions of HParticleTracksorter
// are used (for the documentation of the selection functions see
// HParticleTrackSorter). In addition the user can specify a selectEvent()
// function function which allows to check the full event.
//
// The program works the following way:
//
// I.    The original flags of the HParticleCand objects and HParticleSorter
//       are backupted. This alows to run the task without affecting other tasks.
// II.   In the input HParticleCand objects are flagged using the user specified
//       selection functions or default functions.
// III.  If the user provides a Bool_t selectEvent(TObjArray*) function the
//       properties of the full event can be evaluated. If the function returns
//       kFALSE an empty event ist stored to preserve the 1 to 1 correlation
//       with the dst file.
//       With setSkipEmptyEvents (kTRUE) skipping of not selected events can be forced
// III.  From HParticleCand objects which are flagged kIsUsed the selection
//       procedure starts. Only flagged objects are considered and the detector
//       hits belonging to the candidates are stored in the output too if enabled.
//       With setSkipTracks(kFALSE) all tracks will forced to the output.
//       HParticleEvtInfo, HStart2Cal, HStartHit and HWallHit are fully copied
//       if enabled (see list below).
//       Alternatively (not at the same time!) to using flagged candidates as
//       selection the user can provide a function to select the kept particles by
//       setUserkeepTrack(Bool_t (*function)(HParticleCand* )). The function
//       has to return kTRUE for selected candidates.
//       All indices are resorted. The new event will look like the old one,
//       but containing only "good" candidates. The full features of the analysis
//       of the DSTs are preserved.

Bool_t selectProtonsTreeFilter(HParticleCand* pcand){
   // we want only protons in this case
   Bool_t test = kFALSE;
   
   test = pcand->isFlagAND(4,
   			Particle::kIsAcceptedHitInnerMDC,
   			Particle::kIsAcceptedHitOuterMDC,
   			Particle::kIsAcceptedHitMETA,
   			Particle::kIsAcceptedRK
   		       );
   
   if(test) test = (pcand->getPID() == 14)  ? kTRUE: kFALSE; // only identified protons
   
   return test;
}

Bool_t selectPionsTreeFilter(HParticleCand* pcand){
   // we want only leptons in this case
   //crg271117 my filter
   Bool_t test = kFALSE;
   
   //Double_t beta = pcand->getBeta();
   //Double_t mom  = pcand->getMomentum();
   //Int_t system  = pcand->getSystemUsed();
   
   //if (mom > 1500) return test;
   
   test = pcand->isFlagAND(4,
   			Particle::kIsAcceptedHitInnerMDC,
   			Particle::kIsAcceptedHitOuterMDC,
   			Particle::kIsAcceptedHitMETA,
   			Particle::kIsAcceptedRK
   		       );
   
   
   //if(test) test = pcand->getChi2() < 1000 ? kTRUE: kFALSE;
   if(test) test = (pcand->getPID() == 8 || pcand->getPID() == 9)  ? kTRUE: kFALSE; // only identified pions
   
   //if(!test) return kFALSE;
   
   //if(test) test = pcand->getBeta() > 0.9          ? kTRUE : kFALSE ;
   //if(test) test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE ;
   
   return test;
}

Bool_t selectLeptonsTreeFilter(HParticleCand* pcand)
{
    // build in selection function for lepton candidates.
    // Requires besides an RICH hit, RK + META and fitted
    // inner+outer segment.

    Bool_t test = kFALSE;

    

    Double_t beta = pcand->getBeta();
    Double_t mom  = pcand->getMomentum();
    Int_t system  = pcand->getSystemUsed();

    if (mom > 1500) return test;

    if((system == 0 && beta > 0.94 &&  mom < 200) ||
       (system == 1 && beta > 0.96 &&  mom < 250)
      )
    {   // for low momenta do not use RICH
	test = pcand->isFlagAND(4,
				Particle::kIsAcceptedHitInnerMDC,
				Particle::kIsAcceptedHitOuterMDC,
				Particle::kIsAcceptedHitMETA,
				Particle::kIsAcceptedRK
			       );
    } else { // for large momenta check RICH

	test = pcand->isFlagAND(5,
				Particle::kIsAcceptedHitRICH,
				Particle::kIsAcceptedHitInnerMDC,
				Particle::kIsAcceptedHitOuterMDC,
				Particle::kIsAcceptedHitMETA,
				Particle::kIsAcceptedRK
			       );

	if(test) test = fabs(pcand->getDeltaTheta()) < 4 ? kTRUE : kFALSE ;
	if(test) test = fabs(pcand->getDeltaPhi())   < 4 ? kTRUE : kFALSE ;
    }

    if(test) test = pcand->getChi2() < 1000 ? kTRUE: kFALSE;

    if(!test) return kFALSE;

    if(test) test = pcand->getBeta() > 0.9          ? kTRUE : kFALSE ;
    if(test) test = pcand->getMetaMatchQuality() < 3 ? kTRUE : kFALSE ;

    return test;
}

Bool_t selectLeptonsTreeFilter_old(HParticleCand* pcand){

    //  selection function for lepton candidates.
    Bool_t select = kFALSE;
    if(pcand->isFlagAND(5,
			Particle::kIsAcceptedHitRICH,
			Particle::kIsAcceptedHitInnerMDC,
			Particle::kIsAcceptedHitOuterMDC,
			Particle::kIsAcceptedHitMETA,
			Particle::kIsAcceptedRK)
       &&
       pcand->getBeta() > 0.9
       &&
       pcand->getChi2() < 1000) select = kTRUE;

    if(select) select = fabs(pcand->getDeltaTheta()) < 4    ? kTRUE : kFALSE ;
    if(select) select = fabs(pcand->getDeltaPhi())   < 4    ? kTRUE : kFALSE ;
    if(select) select = pcand->getMomentum()         < 1500 ? kTRUE : kFALSE ;


   return select;
}

Bool_t selectLeptonsTreeFilter_new(HParticleCand* pcand){
   //crg 271117
   //no leptons
   return kFALSE;
}

Bool_t selectEventTreeFilter(TObjArray* ar)
{

    HCategory* catCand = gHades->getCurrentEvent()->getCategory(catParticleCand);
    if(catCand)
    {
        Int_t n=catCand->getEntries();
        HParticleCand* cand=0;
        Int_t nProton = 0;

	    for(Int_t i = 0; i < n; ++i)
        {
            cand = HCategoryManager::getObject(cand,catCand,i);
            if(cand)
                if(cand->getPID() == 14) 
                    ++nProton;
	    }
	    if(nProton < 1) 
            return kFALSE;  // check if potential pions are there

    } 
    else 
        return kFALSE;

    return kTRUE;
}

void addFilter(HTaskSet *masterTaskSet, TString inFile,TString outDir)
{

    TString dir    = gSystem->DirName (inFile);
    //TString file   = HTime::stripFileName(inFile);
    TString file   = gSystem->BaseName (inFile);

    TString fileWoPath = file;
    //fileWoPath.ReplaceAll(".hld",""); // real data
    fileWoPath.ReplaceAll(".root",""); // sim data
    Int_t evtBuild = HTime::getEvtBuilderFileName(fileWoPath,kFALSE);

    HParticleTree* parttree = new HParticleTree("proton_tree","proton_tree");

//       SUPPORTED CATEGORIES:
//
//       catParticleCand
//       catParticleEvtInfo (fullcopy)
//       catStart2Hit       (fullcopy)
//       catStart2Cal       (fullcopy)
//       catTBoxChan        (fullcopy)
//       catWallHit         (fullcopy)
//       catPionTrackerRaw  (fullcopy)
//       catPionTrackerCal  (fullcopy)
//       catPionTrackerHit  (fullcopy)
//       catPionTrackerTrack(fullcopy)
//       catTofHit
//       catTofCluster
//       catRpcCluster
//       catShowerHit
//       catRichHit
//       catRichDirClus     (fullcopy
//       catRichCal         (fullcopy
//       catMdcSeg
//       catMdcHit
//       catMdcCal1
//       catMdcClus
//       catMdcClusInf
//       catMdcClusFit
//       catMdcWireFit
//       catGeantKine             (fullcopy)
//       catMdcGeantRaw           (fullcopy)
//       catTofGeantRaw           (fullcopy)
//       catRpcGeantRaw           (fullcopy)
//       catShowerGeantRaw        (fullcopy)
//       catWallGeantRaw          (fullcopy)
//       catRichGeantRaw (+1,+2)  (fullcopy)


    //------Long-List-Of-Categories-which-are-needed-in-this-DST-production------//
    //----------------------------------------------------------------------------//
    Cat_t pPersistentCatAll[] =
    {
	//catRichDirClus,catRichHit, catRichCal,  // full copy
	catMdcCal1, catMdcSeg, catMdcHit,
	//catShowerHit,
	//catTofHit, catTofCluster,
	//catRpcHit, catRpcCluster,
	catParticleCand, catParticleEvtInfo,
	catWallHit,
	catStart2Cal, catStart2Hit,
	catTBoxChan
    };

    //------Short-List-Of-Categories-which-are-needed-in-this-DST-production------//
    //----------------------------------------------------------------------------//
    Cat_t pPersistentCat[] =
    {
	//catRichDirClus,catRichHit, catRichCal,  // full copy
	catMdcCal1, catMdcSeg, catMdcHit,
	//catShowerHit,
	//catTofHit, catTofCluster,
	//catRpcHit, catRpcCluster,
	catParticleCand, catParticleEvtInfo,
	catWallHit,
	catStart2Cal, catStart2Hit,
	catTBoxChan
    };

    //------lsit of full copy categories------------------------------------------//
    //----------------------------------------------------------------------------//
    Cat_t pPersistentCatFull[] =
    {
	//catRichDirClus,catRichHit,catRichCal
	//catMdcCal1, 
	catMdcSeg,
	// catMdcHit,
	//catParticleCand,
	catParticleEvtInfo,
	//catParticleCand, catParticleEvtInfo,
	catStart2Cal, catStart2Hit,
	//catWallHit
    };


    if(evtBuild == 1) {
	   parttree->setEventStructure(sizeof(pPersistentCatAll)/sizeof(Cat_t) ,pPersistentCatAll);
       parttree->setEventStructure(sizeof(pPersistentCatFull)/sizeof(Cat_t),pPersistentCatFull,kTRUE);     // add more
    } else {
	parttree->setEventStructure(sizeof(pPersistentCat)/sizeof(Cat_t)    ,pPersistentCat);
        //parttree->setEventStructure(sizeof(pPersistentCatFull)/sizeof(Cat_t),pPersistentCatFull,kTRUE);
    }
    
    parttree->setSkipEmptyEvents(kTRUE);
    parttree->setUserSelectionHadrons(selectProtonsTreeFilter);
    //parttree->setUserSelectionHadrons(selectHadronsTreeFilter);
    parttree->setUserSelectionLeptons(selectLeptonsTreeFilter_new);
    
    //crg011217 TString outfile = Form("%s/%s_filtered.root",outDir.Data(),file.Data());
    
    TString rootSuffix =".root";
    //TString outFile   = gSystem->BaseName(outDir.Data());
    if(outDir.EndsWith(rootSuffix)) outDir.ReplaceAll(rootSuffix,"");
    
    TString outfile = outDir + "_filtered" + rootSuffix; //Form("%s/%s_filtered.root",outDir.Data(),file.Data());
    outfile.ReplaceAll("//","/");
    parttree->setUserSelectionEvent(selectEventTreeFilter,NULL);
    parttree->setOutputFile(outfile,"Filter","RECREATE",2 );
    masterTaskSet->add(parttree);
}












