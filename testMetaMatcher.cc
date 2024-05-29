#include "hades.h"
#include "hloop.h"
#include "hdst.h"
#include "hcategory.h"
#include "hcategorymanager.h"
#include "htaskset.h"
#include "hgeomvector.h"
//
//
#include "hparticlemetamatcher.h"
//
#include "hparticlecand.h"
#include "hparticleevtinfo.h"
#include "htofcluster.h"
#include "hrpccluster.h"
#include "hemccluster.h"
#include "emcdef.h"
//
#include "hparticledef.h"
//
#include "TString.h"
//
#include <iostream>
#include <iomanip>
//
using namespace std;
//
Int_t testMetaMatcher(TString infileList="/lustre/hades/dst/apr12/gen9/120/root/be1212000004001.hld_dst_apr12.root",TString outfile="test.root",Int_t nEvents=50)
{
    HLoop loop(kTRUE);
    //-------------------------------------------------
    if(1)
    {
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
        TString paramSource      = "root"; 
        TString paramrelease     = "APR12_dst_gen9"; 
        HDst::setupSpectrometer(beamtime,mdcMods,"rich,mdc,tof,rpc,shower,wall,start,tbox");
        HDst::setupParameterSources(paramSource,asciiParFile,rootParFile,paramrelease);
    }
    //-------------------------------------------------
    Bool_t ret =kFALSE;
    if(infileList.Contains(","))
    {
	    ret = loop.addMultFiles(infileList); 
    } 
    else
    {
	    ret = loop.addFiles(infileList);
    }

    if(ret == 0) 
    {
        cout<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
        return 1;
    }

    if(!loop.setInput("")) 
    {
        cout<<"READBACK: ERROR : cannot read input !"<<endl;
        exit(1);
    }
    loop.printCategories();
    loop.printChain();

    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);

    //-------------------------------------------------
    HTaskSet *masterTaskSet = gHades->getTaskSet("all");
    HParticleMetaMatcher* matcher = new HParticleMetaMatcher();
    matcher->setDebug();
    matcher->setUseEMC(kFALSE);
    masterTaskSet->add(matcher);
    //-------------------------------------------------

    Int_t entries = loop.getEntries();
    if(nEvents < entries && nEvents >= 0 ) 
        entries = nEvents;

    for (Int_t i = 0; i < entries; i++)
    {
        Int_t nbytes =  loop.nextEvent(i);             
        if(nbytes <= 0) 
            { cout<<nbytes<<endl; break; }

        if(i%1000 == 0) 
            cout<<"event "<<i<<endl;

        HParticleEvtInfo* evtInfo=0;
        evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );
        if(evtInfo&&!evtInfo->isGoodEvent(Particle::kGoodTRIGGER|          
                        Particle::kGoodVertexClust|
                        Particle::kGoodVertexCand|
                        Particle::kGoodSTART|
                        Particle::kNoPileUpSTART|
                        Particle::kNoVETO|
                        Particle::kGoodSTARTVETO|
                        Particle::kGoodSTARTMETA
                        )) continue;
        if(candCat)
        {
            Int_t size = candCat->getEntries();
            HParticleCand* cand=0;
            for(Int_t j = 0; j < size; j++) 
            {
                cand = HCategoryManager::getObject(cand,candCat,j);
                if(cand) 
                {
                    if(!cand->isFlagBit(kIsUsed)) 
                        { continue; } 
                    if(0&cand->isRpcClstUsed())
                    {
                        const HRpcCluster* rpc = matcher->recalcRpc(cand);
                        Int_t s,col0,col1,cell0,cell1;
                        HGeomVector hit0,hit1;
                        matcher->predictRpcCell(cand,hit0,hit1,s,col0,cell0,col1,cell1); 
                    }
                    if(0&&(cand->isTofClstUsed() || cand->isTofHitUsed()) )
                    {
                        const HTofCluster* tof = matcher->recalcTof(cand);
                        Int_t s,mod0,mod1,cell0,cell1;
                        HGeomVector hit0,hit1;
                        matcher->predictTofCell(cand,hit0,hit1,s,mod0,cell0,mod1,cell1); 
                    }
                    if(0&&cand->getEmcInd()>-1)
                    {
                        const HEmcCluster* emc = matcher->recalcEmc(cand);
                        Int_t s,pos,cell;
                        HGeomVector hit;
                        vector<Int_t> vcells;
                        matcher->predictEmcCell(cand,hit,s,pos,cell); 
                    }
                    //------------------------------------------------------------------------
                    HParticleWireManager& wM = matcher->getWireManager();
                    wM.setDebug();
                    wM.setWireRange(2);
                    //Int_t nw     = wM.getNWires    ()           { return nWires;}     number of wires in event
                    //Int_t nwused = wM.getNWiresUsed()           { return nWiresUsed;}
                    //Int_t n      = wM.isUsedNtimes(Int_t s,Int_t m,Int_t l,Int_t c);  how many time this cell has been used ? size 6:4:6:220
                    //std::cout << "Current wire info size: " <<  wM.getWireInfo().size() << "\n";
                    //Bool_t wM.getWireInfo(UInt_t i,HParticleWireInfo& w,HParticleCand* c = 0);  index of HParticleCand -> WireInfo
                    //Bool_t wM.sharedWires(HParticleWireInfo& w1     ,HParticleWireInfo& w2     ,Int_t io=0,Int_t range = -1); compare 2 WireInfo objects if they shared wires
                    //Bool_t wM.sharedWires(HParticleCand* c1,HParticleCand* c2,Int_t io=0,Int_t range = -1); compare 2 HParticleCand objects if they shared wires
                    if(1)
                    {
                        HGeomVector hitmdc;
                        HGeomVector hitmdcsec;
                        HGeomVector hitmdclab;
                        Int_t module = 0;
                        matcher->predictMdcHit   (cand,hitmdc   ,module); 
                        matcher->predictMdcHitSec(cand,hitmdcsec,module); 
                        matcher->predictMdcHitLab(cand,hitmdclab,module);
                    }
                    //------------------------------------------------------------------------
                }
            }
        }
    }
    //-------------------------------------------------
    if(gHades)
        gHades->finalizeTasks();
    //-------------------------------------------------
    delete gHades;
    return 0;
}