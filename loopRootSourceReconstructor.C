#include "hades.h"
#include "hrootsource.h"
#include "hevent.h"
#include "htaskset.h"


#include "FemtoAnalysis.hxx"


#include <iostream>
using namespace std;


int loopRootSourceReconstructor(Int_t nEvents    = 100,
	    TString inputDir = "/misc/hadesprojects/analysis_workshop/2022/inputfile",
	    TString inFile   = "Au_Au_1230MeV_1000evts_1_1_dst_apr12.root")
{
    Hades* myHades = new Hades;
    gHades->setQuietMode(2);

    HRootSource* source = new HRootSource();
    source->setDirectory(((Text_t *)inputDir.Data()));
    source->addFile((Text_t *)inFile.Data());

    myHades->setDataSource(source);


    HTaskSet* tasks = gHades->getTaskSet("all");
    tasks->add(new FemtoAnalysis("MyReco","MyReco","myreco.root"));

    if(!myHades->init()){
	cout<<"Hades Init() failed!"<<endl;
	exit(1);
    }

    gHades->eventLoop(nEvents);

    delete myHades;

    return 0;
}
