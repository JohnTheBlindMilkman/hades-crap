#include "hades.h"
#include "hloop.h"
#include "hcategory.h"
#include "hcategorymanager.h"
#include "hparticlecand.h"
#include "hparticledef.h"

#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"

//#include "../newFemtoAnalysis.cc"
#include "../newQaAnalysis.cc"
#include <iostream>

int main(int argc, char **argv)
{
    TROOT Analysis("Analysis","compiled analysis macro");

    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1
    std::cout << argc << " arguments " << std::endl;
    if(argc>1) std::cout << "arg1 =" << argv[1] << std::endl;
    if(argc>2) std::cout << "arg2 =" << argv[2] << std::endl;
    if(argc>3) std::cout << "arg3 =" << argv[3] << std::endl;

    TString number;
    TString nevts;
    switch (argc)
    {
        case 4:       // just inputfile name + nEvents
            nevts  = argv[3];
            return newQaAnalysis(TString(argv[1]),TString(argv[2]),nevts.Atoi());

        default:
            cerr<<"ERROR: analysis() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="", nevents=1000"<<endl;
            return 1; // fail
    }
}
