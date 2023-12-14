#ifndef FemtoAnalysis_hxx
    #define FemtoAnalysis_hxx

    // base functionality
    #include "hreconstructor.h"
    #include "hcategory.h"
    #include "hcategorymanager.h"

    // data objects
    #include "hparticledef.h"
    #include "hparticleevtinfo.h"

    // ROOT specific
    //#include "TString.h"

    // C++ specific
    #include <iostream>

    class FemtoAnalysis : public HReconstructor
    {
    private:
        void Greet();

        TString hello;

    public:
        FemtoAnalysis();
        FemtoAnalysis(const Text_t *name = "",const Text_t *title ="",TString outfile="myreco.root");
        ~FemtoAnalysis();

        Bool_t init();
        Bool_t reinit();
        Int_t execute();
        Bool_t finalize();

        ClassDef(FemtoAnalysis,0);
    };

#endif