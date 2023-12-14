#include "FemtoAnalysis.hxx"

ClassImp(FemtoAnalysis)

FemtoAnalysis::FemtoAnalysis() : hello("Hello World!")
{
}

FemtoAnalysis::FemtoAnalysis(const Text_t *name = "",const Text_t *title ="",TString outfile="myreco.root") : HReconstructor(name,title), hello("Hello World!")
{

}

FemtoAnalysis::~FemtoAnalysis()
{
}

Bool_t FemtoAnalysis::init()
{
    Greet();
    return kTRUE;
}

Bool_t FemtoAnalysis::reinit()
{
    return kTRUE;
}

Int_t FemtoAnalysis::execute()
{
    return 0;
}

Bool_t FemtoAnalysis::finalize()
{
    return kTRUE;
}

void FemtoAnalysis::Greet()
{
    std::cout << hello << std::endl;
}