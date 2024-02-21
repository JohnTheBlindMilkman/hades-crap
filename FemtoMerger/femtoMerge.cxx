#include "indicators.hpp"
#include "HelperFunctions.hxx"
#include "TFile.h"

int main(int argc, char *argv[])
{
    // verify if input is valid
    // create the main object
    // perform first job
    // perform second job
    // finish

    JJFemtoMerger::HelperFunctions helper;

    if (argc < 7)
    {
        helper.PrintHelp();
        return 0;
    }
    
    TFile *inpFile = TFile::Open(argv[1]);
    if (!inpFile->IsOpen() || inpFile->IsZombie())
    {
        helper.PrintHelp();
        return 0;
    }

    TString signName(argv[2]);
    TString bckgName(argv[3]);
    int ktMax = atoi(argv[4]);
    int yMax = atoi(argv[5]);
    int psiMax = atoi(argv[6]);

    if (ktMax < 2 || yMax < 2 || psiMax < 2)
    {
        helper.PrintHelp();
        return 0;
    }

    indicators::ProgressBar bar{indicators::option::BarWidth{50},
                   indicators::option::ForegroundColor{indicators::Color::yellow},
                   indicators::option::ShowElapsedTime{true},
                   indicators::option::ShowRemainingTime{true},
                   indicators::option::PrefixText{"Reading histograms "},
                   indicators::option::FontStyles{
                        std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{40}};

    indicators::ProgressSpinner spinner{
                    indicators::option::PostfixText{"Merging results"},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::SpinnerStates{
                        std::vector<std::string>{"⠈", "⠐", "⠠", "⢀", "⡀", "⠄", "⠂", "⠁"}},
                    indicators::option::FontStyles{
                        std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
}