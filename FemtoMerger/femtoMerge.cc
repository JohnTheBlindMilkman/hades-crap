#include "indicators.hpp"
#include "HelperFunctions.hxx"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"

int main(int argc, char **argv)
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
        return 1;
    }
    
    TFile *inpFile = TFile::Open(argv[1]);
    if (!inpFile->IsOpen() || inpFile->IsZombie())
    {
        helper.PrintHelp();
        return 1;
    }

    const TString signName(argv[2]);
    const TString bckgName(argv[3]);
    const int ktMax = atoi(argv[4]);
    const int yMax = atoi(argv[5]);
    const int psiMax = atoi(argv[6]);

    if (ktMax < 2 || yMax < 2 || psiMax < 2)
    {
        helper.PrintHelp();
        return 1;
    }

    const int maxHistos = ktMax * yMax * psiMax * 2;

    indicators::ProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::PrefixText{"Reading histograms "},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{maxHistos}};

    indicators::ProgressSpinner spinner{
                    indicators::option::PostfixText{"Merging results"},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::SpinnerStates{std::vector<std::string>{"⠈", "⠐", "⠠", "⢀", "⡀", "⠄", "⠂", "⠁"}},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};

    indicators::show_console_cursor(false);

    if (signName.Contains("inv") && bckgName.Contains("inv"))
    {
        std::vector<std::vector<std::vector<TH1D*> > > 
        hSign(ktMax,std::vector<std::vector<TH1D*> >(yMax,std::vector<TH1D*>(psiMax,nullptr))), 
        hBckg(ktMax,std::vector<std::vector<TH1D*> >(yMax,std::vector<TH1D*>(psiMax,nullptr)));

        auto job1 = [&bar]()
        {
            while (true) 
            {
                if (bar.is_completed()) 
                {
                    bar.set_option(indicators::option::ForegroundColor{indicators::Color::green});
                    bar.set_option(indicators::option::PrefixText{"✔"});
                    /*bar.set_option(indicators::option::HideBarWhenComplete{true});*/
                    bar.set_option(indicators::option::ShowPercentage{false});
                    bar.set_option(indicators::option::PostfixText{"Finished!"});
                    bar.mark_as_completed();
                    break;
                }
            }
        };

        std::thread thread1(job1);

        for (int kt = 1; kt <= ktMax; ++kt)
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    hSign[kt-1][y-1][psi-1] = inpFile->Get<TH1D>(TString::Format("%s%d%d%d",signName.Data(),kt,y,psi));
                    if (hSign[kt-1][y-1][psi-1] != nullptr)
                        hSign[kt-1][y-1][psi-1]->Sumw2();
                    bar.tick();

                    hBckg[kt-1][y-1][psi-1] = inpFile->Get<TH1D>(TString::Format("%s%d%d%d",bckgName.Data(),kt,y,psi));
                    if (hBckg[kt-1][y-1][psi-1] != nullptr)
                        hBckg[kt-1][y-1][psi-1]->Sumw2();
                    bar.tick();
                }

        bar.mark_as_completed();

        auto job2 = [&spinner]()
        {
            while (true)
            {
                if (spinner.is_completed())
                {
                    spinner.set_option(indicators::option::PostfixText{"Authenticated!"});
                    spinner.mark_as_completed();
                    break;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(40));
            }  
        };

        std::thread thread2(job2);

        std::this_thread::sleep_for(std::chrono::seconds(5));
	spinner.mark_as_completed();

	std::cout << "Waiting for helper threads to finish...\n";
        thread1.join();
        thread2.join();
    }
    else if (signName.Contains("osl") && bckgName.Contains("osl"))
    {
        std::vector<std::vector<std::vector<TH3D*> > > 
        hSign(ktMax,std::vector<std::vector<TH3D*> >(yMax,std::vector<TH3D*>(psiMax,nullptr))), 
        hBckg(ktMax,std::vector<std::vector<TH3D*> >(yMax,std::vector<TH3D*>(psiMax,nullptr)));
    }
    else
    {
        helper.PrintHelp();
        indicators::show_console_cursor(true);
        return 1;
    }
    indicators::show_console_cursor(true);
}
