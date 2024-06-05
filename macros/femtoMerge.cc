#include "../FemtoMerger/indicators.hpp"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"

void femtoMerge(TString fileName = "/u/kjedrzej/hades-crap/slurmOutput/apr12ana_quarter_24_06_05.root", TString signName = "hQoslSign_", TString bckgName = "hQoslBckg_", int ktMax = 5, int yMax = 3, int psiMax = 8)
{ 
    const int maxHistos = ktMax * yMax * psiMax * 2;

    TFile *inpFile = TFile::Open(fileName);

    indicators::BlockProgressBar bar{indicators::option::BarWidth{80},
                    indicators::option::ForegroundColor{indicators::Color::yellow},
                    indicators::option::ShowPercentage{true},
                    indicators::option::PrefixText{"Reading histograms "},
                    indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
                    indicators::option::MaxProgress{maxHistos}};

    indicators::ProgressSpinner spinner{
                    indicators::option::PostfixText{"Merging results"},
                    indicators::option::ShowPercentage{false},
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
                    std::this_thread::sleep_for(std::chrono::milliseconds(40)); // adding small delay because it is too fast
                    std::cout << termcolor::bold << termcolor::green << "Finished!\n" << termcolor::reset;
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

        if(!bar.is_completed())
            bar.mark_as_completed();
        thread1.join();

        auto job2 = [&spinner]()
        {
            while (true)
            {
                if (spinner.is_completed())
                {
                    std::cout << termcolor::bold << termcolor::green << "Merging complete!\n" << termcolor::reset;
                    break;
                }
                else
                    spinner.tick();
                std::this_thread::sleep_for(std::chrono::milliseconds(40));
            }  
        };

        std::thread thread2(job2);

        std::vector<TH1D*> hSignKt(ktMax,nullptr), hBckgKt(ktMax,nullptr);
        for (int kt = 1; kt <= ktMax; ++kt)
        {
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignKt[kt-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgKt[kt-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignKt[kt-1] = new TH1D(*hSign[kt-1][y-1][psi-1]);
                        hSignKt[kt-1]->SetName(TString::Format("hQinvSignKt%d",kt));
                        hBckgKt[kt-1] = new TH1D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgKt[kt-1]->SetName(TString::Format("hQinvBckgKt%d",kt));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignKt[kt-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgKt[kt-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        std::vector<TH1D*> hSignY(yMax,nullptr), hBckgY(yMax,nullptr);
        for (int y = 1; y <= yMax; ++y)
        {
            for (int kt = 1; kt <= ktMax; ++kt)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignY[y-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgY[y-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignY[y-1] = new TH1D(*hSign[kt-1][y-1][psi-1]);
                        hSignY[y-1]->SetName(TString::Format("hQinvSignY%d",y));
                        hBckgY[y-1] = new TH1D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgY[y-1]->SetName(TString::Format("hQinvBckgY%d",y));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignY[y-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgY[y-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        std::vector<TH1D*> hSignPsi(psiMax,nullptr), hBckgPsi(psiMax,nullptr);
        for (int psi = 1; psi <= psiMax; ++psi)
        {
            for (int kt = 1; kt <= ktMax; ++kt)
                for (int y = 1; y <= yMax; ++y)
                {
                    if (hSignPsi[psi-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgPsi[psi-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignPsi[psi-1] = new TH1D(*hSign[kt-1][y-1][psi-1]);
                        hSignPsi[psi-1]->SetName(TString::Format("hQinvSignPsi%d",psi));
                        hBckgPsi[psi-1] = new TH1D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgPsi[psi-1]->SetName(TString::Format("hQinvBckgPsi%d",psi));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignPsi[psi-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgPsi[psi-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        TH1D *hSignInteg = nullptr, *hBckgInteg = nullptr;
        for (int kt = 1; kt <= ktMax; ++kt)
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignInteg == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgInteg == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignInteg = new TH1D(*hSign[kt-1][y-1][psi-1]);
                        hSignInteg->SetName("hQinvSignInteg");
                        hBckgInteg = new TH1D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgInteg->SetName("hQinvBckgInteg");
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignInteg->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgInteg->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }

        TFile *otpFile = TFile::Open(fileName.Insert(fileName.First('.'),"_processed"),"recreate");
        for (int kt = 1; kt <= ktMax; ++kt)
        {
            hSignKt[kt-1]->Write();
            hBckgKt[kt-1]->Write();
        }
        for (int y = 1; y <= yMax; ++y)
        {
            hSignY[y-1]->Write();
            hBckgY[y-1]->Write();
        }
        for (int psi = 1; psi <= psiMax; ++psi)
        {
            hSignPsi[psi-1]->Write();
            hBckgPsi[psi-1]->Write();
        }
        hSignInteg->Write();
        hBckgInteg->Write();

        otpFile->Close();
        inpFile->Close();

        spinner.mark_as_completed();
        thread2.join();
    }
    else if (signName.Contains("osl") && bckgName.Contains("osl"))
    {
        std::vector<std::vector<std::vector<TH3D*> > > 
        hSign(ktMax,std::vector<std::vector<TH3D*> >(yMax,std::vector<TH3D*>(psiMax,nullptr))), 
        hBckg(ktMax,std::vector<std::vector<TH3D*> >(yMax,std::vector<TH3D*>(psiMax,nullptr)));

        auto job1 = [&bar]()
        {
            while (true) 
            {
                if (bar.is_completed()) 
                {
                    std::this_thread::sleep_for(std::chrono::milliseconds(40)); // adding small delay because it is too fast
                    std::cout << termcolor::bold << termcolor::green << "Finished!\n" << termcolor::reset;
                    break;
                }
            }
        };

        std::thread thread1(job1);

        for (int kt = 1; kt <= ktMax; ++kt)
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    hSign[kt-1][y-1][psi-1] = inpFile->Get<TH3D>(TString::Format("%s%d%d%d",signName.Data(),kt,y,psi));
                    if (hSign[kt-1][y-1][psi-1] != nullptr)
                        hSign[kt-1][y-1][psi-1]->Sumw2();
                    bar.tick();

                    hBckg[kt-1][y-1][psi-1] = inpFile->Get<TH3D>(TString::Format("%s%d%d%d",bckgName.Data(),kt,y,psi));
                    if (hBckg[kt-1][y-1][psi-1] != nullptr)
                        hBckg[kt-1][y-1][psi-1]->Sumw2();
                    bar.tick();
                }

        if(!bar.is_completed())
            bar.mark_as_completed();
        thread1.join();

        auto job2 = [&spinner]()
        {
            while (true)
            {
                if (spinner.is_completed())
                {
                    std::cout << termcolor::bold << termcolor::green << "Merging complete!\n" << termcolor::reset;
                    break;
                }
                else
                    spinner.tick();
                std::this_thread::sleep_for(std::chrono::milliseconds(40));
            }  
        };

        std::thread thread2(job2);

        std::vector<TH3D*> hSignKt(ktMax,nullptr), hBckgKt(ktMax,nullptr);
        for (int kt = 1; kt <= ktMax; ++kt)
        {
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignKt[kt-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgKt[kt-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignKt[kt-1] = new TH3D(*hSign[kt-1][y-1][psi-1]);
                        hSignKt[kt-1]->SetName(TString::Format("hQoslSignKt%d",kt));
                        hBckgKt[kt-1] = new TH3D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgKt[kt-1]->SetName(TString::Format("hQoslBckgKt%d",kt));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignKt[kt-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgKt[kt-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        std::vector<TH3D*> hSignY(yMax,nullptr), hBckgY(yMax,nullptr);
        for (int y = 1; y <= yMax; ++y)
        {
            for (int kt = 1; kt <= ktMax; ++kt)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignY[y-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgY[y-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignY[y-1] = new TH3D(*hSign[kt-1][y-1][psi-1]);
                        hSignY[y-1]->SetName(TString::Format("hQoslSignY%d",y));
                        hBckgY[y-1] = new TH3D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgY[y-1]->SetName(TString::Format("hQoslBckgY%d",y));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignY[y-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgY[y-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        std::vector<TH3D*> hSignPsi(psiMax,nullptr), hBckgPsi(psiMax,nullptr);
        for (int psi = 1; psi <= psiMax; ++psi)
        {
            for (int kt = 1; kt <= ktMax; ++kt)
                for (int y = 1; y <= yMax; ++y)
                {
                    if (hSignPsi[psi-1] == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgPsi[psi-1] == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignPsi[psi-1] = new TH3D(*hSign[kt-1][y-1][psi-1]);
                        hSignPsi[psi-1]->SetName(TString::Format("hQoslSignPsi%d",psi));
                        hBckgPsi[psi-1] = new TH3D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgPsi[psi-1]->SetName(TString::Format("hQoslBckgPsi%d",psi));
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignPsi[psi-1]->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgPsi[psi-1]->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }
        }

        TH3D *hSignInteg = nullptr, *hBckgInteg = nullptr;
        for (int kt = 1; kt <= ktMax; ++kt)
            for (int y = 1; y <= yMax; ++y)
                for (int psi = 1; psi <= psiMax; ++psi)
                {
                    if (hSignInteg == nullptr && hSign[kt-1][y-1][psi-1] != nullptr && hBckgInteg == nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignInteg = new TH3D(*hSign[kt-1][y-1][psi-1]);
                        hSignInteg->SetName("hQoslSignInteg");
                        hBckgInteg = new TH3D(*hBckg[kt-1][y-1][psi-1]);
                        hBckgInteg->SetName("hQoslBckgInteg");
                    }
                    else if(hSign[kt-1][y-1][psi-1] != nullptr && hBckg[kt-1][y-1][psi-1] != nullptr)
                    {
                        hSignInteg->Add(hSign[kt-1][y-1][psi-1]);
                        hBckgInteg->Add(hBckg[kt-1][y-1][psi-1]);
                    }
                }

        TFile *otpFile = TFile::Open(fileName.Insert(fileName.First('.'),"_processed"),"recreate");
        for (int kt = 1; kt <= ktMax; ++kt)
        {
            hSignKt[kt-1]->Write();
            hBckgKt[kt-1]->Write();
        }
        for (int y = 1; y <= yMax; ++y)
        {
            hSignY[y-1]->Write();
            hBckgY[y-1]->Write();
        }
        for (int psi = 1; psi <= psiMax; ++psi)
        {
            hSignPsi[psi-1]->Write();
            hBckgPsi[psi-1]->Write();
        }
        hSignInteg->Write();
        hBckgInteg->Write();

        otpFile->Close();
        inpFile->Close();

        spinner.mark_as_completed();
        thread2.join();
    }

    indicators::show_console_cursor(true);
}
