#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

#include "../FemtoMixer/PairUtils.hxx"


void femtoMerge1D()
{
    const TString fileName = "/u/kjedrzej/hades-crap/slurmOutput/apr12ana_all_25_09_24.root";
    const TString SignNameBase = "hQinvSign_";
    const TString BckgNameBase = "hQinvBckg_";
    constexpr auto ktIndexing = Mixing::PairGrouping{}.GetKtIndexSequence();
    constexpr auto rapIndexing = Mixing::PairGrouping{}.GetRapIndexSequence();

    std::array<std::array<TH1D*,rapIndexing.size()>,ktIndexing.size()> hSign, hBckg;
    std::array<TH1D*,ktIndexing.size()> hSignKt, hBckgKt;
    std::array<TH1D*,rapIndexing.size()> hSignRap, hBckgRap;
    TH1D *hSignInteg = nullptr, *hBckgInteg = nullptr;

    // kt-differential
    for (auto &hist : hSignKt) 
        hist = nullptr;
    for (auto &hist : hBckgKt) 
        hist = nullptr;
    // y-differential
    for (auto &hist : hSignRap) 
        hist = nullptr;
    for (auto &hist : hBckgRap) 
        hist = nullptr;
    // kt- and y- differential
    for (auto &signKt : hSign)
        for (auto &hist : signKt)
            hist = nullptr;
    for (auto &bckgKt : hBckg)
        for (auto &hist : bckgKt)
            hist = nullptr;
    
    TFile *inpFile = TFile::Open(fileName);

    for (const std::size_t kt : ktIndexing)
        for (const std::size_t y : rapIndexing)
        {
            hSign[kt-1][y-1] = inpFile->Get<TH1D>(TString::Format("%s%02ld%02ld",SignNameBase.Data(),kt,y));
            if (hSign[kt-1][y-1] != nullptr)
            {
                hSign[kt-1][y-1]->SetName(TString::Format("hQinvSignKt%ldY%ld",kt,y));
                hSign[kt-1][y-1]->Sumw2();
            }
            hBckg[kt-1][y-1] = inpFile->Get<TH1D>(TString::Format("%s%02ld%02ld",BckgNameBase.Data(),kt,y));
            if (hBckg[kt-1][y-1] != nullptr)
            {
                hBckg[kt-1][y-1]->SetName(TString::Format("hQinvBckgKt%ldY%ld",kt,y));
                hBckg[kt-1][y-1]->Sumw2();
            }

            if (hSignInteg == nullptr && hSign[kt-1][y-1] != nullptr && hBckgInteg == nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignInteg = new TH1D(*hSign[kt-1][y-1]);
                hSignInteg->SetName("hQinvSignInteg");
                hBckgInteg = new TH1D(*hBckg[kt-1][y-1]);
                hBckgInteg->SetName("hQinvBckgInteg");
            }
            else if (hSign[kt-1][y-1] != nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignInteg->Add(hSign[kt-1][y-1]);
                hBckgInteg->Add(hBckg[kt-1][y-1]);
            }
        }

    for (const std::size_t kt : ktIndexing)
        for (const std::size_t y : rapIndexing)
        {
            if (hSignKt[kt-1] == nullptr && hSign[kt-1][y-1] != nullptr && hBckgKt[kt-1] == nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignKt[kt-1] = new TH1D(*hSign[kt-1][y-1]);
                hSignKt[kt-1]->SetName(TString::Format("hQinvSignKt%ld",kt));
                hBckgKt[kt-1] = new TH1D(*hBckg[kt-1][y-1]);
                hBckgKt[kt-1]->SetName(TString::Format("hQinvBckgKt%ld",kt));
            }
            else if (hSign[kt-1][y-1] != nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignKt[kt-1]->Add(hSign[kt-1][y-1]);
                hBckgKt[kt-1]->Add(hBckg[kt-1][y-1]);
            }
        }

    for (const std::size_t y : rapIndexing)
        for (const std::size_t kt : ktIndexing)
        {
            if (hSignRap[y-1] == nullptr && hSign[kt-1][y-1] != nullptr && hBckgRap[y-1] == nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignRap[y-1] = new TH1D(*hSign[kt-1][y-1]);
                hSignRap[y-1]->SetName(TString::Format("hQinvSignY%ld",y));
                hBckgRap[y-1] = new TH1D(*hBckg[kt-1][y-1]);
                hBckgRap[y-1]->SetName(TString::Format("hQinvBckgY%ld",y));
            }
            else if (hSign[kt-1][y-1] != nullptr && hBckg[kt-1][y-1] != nullptr)
            {
                hSignRap[y-1]->Add(hSign[kt-1][y-1]);
                hBckgRap[y-1]->Add(hBckg[kt-1][y-1]);
            }
        }

    TString otpFileName = fileName;
    TFile *otpFile = TFile::Open(otpFileName.Insert(otpFileName.First('.'),"_processed"),"recreate");

    // kt- and y- integrated
    if (hSignInteg != nullptr)
        hSignInteg->Write();
    if (hBckgInteg != nullptr)
        hBckgInteg->Write();

    // kt-differential
    for (const auto &hist : hSignKt) 
        if (hist != nullptr)
            hist->Write();
    for (const auto &hist : hBckgKt) 
        if (hist != nullptr)
            hist->Write();
    // y-differential
    for (const auto &hist : hSignRap) 
        if (hist != nullptr)
            hist->Write();
    for (const auto &hist : hBckgRap) 
        if (hist != nullptr)
            hist->Write();
    
    // kt- and y- differential
    for (const auto &signKt : hSign)
        for (const auto &hist : signKt)
            if (hist != nullptr)
                hist->Write();
    for (const auto &bckgKt : hBckg)
        for (const auto &hist : bckgKt)
            if (hist != nullptr)
                hist->Write();
}