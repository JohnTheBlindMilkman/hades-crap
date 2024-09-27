#include <iostream>
#include <fstream>
#include <regex>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "THashList.h"

#include "../Externals/indicators.hpp"

constexpr std::size_t nSectors{6}, nPlanes{4}, nChannels{12};

// find bin in histogram which has the label corresponidng to the currenlty processed hld filepath
int GetBin(TH1F *hist, std::string name)
{
    const int nbins =  hist->GetNbinsX();
    for (int bin = 1; bin <= nbins; ++bin)
        if (name == std::string(hist->GetXaxis()->GetLabels()->At(bin-1)->GetName()))
            return bin;

    return 0;
}

// check if at any histogram the voltage is less than what we allow (technicaly this also checks if the voltage is over nominal)
bool IsVoltageBelow(std::array<std::array<std::array<TH1F*,nChannels>,nPlanes>,nSectors> histArr, int day, int bin, float diff)
{
    static const std::map<int,std::array<std::array<float,nPlanes>,nSectors>> nominalVoltages
    {
        {95,{{{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700}}}},
        {96,{{{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700},{1800,1375,1500,1700}}}},
        {97,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1725,1375,1500,1700},{1750,1325,1500,1700},{1750,1375,1500,1700}}}},
        {98,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1300,1500,1700},{1750,1375,1500,1700}}}},
        {99,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {100,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {101,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {102,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1200,1500,1700},{1750,800,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {103,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1200,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {104,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1050,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {105,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,1000,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {106,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {107,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {108,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {109,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {110,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {111,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {112,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {113,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {114,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {115,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1375,1500,1700},{1750,1375,1500,1700}}}},
        {116,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {117,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {118,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {119,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {120,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {121,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {122,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {123,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {124,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {125,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,500,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}},
        {126,{{{1750,1375,1500,1700},{1750,1375,1500,1700},{1750,890,1500,1700},{1750,1200,1500,1700},{1750,1350,1500,1700},{1750,1375,1500,1700}}}}
    };

    for (std::size_t sec = 0; sec < nSectors; ++sec)
        for (std::size_t pl = 0; pl < nPlanes; ++pl)
        {
            if (std::abs(nominalVoltages.at(day)[sec][pl] - histArr[sec][pl][0]->GetBinContent(bin)) > diff && 
                histArr[sec][pl][0]->GetBinContent(bin) > 20)
            {
                return true;
            }
        }

    return false;
}

void markBadRuns()
{
    constexpr float voltageDiff{5.0};
    constexpr int firstDay{95},lastDay{126}, maxCount{10};
    const std::string qaFileBase{"/u/kjedrzej/hades-crap/QaTtreeAnalysis/output/"};
    const std::string listFileBase{"/u/kjedrzej/hades-crap/loopDST/listsApr12/"};
    const std::string outputFileBase{"/u/kjedrzej/hades-crap/loopDST/listsApr12good/"};
    const std::regex rgxp{"be12[0-9]{11}"}; // regexp for finding hld file names, e.g. be1210204502608 ("string which starts with be12 and has 11 numbers following it")

    // just a progress bar (in this case, technically, just a fancy way of showing which line we currenlty read, not the actual progress)
    indicators::IndeterminateProgressBar spinner{
        indicators::option::BarWidth{20},
        indicators::option::Start{"["},
        indicators::option::Fill{"·"},
        indicators::option::Lead{"<>"},
        indicators::option::End{" ]"},
        indicators::option::PostfixText{"Reading line: 0"},
        indicators::option::ForegroundColor{indicators::Color::yellow},
        indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};

    std::array<std::array<std::array<TH1F*,nChannels>,nPlanes>,nSectors> histVec;
    std::size_t counter = 0;
    std::size_t progressCounter = maxCount;

    std::ifstream istream;
    std::ofstream ostream;
    std::string tmp,strIn,strOut;
    std::smatch nameMatch;

    TFile *inpRoot;
    TFile *otpRoot;

    indicators::show_console_cursor(false);

    for (int i = firstDay; i <= lastDay; ++i)
    {
        // create day-of-the-year string, e.g. 094, 121, etc.
        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << i;
        strIn = qaFileBase + "QAHistogram_day_" + ss.str() + ".root";
        strOut = qaFileBase + "MDC_min_combined_day_" + ss.str() + ".root";

        // load all histograms of interest into memory
        inpRoot = TFile::Open(strIn.c_str());
        for (std::size_t sec = 0; sec < nSectors; ++sec)
            for (std::size_t pl = 0; pl < nPlanes; ++pl)
                for (std::size_t lay = 0; lay < nChannels; ++lay)
                {
                    histVec[sec][pl][lay] = inpRoot->Get<TH1F>(TString::Format("MDC/hHVMinMdc[%ld][%ld][%ld]",sec,pl,lay));
                }

        istream.open(listFileBase + "day_" + ss.str() + ".list");
        ostream.open(outputFileBase + "day_" + ss.str() + ".list");

        // for each filepath: 
        //      1. check if matches regexp, 
        //      2. find which bin coresponds to the regexp (assuming all histograms have the same binning),
        //      3. check if at any layer the min HV was below nominal HV by a certain value,
        //      4. if all is ok, save filepath to new file list (we have now no "bad runs").
        while(istream >> tmp)
        {
            ++counter;
            if (--progressCounter == 0)
            {
                spinner.set_option(indicators::option::PostfixText{"Reading line: " + std::to_string(counter)});
                spinner.tick();
                progressCounter = maxCount;
            }
            
            if (std::regex_search(tmp,nameMatch,rgxp))
            {
                int runBin = GetBin(histVec[0][0][0],nameMatch[0]);
                if (runBin > 0)
                {
                    if (!IsVoltageBelow(histVec,i,runBin,voltageDiff))
                        ostream << tmp << "\n";
                }
            }
        }
        ostream.close();
        istream.close();

        /* otpRoot = TFile::Open(strOut.c_str(),"recreate");
        std::array<TCanvas*,4> cSectors;
        for (std::size_t pl = 0; pl < nPlanes; ++pl)
        {
            cSectors[pl] = new TCanvas(TString::Format("day_%d_c_pl%ld",i,pl),"",1600,900);
            for (std::size_t sec = 0; sec < nSectors; ++sec)
            {
                if (sec == 0)
                    histVec[sec][pl][0]->Draw("p pmc");
                else
                    histVec[sec][pl][0]->Draw("p pmc same");
            }
            cSectors[pl]->BuildLegend();
            cSectors[pl]->Write();
        }
        otpRoot->Close(); */
    }

    spinner.mark_as_completed();
    indicators::show_console_cursor(true);
}