
#ifndef __CINT__

#include <iostream>
#include <array>
#include <sstream>
#include <iomanip>

#include "TreeMacro.h"

#endif

Int_t TreeAnalysis()
{
    constexpr int firstDay{106};//95
    constexpr int lastDay{126};//126
    const std::string dirBase{"/lustre/hades/dst/apr12/gen9/"};

    TChain* Chain;
    TFile *outputfile;
    
    std::string str;
    for (int i = firstDay; i <= lastDay; ++i)
    {
        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << i;
        str = dirBase + ss.str() + "/qa/*_Tree.root";
        std::cout << "Adding flies from " << str << "\n";

        Chain = new TChain("T");
        Chain->Add(str.c_str());

        TString ofilename("");
        ofilename.Append("./output/QAHistogram_day_");
        ofilename.Append(ss.str().c_str());
        ofilename.Append(".root");

        outputfile = new TFile(ofilename,"RECREATE");

        AnalysisMacro AM(outputfile,Chain);
        AM.Loop();
        AM.finalize(outputfile);

        delete Chain;
    }
}

#ifndef __CINT__
int main(int argc, char **argv)
{
    return TreeAnalysis();
}
#endif
