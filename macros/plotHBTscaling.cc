#include <array>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"

#include "../Externals/Palettes.hxx"

void MakeNice(TH1D* hist)
{
    hist->GetXaxis()->SetTitleOffset(); // invoking this functione becasue the side direction title got wonky
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetXaxis()->SetNdivisions(506);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetNdivisions(506);
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(JJColor::fWutAllColors[1]); // navy WUT
    hist->SetLineColor(JJColor::fWutSecondaryColors[1]); // blue WUT
}

void plotHBTscaling()
{
    gStyle->SetOptStat(0);

    // setting up constants
    constexpr float beamRapidity{0.74};
    constexpr std::size_t centralitySize{1};
    constexpr unsigned long ktArrSize{9}; 
    float ktVal[] = {200,400,600,800,1000,1200,1400,1600,1800,2000};
    constexpr std::array<const char *,centralitySize> centTitle
    {
        "0-10%",
        /* "10-20%",
        "20-30%",
        "30-40%" */
    };
    constexpr std::array<std::array<float,ktArrSize>,centralitySize> ktRinvVal
    {{
        {3.53619,3.49615,3.42098,3.27191,3.15268,3.07505,2.99654,2.93566,2.93476},
        /*{2.38104,2.47855,1.92715,1.60607,1.39969,1.25969,1.33657},
        {2.33543,2.13086,1.6883,1.46947,1.28822,1.25409,1.25248},
        {2.20258,1.90858,1.56513,1.46191,1.34012,1.29267,1.30359} ,
        {1.82076,1.5506,1.30417,1.38577,1.34865,1.25208,1.12921},
        {2.76133,1.37655,0.8000,0.8000,0.824145,0.904802,0.8} */
    }};
    constexpr std::array<std::array<float,ktArrSize>,centralitySize> ktRinvErr
    {{
        {0.096581,0.037686,0.031085,0.032110,0.036641,0.043290,0.050734,0.086118,0.131065},
        /*{0.321332,0.073244,0.041912,0.033537,0.034145,0.042166,0.054002},
        {0.309059,0.068947,0.040233,0.034590,0.037837,0.047166,0.066971},
        {0.367966,0.082206,0.048236,0.043702,0.051008,0.067307,0.101443} ,
        {0.459693,0.111709,0.066131,0.061287,0.073851,0.106410,0.177215},
        {0.890248,0.282144,0.020086,0.059329,1.077467,0.431814,0.361839} */
    }};
    // constexpr std::array<std::array<float,ktArrSize>,centralitySize> ktLambdaVal
    // {{
    //     {0.369406,0.283478,0.249706,0.223526,0.185309,0.198803,0.225417},
    //     {0.21512,0.268323,0.178499,0.162333,0.156201,0.150104,0.186897},
    //     {0.2693,0.226168,0.164863,0.166087,0.164345,0.185351,0.199688},
    //     {0.29207,0.21334,0.171763,0.197989,0.205483,0.223435,0.244374}/* ,
    //     {0.266259,0.171914,0.150542,0.213109,0.244485,0.241623,0.213187},
    //     {1.00000,0.209352,0.114386,0.113944,0.127891,0.153874,0.136215} */
    // }};
    // constexpr std::array<std::array<float,ktArrSize>,centralitySize> ktLambdaErr
    // {{
    //     {0.031888,0.033432,0.014318,0.010899,0.00958265,0.0113632,0.016378},
    //     {0.085105,0.024931,0.010373,0.008161,0.008419,0.010343,0.016617},
    //     {0.102378,0.020603,0.009601,0.008818,0.009961,0.014326,0.022186},
    //     {0.134009,0.023990,0.012268,0.013304,0.016687,0.024490,0.040913}/* ,
    //     {0.158196,0.027692,0.015212,0.061287,0.028433,0.041484,0.062087},
    //     {0.939921,0.092389,0.010942,0.013428,0.073720,0.116839,0.056353} */
    // }};

    constexpr unsigned long yArrSize{9}; 
    float yVal[] = {0.09 - beamRapidity,
                    0.19 - beamRapidity,
                    0.29 - beamRapidity,
                    0.39 - beamRapidity,
                    0.49 - beamRapidity,
                    0.59 - beamRapidity,
                    0.69 - beamRapidity,
                    0.79 - beamRapidity,
                    0.89 - beamRapidity,
                    0.99 - beamRapidity,
                    1.09 - beamRapidity,
                    1.19 - beamRapidity,
                    1.29 - beamRapidity,
                    1.39 - beamRapidity};
    
    constexpr std::array<std::array<float,yArrSize>,centralitySize> yRinvVal
    {{
        {3.39828,3.32912,3.28655,3.29742,3.3066,3.24906,3.19982,3.20594,3.17633},
        /*{1.07771,1.84382,2.29573,2.10163},
        {1.02425,1.73469,2.04458,1.9476},
        {0.874663,1.69724,1.88363,1.81764} ,
        {0.922745,1.58081,1.76314,1.68245},
        {0.597963,1.48342,1.36013,1.33193} */
    }};
    constexpr std::array<std::array<float,yArrSize>,centralitySize> yRinvErr
    {{
        {0.031204,0.033211,0.032388,0.041213,0.049352,0.049081,0.059002,0.071239,0.109666},
        /*{0.032909,0.032433,0.044305,0.066979},
        {0.036720,0.033131,0.043927,0.068140},
        {0.044542,0.042394,0.055940,0.081756} ,
        {0.058074,0.059165,0.082453,0.123917},
        {0.190417,0.159474,0.263952,0.313245} */
    }};
    // constexpr std::array<std::array<float,yArrSize>,centralitySize> yLambdaVal
    // {{
    //     {0.117061,0.255825,0.415177,0.341175},
    //     {0.0799886,0.206251,0.380511,0.287582},
    //     {0.0821282,0.225987,0.378103,0.312455},
    //     {0.07338,0.263113,0.398375, 0.359604}/* ,
    //     {0.0920215,0.290457,0.446639,0.395947},
    //     {0.0825105,0.382272,0.362114,0.400213} */
    // }};
    // constexpr std::array<std::array<float,yArrSize>,centralitySize> yLambdaErr
    // {{
    //     {0.004795,0.010966,0.023724,0.031424},
    //     {0.003988,0.009658,0.044305,0.027005},
    //     {0.004557,0.011042,0.023155,0.030659},
    //     {0.005111,0.016630,0.032076,0.043424}/* ,
    //     {0.008306,0.026099,0.054005,0.073745},
    //     {0.022266,0.100886,0.157562,0.209947} */
    // }};

    std::vector<TH1D*> hKtRinvArr;
    std::vector<TH1D*> hKtLambdaArr;

    for (std::size_t cent = 0; cent < centralitySize; ++cent)
    {
        TH1D *hRinvTmp = new TH1D(TString::Format("hKtRinv_%ld",cent+1),";k_{T} [MeV/c];R_{inv} [fm]",ktArrSize,ktVal);
        //TH1D *hLambdaTmp = new TH1D(TString::Format("hKtLambda_%ld",cent+1),";k_{T} [MeV/c];#lambda [a.u.]",ktArrSize,ktVal);
        for (std::size_t kt = 1; kt <= ktArrSize; ++kt)
        {
            hRinvTmp->SetBinContent(kt,ktRinvVal[cent][kt-1]);
            hRinvTmp->SetBinError(kt,ktRinvErr[cent][kt-1]);

            // hLambdaTmp->SetBinContent(kt,ktLambdaVal[cent][kt-1]);
            // hLambdaTmp->SetBinError(kt,ktLambdaErr[cent][kt-1]);
        }

        MakeNice(hRinvTmp);
        //MakeNice(hLambdaTmp);

        hKtRinvArr.push_back(hRinvTmp);
        //hKtLambdaArr.push_back(hLambdaTmp);
    }

    JJColor::CreatePrimaryWutGradient();

    TCanvas *canvKtRinv = new TCanvas("canvKtRinv","",800,450);
    canvKtRinv->SetMargin(0.15,0.05,0.15,0.05);
    for (std::size_t cent = 0; cent < centralitySize; ++cent)
    {
        hKtRinvArr[cent]->SetMarkerColor(JJColor::fWutSecondaryColors[cent+1]);
        hKtRinvArr[cent]->SetLineColor(JJColor::fWutSecondaryColors[cent+1]);
        hKtRinvArr[cent]->SetTitle(centTitle[cent]);
        if (cent == 0)
        {
            hKtRinvArr[cent]->GetYaxis()->SetRangeUser(2.5,4.5);
            hKtRinvArr[cent]->Draw("p e");
        }
        else
        {
            hKtRinvArr[cent]->Draw("p e same");
        }
    }
    //canvKtRinv->BuildLegend();

    // TCanvas *canvKtLambda = new TCanvas("canvKtLambda","",800,450);
    // canvKtLambda->SetMargin(0.2,0.02,0.15,0.02);
    // for (std::size_t cent = 0; cent < centralitySize; ++cent)
    // {
    //     hKtLambdaArr[cent]->SetMarkerColor(JJColor::fWutSecondaryColors[cent+1]);
    //     hKtLambdaArr[cent]->SetLineColor(JJColor::fWutSecondaryColors[cent+1]);
    //     hKtLambdaArr[cent]->SetTitle(centTitle[cent]);
    //     if (cent == 0)
    //     {
    //         hKtLambdaArr[cent]->GetYaxis()->SetRangeUser(0.,1.);
    //         hKtLambdaArr[cent]->Draw("p e1 x0");
    //     }
    //     else
    //     {
    //         hKtLambdaArr[cent]->Draw("p e1 x0 same");
    //     }
    // }
    // canvKtLambda->BuildLegend();

    std::vector<TH1D*> hYRinvArr,hYLambdaArr;
    
    for (std::size_t cent = 0; cent < centralitySize; ++cent)
    {
        TH1D *hYRinv = new TH1D(TString::Format("hYRinv_%ld",cent),";y^{1,2}_{c.m.};R_{inv} [fm]",yArrSize,yVal);
        //TH1D *hYLambda = new TH1D(TString::Format("hYLambda_%ld",cent),";y_{c.m.} [a.u.];#lambda [a.u.]",yArrSize,yVal);
        for (std::size_t y = 1; y <= yArrSize; ++y)
        {
            hYRinv->SetBinContent(y,yRinvVal[cent][y-1]);
            hYRinv->SetBinError(y,yRinvErr[cent][y-1]);

            // hYLambda->SetBinContent(y,yLambdaVal[cent][y-1]);
            // hYLambda->SetBinError(y,yLambdaErr[cent][y-1]);
        }

        MakeNice(hYRinv);
        //MakeNice(hYLambda);

        hYRinvArr.push_back(hYRinv);
        //hYLambdaArr.push_back(hYLambda);
    }

    TCanvas *canvYRinv = new TCanvas("canvYRinv","",800,450);
    canvYRinv->SetMargin(0.15,0.05,0.15,0.05);
    for (std::size_t cent = 0; cent < centralitySize; ++cent)
    {
        hYRinvArr[cent]->SetMarkerColor(JJColor::fWutSecondaryColors[cent+1]);
        hYRinvArr[cent]->SetLineColor(JJColor::fWutSecondaryColors[cent+1]);
        hYRinvArr[cent]->SetTitle(centTitle[cent]);
        if(cent == 0)
        {
            hYRinvArr[cent]->GetYaxis()->SetRangeUser(3,3.8);
            hYRinvArr[cent]->Draw("p e");
        }
        else
        {
            hYRinvArr[cent]->Draw("p e same");
        }
    }
    //canvYRinv->BuildLegend();
    
    // TCanvas *canvYLambda = new TCanvas("canvYLambda","",800,450);
    // canvYLambda->SetMargin(0.2,0.02,0.15,0.02);
    // for (std::size_t cent = 0; cent < centralitySize; ++cent)
    // {
    //     hYLambdaArr[cent]->SetMarkerColor(JJColor::fWutSecondaryColors[cent+1]);
    //     hYLambdaArr[cent]->SetLineColor(JJColor::fWutSecondaryColors[cent+1]);
    //     hYLambdaArr[cent]->SetTitle(centTitle[cent]);
    //     if(cent == 0)
    //     {
    //         hYLambdaArr[cent]->GetYaxis()->SetRangeUser(0.,1.);
    //         hYLambdaArr[cent]->Draw("p e1 x0");
    //     }
    //     else
    //     {
    //         hYLambdaArr[cent]->Draw("p e1 x0 same");
    //     }
    // }
    // canvYLambda->BuildLegend();

    TFile *outputFile = TFile::Open("../output/1Dfit_0_10.root","RECREATE");
    canvKtRinv->Write();
    //canvKtLambda->Write();
    canvYRinv->Write();
    //canvYLambda->Write();

    canvKtRinv->SaveAs("canvKtRinv.svg");
    //canvKtLambda->SaveAs("canvKtLambda.svg");
    canvYRinv->SaveAs("canvYRinv.svg");
    //canvYLambda->SaveAs("canvYLambda.svg");

    for (auto hist : hKtRinvArr)
        hist->Write();
    // for (auto hist : hKtLambdaArr)
    //     hist->Write();
    for(auto hist : hYRinvArr)
        hist->Write();
    // for(auto hist : hYLambdaArr)
    //     hist->Write();
}