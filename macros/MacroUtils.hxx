#ifndef MacroUtils_hxx
    #define MacroUtils_hxx

    #include "TH1.h"
    #include "TH2.h"
    #include "TCanvas.h"
    #include "TString.h"
    #include "TPavesText.h"

    namespace JJUtils
    {
        namespace CF
        {
            /**
             * @brief Returns a value which can be used to normalise the correlation function at given range to unity.
             * 
             * @param hist correlation function which you want to normalise
             * @param xMin lower edge of the normalisation range
             * @param xMax upper edge of the normalisation range
             * @return norm for the histogram or 0 if no bins were found
             */
            double GetNormByRange(const TH1 *hist, double xMin, double xMax)
            {
                int nBins = 0;
                double val = 0., xVal;
                for (int i = 1; i < hist->GetNbinsX(); i++)
                {
                    xVal = hist->GetBinCenter(i);
                    if (xVal >= xMin && xVal <= xMax)
                    {
                        val += hist->GetBinContent(i);
                        nBins++;
                    } 
                }

                return (nBins > 0) ? val / nBins : 0.;
            }

            /**
             * @brief Returns a value which can be used to normalise the correlation function to unity.
             * 
             * @param hNum numerator of the correlation function
             * @param hDen denominator of the correlation function
             * @return norm for the histogram or 0 if no entries were found
             */
            double GetNormByEntires(const TH1 *hNum, const TH1 *hDen)
            {
                double entriesNum = hNum->GetEntries();
                return (entriesNum > 0.) ? hDen->GetEntries() / entriesNum : 0.;
            }

        } // namespace CF

        namespace Generic
        {
            namespace Detail
            {
                /**
                 * @brief Returns Pearson's Correlation Coefficient. Implemented from here:  https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
                 * 
                 * @param hX Histogram of data from one sample
                 * @param hY Histogram of data from another sample
                 * @return Value of the PCC
                 */
                double GetPCC(const TH1 *hX, const TH1 *hY)
                {
                    int nBinsX = hX->GetNbinsX();
                    int nBinsY = hY->GetNbinsX();

                    TH2D hNew("hNew","",nBinsX,hX->GetXaxis()->GetXmin(),hX->GetXaxis()->GetXmax(),nBinsY,hY->GetXaxis()->GetXmin(),hY->GetXaxis()->GetXmax());

                    for (int i = 1; i <= nBinsX; ++i)
                        for (int j = 1; j <= nBinsY; ++j)
                            hNew.SetBinContent(i,j,hX->GetBinContent(i) + hY->GetBinContent(j));

                    return hNew.GetCovariance() / (hX->GetStdDev()*hY->GetStdDev());
                }
            } // namespace Detail

            /**
             * @brief Set the errors of hout to include correlation between the numerator and denominator for any distribution such that hout = hNum / hDen. Yes, hNum and hDen should be const but for some reason the method TH1::FindBin is not const...
             * 
             * @param hout 
             * @param hNum 
             * @param hDen 
             */
            void SetErrors(TH1 *hout, TH1 *hNum, TH1 *hDen)
            {
                const int iterMax = hout->GetNbinsX();
                double vErr = 0, vNum = 0, vDen = 0, eNum = 0, eDen = 0;
                for (int i = 1; i <= iterMax; i++)
                {
                    vErr = 0;
                    vNum = hNum->GetBinContent(i);
                    eNum = hNum->GetBinError(i);
                    vDen = hDen->GetBinContent(hDen->FindBin(hNum->GetBinCenter(i)));
                    eDen = hDen->GetBinError(hDen->FindBin(hNum->GetBinCenter(i)));

                    // propagation of uncertainty for a function num/den with inclusion of the correlation between the constituents
                    if (fabs(vDen) > 0.)
                        vErr = std::sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - Detail::GetPCC(hNum,hDen)*(2*vNum*eNum*eDen)/(vDen*vDen*vDen));
                    hout->SetBinError(i,vErr);
                }
            }
            
        } // namespace Generic

        namespace Drawing
        {
            /**
             * @brief Struct representing the size and position of rectangle with diagonal between (x1,y1) and (x2,y2)
             * 
             */
            struct Position
            {
                double x1,y1,x2,y2;
            };

            /**
             * @brief Add TPaveText with given text at the defined position. No border, no fill, aligned to center.
             * 
             * @param text text you want to add
             * @param pos position of the TPaveText
             * @return new TPaveText with the newly assigned text
             */
            TPaveText* AddTextBox(const TString &text, Position &&pos = {0.2,0.8,0.5,0.9})
            {
                TPaveText *pt = new TPaveText(0.2,0.8,0.5,0.9,"NB");
                pt->AddText(text);
                pt->SetTextAlign(kHAlignCenter+kVAlignCenter);
                pt->SetFillStyle(1);
                pt->SetBorderSize(0);
                pt->SetTextSize(0.05);

                return pt;
            }

            /**
             * @brief Enclose TCanvas within another canvas and add the text specified in the title
             * 
             * @param c canvas which you want to draw
             * @param title text that will be added at the top of the canvas
             * @return new canvas, containig the old one with added text on top
             */
            TCanvas* AddFrame(TCanvas *c, const TString &title)
            {
                const float topMargin = 0.075;
                const float whModifier = 1 + topMargin;
                TCanvas *newCanv = new TCanvas(c->GetName(),"",c->GetWw(),c->GetWh()*whModifier);
                newCanv->SetMargin(0.,0.,0.,0.);

                TPad *top = new TPad("top","",0.,1.-topMargin,1.,1.);
                TPad *bottom = new TPad("bottom","",0.,0.,1.,1.-topMargin);

                TPaveText *pt = new TPaveText(c->GetLeftMargin(),0.01,(1. - c->GetRightMargin()),0.9,"NB");
                pt->AddText(title);
                pt->SetTextFont(62);
                pt->SetTextAlign(kHAlignRight+kVAlignCenter);
                pt->SetFillStyle(1);
                pt->SetBorderSize(0);

                top->cd();
                pt->Draw();

                bottom->cd();
                c->DrawClonePad();

                newCanv->cd();
                top->Draw();
                bottom->Draw();

                return newCanv;
            }
        } // namespace Drawing
        
    } // namespace JJUtils

#endif