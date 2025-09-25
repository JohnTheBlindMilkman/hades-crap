#ifndef MacroUtils_hxx
    #define MacroUtils_hxx

    #include "TH1.h"
    #include "TH2.h"
    #include "TCanvas.h"
    #include "TString.h"
    #include "TPavesText.h"
    #include "TMath.h"

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

            std::vector<TH1D*> MergeHistograms(const std::vector<TH1D*> &inpHistos, const std::vector<std::vector<std::size_t> > &groups)
            {
                std::vector<TH1D*> otpHistos(groups.size(),nullptr);
                
                for (std::size_t i = 0; i < groups.size(); ++i)
                {
                    for (const auto &j : groups.at(i))
                    {
                        if (j == groups.at(i).front())
                        {
                            otpHistos.at(i) = inpHistos.at(j);
                        }
                        else
                        {
                            otpHistos.at(i)->Add(inpHistos.at(j));
                        }
                    }
                }

                return otpHistos;
            }

        } // namespace CF

        namespace Generic
        {
            /**
             * @brief Set the errors of hout to include full correlation between the numerator and denominator for any distribution such that hout = hNum / hDen. Yes, hNum and hDen both should be const but for some reason the method TH1::FindBin is not const...
             * 
             * @param hout 
             * @param hNum 
             * @param hDen 
             */
            void SetErrorsDivide(TH1 *hout, const TH1 *hNum, TH1 *hDen)
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
                    vErr = std::sqrt((eNum*eNum)/(vDen*vDen) + ((vNum*vNum)*(eDen*eDen))/(vDen*vDen*vDen*vDen) - (2*vNum*eNum*eDen)/(vDen*vDen*vDen));
                    
                    (TMath::IsNaN(vErr)) ? hout->SetBinError(i,0) : hout->SetBinError(i,vErr);
                }
            }

            /**
             * @brief Set the errors of hout to include full correlation between the multiplicaton for any distribution such that hout = hLhs * hRhs. Yes, hLhs and hRhs both should be const but for some reason the method TH1::FindBin is not const...
             * 
             * @param hout 
             * @param hLhs 
             * @param hRhs 
             */
            void SetErrorsMultiply(TH1 *hout, const TH1 *hLhs, TH1 *hRhs)
            {
                const int iterMax = hout->GetNbinsX();
                double vErr = 0, vLhs = 0, vRhs = 0, eLhs = 0, eRhs = 0;
                for (int i = 1; i <= iterMax; i++)
                {
                    vErr = 0;
                    vLhs = hLhs->GetBinContent(i);
                    eLhs = hLhs->GetBinError(i);
                    vRhs = hRhs->GetBinContent(hRhs->FindBin(hLhs->GetBinCenter(i)));
                    eRhs = hRhs->GetBinError(hRhs->FindBin(hLhs->GetBinCenter(i)));

                    // propagation of uncertainty for a function lhs*rhs with inclusion of the correlation between the constituents
                    vErr = std::sqrt(std::pow(vRhs * eLhs,2) + std::pow(vLhs * eRhs,2) + (2 * vLhs * vRhs * eLhs * eRhs));
                    hout->SetBinError(i,vErr);
                }
            }
        } // namespace Generic

        namespace Drawing
        {
            /**
             * @brief Struct representing the size and position of rectangle with a diagonal between (x1,y1) and (x2,y2)
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
                TPaveText *pt = new TPaveText(pos.x1,pos.y1,pos.x2,pos.y2,"NB");
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

            /**
             * @brief TCanvas partitioning function. Allows to divide a TCanvas for applications when we want to draw pads with the same X and Y axis, but only tle left most and bottom plots should have axis lablels ant ticks. Taken from https://root.cern/doc/v636/canvas2_8C.html
             * 
             * @param C TCanvas we want to divide
             * @param Nx Number of divisions in the X axis
             * @param Ny Number of divisions in the Y axis
             * @param lMargin Left margin of the canvas
             * @param rMargin Right margin of the canvas
             * @param bMargin Bottom margin of the canvas
             * @param tMargin Top margin of the canvas
             */
            void CanvasPartition(TCanvas *C, const Int_t Nx, const Int_t Ny, Float_t lMargin, Float_t rMargin, Float_t bMargin, Float_t tMargin)
            {
                if (!C)
                    return;
                
                // Setup Pad layout:
                Float_t vSpacing = 0.0;
                Float_t vStep = (1. - bMargin - tMargin - (Ny - 1) * vSpacing) / Ny;
                
                Float_t hSpacing = 0.0;
                Float_t hStep = (1. - lMargin - rMargin - (Nx - 1) * hSpacing) / Nx;
                
                Float_t vposd, vposu, vmard, vmaru, vfactor;
                Float_t hposl, hposr, hmarl, hmarr, hfactor;
                
                for (Int_t i = 0; i < Nx; i++) {
                
                    if (i == 0) {
                        hposl = 0.0;
                        hposr = lMargin + hStep;
                        hfactor = hposr - hposl;
                        hmarl = lMargin / hfactor;
                        hmarr = 0.0;
                    } else if (i == Nx - 1) {
                        hposl = hposr + hSpacing;
                        hposr = hposl + hStep + rMargin;
                        hfactor = hposr - hposl;
                        hmarl = 0.0;
                        hmarr = rMargin / (hposr - hposl);
                    } else {
                        hposl = hposr + hSpacing;
                        hposr = hposl + hStep;
                        hfactor = hposr - hposl;
                        hmarl = 0.0;
                        hmarr = 0.0;
                    }
                
                    for (Int_t j = 0; j < Ny; j++) {
                
                        if (j == 0) {
                            vposd = 0.0;
                            vposu = bMargin + vStep;
                            vfactor = vposu - vposd;
                            vmard = bMargin / vfactor;
                            vmaru = 0.0;
                        } else if (j == Ny - 1) {
                            vposd = vposu + vSpacing;
                            vposu = vposd + vStep + tMargin;
                            vfactor = vposu - vposd;
                            vmard = 0.0;
                            vmaru = tMargin / (vposu - vposd);
                        } else {
                            vposd = vposu + vSpacing;
                            vposu = vposd + vStep;
                            vfactor = vposu - vposd;
                            vmard = 0.0;
                            vmaru = 0.0;
                        }
                
                        C->cd(0);
                
                        auto name = TString::Format("pad_%d_%d", i, j);
                        auto pad = (TPad *)C->FindObject(name.Data());
                        if (pad)
                            delete pad;
                        pad = new TPad(name.Data(), "", hposl, vposd, hposr, vposu);
                        pad->SetLeftMargin(hmarl);
                        pad->SetRightMargin(hmarr);
                        pad->SetBottomMargin(vmard);
                        pad->SetTopMargin(vmaru);
                
                        pad->SetFrameBorderMode(0);
                        pad->SetBorderMode(0);
                        pad->SetBorderSize(0);
                
                        pad->Draw();
                    }
                }
            }
            
            /**
             * @brief Helper function for placing graphics objects in the same position on X axis, regardless of the pad width. Taken from https://root.cern/doc/v636/canvas2_8C.html
             * 
             * @param x X coordinate we want
             * @return X coordinate for currently active pad needed to match x
             */
            double XtoPad(double x)
            {
                double xl, yl, xu, yu;
                gPad->GetPadPar(xl, yl, xu, yu);
                double pw = xu - xl;
                double lm = gPad->GetLeftMargin();
                double rm = gPad->GetRightMargin();
                double fw = pw - pw * lm - pw * rm;
                return (x * fw + pw * lm) / pw;
            }
            
            /**
             * @brief Helper function for placing graphics objects in the same position on Y axis, regardless of the pad height. Taken from https://root.cern/doc/v636/canvas2_8C.html
             * 
             * @param y Y coordinate we want
             * @return Y coordinate for currently active pad needed to match y
             */
            double YtoPad(double y)
            {
                double xl, yl, xu, yu;
                gPad->GetPadPar(xl, yl, xu, yu);
                double ph = yu - yl;
                double tm = gPad->GetTopMargin();
                double bm = gPad->GetBottomMargin();
                double fh = ph - ph * bm - ph * tm;
                return (y * fh + bm * ph) / ph;
            }
        } // namespace Drawing
        
    } // namespace JJUtils

#endif