#ifndef PairCandidate_hxx
    #define PairCandidate_hxx

#include "EventCandidate.hxx"

#include <deque>
#include <random>
#include <iostream>

namespace Selection
{
    class PairCandidate
    {
        public:
            enum class Behaviour{OneUnder,Uniform,Weighted};
            /**
             * @brief Construct a new Pair Candidate object
             * 
             * @param trck1 
             * @param trck2 
             */
            PairCandidate(const TrackCandidate &trck1, const TrackCandidate &trck2) : Particle1(trck1), Particle2(trck2),SplittingLevel(0.),SharedWires(0)
            {
                CFKinematics(trck1,trck2);
                OpeningAngle = CalcOpeningAngle(trck1,trck2);
                Rapidity = (trck1.Rapidity + trck2.Rapidity) / 2.;
                AzimuthalAngle = (trck1.AzimuthalAngle + trck2.AzimuthalAngle) / 2.;
                DeltaPhi = trck1.AzimuthalAngle - trck2.AzimuthalAngle;
                DeltaTheta = trck1.PolarAngle - trck2.PolarAngle;
                pairId = trck1.GetID() + trck2.GetID();
                std::tie(SplittingLevel,SharedWires,BothLayers,MinWireDistance,wireDistances) = CalcSplittingLevelAndSharedWires(trck1,trck2);
                SharedMetaCells = CalcSharedMetaCells(trck1,trck2);
            }
            /**
             * @brief Destroy the Pair Candidate object
             * 
             */
            ~PairCandidate(){}
            /**
             * @brief Perform pair selection based on the fraction of neighbouring wires with certain distance from each other. Allows for a modifiable behaviour of selection
             * 
             * @tparam T Behaviour type: OneUnder, Uniform, and Weighted
             * @param fraction in what fraction of the all hits the merging can occur (values between 0 and 1 are allowed)
             * @param cutoff how close the wires are allowed to be
             * @return true 
             * @return false 
             */
            template<Behaviour T> bool RejectPairByCloseHits(const float &fraction,const unsigned &cutoff) const
            {
                if (fraction < 0. || fraction > 1.)
                {
                    std::cerr << "error<PairCandidate::RejectPairByCloseHits>: provided fraction is outside the bounds (0,1)" << std::endl;
                    std::exit(1);
                }

                return Reject(type<T>{},fraction,cutoff);
            }
            /**
             * @brief Returns unique ID of the pair
             * 
             * @return std::string 
             */
            std::string GetID() const
            {
                return pairId;
            }
            /**
             * @brief Get the transverse component of the average pair momentum
             * 
             * @return float 
             */
            float GetKt() const
            {
                return Kt;
            }
            /**
             * @brief Get the pair average rapidity
             * 
             * @return float 
             */
            float GetRapidity() const
            {
                return Rapidity;
            }
            /**
             * @brief Get the pair average azimuthal angle w.r.t. the EP (if specified in the constructor)
             * 
             * @return float 
             */
            float GetPhi() const
            {
                return AzimuthalAngle;
            }
            /**
             * @brief Get difference of the azimuthal angles between the two pair components
             * 
             * @return float 
             */
            float GetDPhi() const
            {
                return DeltaPhi;
            }
            /**
             * @brief Get difference of the polar angles between the two pair components
             * 
             * @return float 
             */
            float GetDTheta() const
            {
                return DeltaTheta;
            }
            /**
             * @brief Get the Splitting Level object
             * 
             * @return float 
             */
            float GetSplittingLevel() const
            {
                return SplittingLevel;
            }
            /**
             * @brief Get the Shared Wires object
             * 
             * @return unsigned int 
             */
            unsigned int GetSharedWires() const
            {
                return SharedWires;
            }
            /**
             * @brief Get the Both Layers object
             * 
             * @return unsigned 
             */
            unsigned GetBothLayers() const
            {
                return BothLayers;
            }
            /**
             * @brief Get the All Layer Distances object
             * 
             * @return std::array<int,HADES::MDC::WireInfo::numberOfAllLayers> 
             */
            std::array<int,HADES::MDC::WireInfo::numberOfAllLayers> GetAllLayerDistances() const
            {
                return wireDistances;
            }
            /**
             * @brief Get the Min Wire Distance object
             * 
             * @return unsigned 
             */
            unsigned GetMinWireDistance() const
            {
                return MinWireDistance;
            }
            unsigned GetSharedMetaCells() const
            {
                return SharedMetaCells;
            }
            /**
             * @brief Get the Qinv
             * 
             * @return float 
             */
            float GetQinv() const
            {
                return QInv;
            }
            /**
             * @brief Get the Qout, Qside, and Qlong components
             * 
             * @return std::tuple<float,float,float> 
             */
            std::tuple<float,float,float> GetOSL() const
            {
                return std::make_tuple(QOut,QSide,QLong);
            }
        private:
            std::string pairId;
            TrackCandidate Particle1,Particle2;
            unsigned SharedWires, BothLayers, MinWireDistance, SharedMetaCells;
            std::array<int,HADES::MDC::WireInfo::numberOfAllLayers> wireDistances;
            float QInv, QOut, QSide, QLong, Kt, Rapidity, AzimuthalAngle, OpeningAngle, DeltaPhi, DeltaTheta, SplittingLevel;
            template <Behaviour T> struct type {}; // helper struct

            /**
             * @brief Calculates the pair variables in their centre of mass system (here: LCMS)
             * 
             * @param part1 
             * @param part2 
             */
            void CFKinematics(const TrackCandidate &part1, const TrackCandidate &part2)
            {
                // adapted from https://github.com/DanielWielanek/HAL/blob/main/analysis/femto/base/FemtoPairKinematics.cxx
                Double_t tPx = part1.Px + part2.Px;
                Double_t tPy = part1.Py + part2.Py;
                Double_t tPz = part1.Pz + part2.Pz;
                Double_t tE = part1.Energy + part2.Energy;
                Double_t tPt = tPx * tPx + tPy * tPy;
                Double_t tMt = tE * tE - tPz * tPz;  // mCVK;
                tMt = TMath::Sqrt(tMt);
                Kt = TMath::Sqrt(tPt);
                Double_t tBeta  = tPz / tE;
                Double_t tGamma = tE / tMt;

                // Transform to LCMS

                Double_t particle1lcms_pz = tGamma * (part1.Pz - tBeta * part1.Energy);
                Double_t particle1lcms_e  = tGamma * (part1.Energy - tBeta * part1.Pz);
                Double_t particle2lcms_pz = tGamma * (part2.Pz - tBeta * part2.Energy);
                Double_t particle2lcms_e  = tGamma * (part2.Energy - tBeta * part2.Pz);

                // Rotate in transverse plane

                Double_t particle1lcms_px = (part1.Px * tPx + part1.Py * tPy) / Kt;
                Double_t particle1lcms_py = (-part1.Px * tPy + part1.Py * tPx) / Kt;

                Double_t particle2lcms_px = (part2.Px * tPx + part2.Py * tPy) / Kt;
                Double_t particle2lcms_py = (-part2.Px * tPy + part2.Py * tPx) / Kt;

                QOut = abs(particle1lcms_px - particle2lcms_px);
                QSide = abs(particle1lcms_py - particle2lcms_py);
                QLong = abs(particle1lcms_pz - particle2lcms_pz);
                Double_t mDE = particle1lcms_e - particle2lcms_e;
                QInv = TMath::Sqrt(TMath::Abs(QOut * QOut + QSide * QSide + QLong * QLong - mDE * mDE));
            }
            /**
             * @brief Calculates the opening angle in deg between the two tracks (reimplemented from TLorenzVector)
             * 
             * @param part1 
             * @param part2 
             * @return float 
             */
            float CalcOpeningAngle(const TrackCandidate &part1, const TrackCandidate &part2) const
            {
                float ptot = sqrt((part1.Px*part1.Px + part1.Py*part1.Py + part1.Pz*part1.Pz) * (part2.Px*part2.Px + part2.Py*part2.Py + part2.Pz*part2.Pz));
                if (ptot <= 0.)
                {	
                    return 0.;
                }
                else
                {
                    float arg = (part1.Px*part2.Px + part1.Py*part2.Py + part1.Pz*part2.Pz) / ptot;
                    if (arg > 1.) return 1.;
                    if (arg < -1.) return -1.;
                    return acos(arg);
                }
            }
            /**
             * @brief Gets rid of the angle wrap when calculating pair azimuthal angle. Used phi range is (-202.5,157.5] (I need this for asHBT, to have a in-plane and out-of-plane bin)
             * 
             * @param angle 
             * @return float 
             */
            float ConstrainAngle(const float &angle)
            {
                if (angle > 157.5)
                    return angle -360.;
                else if (angle <= -202.5)
                    return angle + 360.;
                else
                    return angle;
            }
            /**
             * @brief Calculates the splitting level (SL) and number of shared wires (SW) between two tracks. Those values are a measure of track merging and splitting
             * 
             * @param part1 
             * @param part2 
             * @return std::pair<float,unsigned int> first = SL, second = SW
             */
            std::tuple<float,unsigned,unsigned,unsigned,std::array<int,HADES::MDC::WireInfo::numberOfAllLayers> > CalcSplittingLevelAndSharedWires(const TrackCandidate &part1, const TrackCandidate &part2) const
            {
                float SL = 0;
                unsigned int SW = 0, MWD = 1000;
                std::array<int,HADES::MDC::WireInfo::numberOfAllLayers> allWD;
                std::fill(allWD.begin(),allWD.end(),-999); // we don't want to initialise with zeros, because it would mean that the unpopulated fields are sharing a wire
                int n0 = 0, n1 = 0, n2 = 0;
                for (const int &layer : HADES::MDC::WireInfo::allLayerIndexing)
                {
                    auto wires1 = part1.GetWires(layer);
                    auto wires2 = part2.GetWires(layer);
                    const std::size_t wires1Size = wires1.size();
                    const std::size_t wires2Size = wires2.size();

                    if (wires1Size > 0 && wires2Size > 0)
                    {
                        // iterate over all elements in both of the containers
                        for (const auto &wire1 : wires1)
                            for (const auto &wire2 : wires2)
                            {
                                if (wire1 == wire2)
                                {
                                    SL += 1.;
                                    ++SW;
                                }
                                else
                                {
                                    SL -= 1.;
                                }

                                if (abs(wire1 - wire2) < MWD)
                                    MWD = abs(wire1 - wire2);
                            }

                        allWD.at(layer) = MWD;
                        ++n0;
                        ++n1;
                        ++n2;
                    }
                    else if (wires1Size > 0 && wires2Size == 0)
                    {
                        SL += 1.;
                        ++n1;
                    }
                    else if (wires1Size == 0 && wires2Size > 0)
                    {
                        SL +=1.;
                        ++n2;
                    }
                }

                SL /= (n1 + n2);
                return std::make_tuple(SL,SW,n0,MWD,allWD); // this should be a struct at this point...
            }
            unsigned CalcSharedMetaCells(const TrackCandidate &part1, const TrackCandidate &part2) const
            {
                unsigned SMC = 0;
                std::size_t i = 0, j = 0;
                auto meta1 = part1.GetMetaCells();
                auto meta2 = part2.GetMetaCells();
                std::size_t meta1size = meta1.size();
                std::size_t meta2size = meta2.size();

                // iterate over all elements in both of the containers
                while (i < meta1size && j < meta2size) // this method gives me O(NlogN) complexity, but is probably overengineered (I have max 2 entries in vector)
                {
                    if (meta1[i] == meta2[j])
                    {
                        ++i;
                        ++j;
                        ++SMC;
                    }
                    else if (meta1[i] < meta2[j])
                    {
                        ++i;
                    }
                    else
                    {
                        ++j;
                    }
                }

                return SMC;
            }
            bool Reject(type<Behaviour::Uniform>,const float &fraction, const unsigned &cutoff) const
            {
                float closeLayers = 0;
                float totalLayers = 0;
                for (const int &dist : wireDistances)
                {
                    if (dist >= 0)
                    {
                        if (dist < cutoff)
                        {
                            closeLayers += 1.;
                        }
                        totalLayers += 1.;
                    }
                }
                
                if (totalLayers > 0)
                    return ( (closeLayers/totalLayers) > fraction) ? true : false;
                else
                    return true;
            }
            bool Reject(type<Behaviour::OneUnder>,const float &fraction, const unsigned &cutoff) const
            {
                float closeLayers = 0;
                float totalLayers = 0;
                for (const int &dist : wireDistances)
                {
                    if (dist >= 0)
                    {
                        if (dist < cutoff)
                        {
                            if (abs(dist - static_cast<int>(cutoff)) > 1) // cast to avoid mismatch type warning, I know it may round up the number, but if you put number big enough for this to happen you already did someting wrong
                            {    
                                return true;
                            }
                            else
                            {
                                closeLayers += 1.;
                            }
                        }
                        totalLayers += 1.;
                    }
                }

                if (totalLayers > 0)
                    return ( (closeLayers/totalLayers) > fraction) ? true : false;
                else
                    return true;
            }
            bool Reject(type<Behaviour::Weighted>,const float &fraction, const unsigned &cutoff) const
            {
                float closeLayers = 0.;
                float totalLayers = 0.;
                for (const int &dist : wireDistances)
                {
                    if (dist >= 0)
                    {
                        if (dist < cutoff)
                        {
                            closeLayers += abs(dist - static_cast<int>(cutoff)); // ibid.
                        }
                        totalLayers += 1.;
                    }
                }

                if (totalLayers > 0)
                    return ( (closeLayers/totalLayers) > fraction) ? true : false; // this fraction could be above 1, not sure how to interpret this
                else
                    return true;
            }
    };
}

#endif
