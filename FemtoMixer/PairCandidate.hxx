#ifndef PairCandidate_hxx
    #define PairCandidate_hxx

    #include "EventCandidate.hxx"

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
                PairCandidate(const std::shared_ptr<TrackCandidate> &trck1, const std::shared_ptr<TrackCandidate> &trck2) : 
                    pairId(trck1->GetID() + trck2->GetID()),
                    Particle1(trck1), 
                    Particle2(trck2), 
                    GeantKinePair((trck1->GeantKineTrack == nullptr || trck1->GeantKineTrack == nullptr) ? nullptr : new PairCandidate(trck1->GeantKineTrack,trck2->GeantKineTrack)), 
                    pairLayers(HADES::MDC::CreatePairLayers(trck1->GetAllWires(),trck2->GetAllWires())), 
                    wireDistances(HADES::MDC::CalculateWireDistances(pairLayers)),
                    SharedWires(HADES::MDC::CalculateSharedWires(pairLayers)), 
                    BothLayers(HADES::MDC::CalculateBothLayers(pairLayers)), 
                    SharedMetaCells(CalcSharedMetaCells(trck1,trck2)), 
                    MinWireDistance(*std::min_element(wireDistances.begin(),wireDistances.end())), 
                    QInv(0.), QOut(0.), QSide(0.), QLong(0.), Kt(0.), 
                    Rapidity((trck1->GetRapidity() + trck2->GetRapidity()) / 2.), 
                    AzimuthalAngle(ConstrainAngle(trck1->GetPhi() + trck2->GetPhi()) / 2.), 
                    OpeningAngle(CalcOpeningAngle(trck1,trck2)), 
                    DeltaPhi(trck1->GetPhi() - trck2->GetPhi()), 
                    DeltaTheta(trck1->GetTheta() - trck2->GetTheta()), 
                    SplittingLevel(HADES::MDC::CalcluateSplittingLevel(pairLayers)), 
                    areSameSector(trck1->GetSector() == trck2->GetSector())
                {
                    CFKinematics(trck1,trck2);
                }
                /**
                 * @brief Perform pair selection based on the fraction of neighbouring wires with certain distance from each other. Allows for a modifiable behaviour of selection
                 * 
                 * @tparam T Behaviour type
                 * @param fraction in what fraction of the all hits the merging can occur (values between 0 and 1 are allowed)
                 * @param cutoff how close the wires are allowed to be
                 * @return true 
                 * @return false 
                 */
                template<Behaviour T> 
                bool RejectPairByCloseHits(const float &fraction, const unsigned &cutoff) const
                {
                    if (fraction < 0. || fraction > 1.)
                    {
                        std::cerr << "error<PairCandidate::RejectPairByCloseHits>: provided fraction is outside the bounds [0,1]" << std::endl;
                        std::exit(1);
                    }

                    return Reject(type<T>{},fraction,cutoff);
                }
                /**
                 * @brief Returns unique ID of the pair
                 * 
                 * @return std::string 
                 */
                std::string GetID() const noexcept
                {
                    return pairId;
                }
                /**
                 * @brief Get the transverse component of the average pair momentum
                 * 
                 * @return float 
                 */
                float GetKt() const noexcept
                {
                    return Kt;
                }
                /**
                 * @brief Get the pair average rapidity
                 * 
                 * @return float 
                 */
                float GetRapidity() const noexcept
                {
                    return Rapidity;
                }
                /**
                 * @brief Get the pair average azimuthal angle
                 * 
                 * @return float 
                 */
                float GetPhi() const noexcept
                {
                    return AzimuthalAngle;
                }
                /**
                 * @brief Get difference of the azimuthal angles between the two pair components
                 * 
                 * @return float 
                 */
                float GetDPhi() const noexcept
                {
                    return DeltaPhi;
                }
                /**
                 * @brief Get difference of the polar angles between the two pair components
                 * 
                 * @return float 
                 */
                float GetDTheta() const noexcept
                {
                    return DeltaTheta;
                }
                /**
                 * @brief Get the Splitting Level object
                 * 
                 * @return float 
                 */
                float GetSplittingLevel() const noexcept
                {
                    return SplittingLevel;
                }
                /**
                 * @brief Get the Shared Wires object
                 * 
                 * @return unsigned int 
                 */
                unsigned GetSharedWires() const noexcept
                {
                    return SharedWires;
                }
                /**
                 * @brief Get the Both Layers object
                 * 
                 * @return unsigned 
                 */
                unsigned GetBothLayers() const noexcept
                {
                    return BothLayers;
                }
                /**
                 * @brief Get the All Layer Distances object
                 * 
                 * @return collection of distances between the wires in each MDC layer
                 */
                HADES::MDC::WireDistances GetAllLayerDistances() const noexcept
                {
                    return wireDistances;
                }
                /**
                 * @brief Get the Min Wire Distance object
                 * 
                 * @return unsigned 
                 */
                HADES::MDC::OptionalDistance<unsigned> GetMinWireDistance() const noexcept
                {
                    return MinWireDistance;
                }
                unsigned GetSharedMetaCells() const noexcept
                {
                    return SharedMetaCells;
                }
                /**
                 * @brief Get the Qinv
                 * 
                 * @return float 
                 */
                float GetQinv() const noexcept
                {
                    return QInv;
                }
                /**
                 * @brief Get the Qout, Qside, and Qlong components
                 * 
                 * @return std::tuple<float,float,float> 
                 */
                std::tuple<float,float,float> GetOSL() const noexcept
                {
                    return std::make_tuple(QOut,QSide,QLong);
                }
                /**
                 * @brief Get the Geant Kine Pair object
                 * 
                 * @return std::shared_ptr<PairCandidate> 
                 */
                [[nodiscard]] const std::shared_ptr<PairCandidate>& GetGeantKinePair() const noexcept
                {
                    return GeantKinePair;
                }
                /**
                 * @brief Check if the tracks are in thew same sector
                 * 
                 * @return true if they are and false otherwise
                 */
                [[nodiscard]] bool AreTracksFromTheSameSector() const noexcept
                {
                    return areSameSector;
                }
            private:
                std::string pairId;
                std::shared_ptr<TrackCandidate> Particle1,Particle2;
                std::shared_ptr<PairCandidate> GeantKinePair;
                HADES::MDC::LayersPair pairLayers;
                HADES::MDC::WireDistances wireDistances;
                unsigned SharedWires, BothLayers, SharedMetaCells;
                HADES::MDC::OptionalDistance<unsigned> MinWireDistance;
                float QInv, QOut, QSide, QLong, Kt, Rapidity, AzimuthalAngle, OpeningAngle, DeltaPhi, DeltaTheta, SplittingLevel;
                bool areSameSector;
                template <Behaviour T> struct type {}; // helper struct

                /**
                 * @brief Calculates the pair variables in their centre of mass system (here: LCMS)
                 * 
                 * @param part1 
                 * @param part2 
                 */
                void CFKinematics(const std::shared_ptr<TrackCandidate> &part1, const std::shared_ptr<TrackCandidate> &part2)
                {
                    // adapted from https://github.com/DanielWielanek/HAL/blob/main/analysis/femto/base/FemtoPairKinematics.cxx
                    double tPx = part1->Px + part2->Px;
                    double tPy = part1->Py + part2->Py;
                    double tPz = part1->Pz + part2->Pz;
                    double tE = part1->Energy + part2->Energy;
                    double tPt = tPx * tPx + tPy * tPy;
                    double tMt = tE * tE - tPz * tPz;  // mCVK;
                    tMt = std::sqrt(tMt);
                    Kt = std::sqrt(tPt);
                    double tBeta  = tPz / tE;
                    double tGamma = tE / tMt;

                    // Transform to LCMS

                    double particle1lcms_pz = tGamma * (part1->Pz - tBeta * part1->Energy);
                    double particle1lcms_e  = tGamma * (part1->Energy - tBeta * part1->Pz);
                    double particle2lcms_pz = tGamma * (part2->Pz - tBeta * part2->Energy);
                    double particle2lcms_e  = tGamma * (part2->Energy - tBeta * part2->Pz);

                    // Rotate in transverse plane

                    double particle1lcms_px = (part1->Px * tPx + part1->Py * tPy) / Kt;
                    double particle1lcms_py = (-part1->Px * tPy + part1->Py * tPx) / Kt;

                    double particle2lcms_px = (part2->Px * tPx + part2->Py * tPy) / Kt;
                    double particle2lcms_py = (-part2->Px * tPy + part2->Py * tPx) / Kt;

                    QOut = std::abs(particle1lcms_px - particle2lcms_px);
                    QSide = std::abs(particle1lcms_py - particle2lcms_py);
                    QLong = std::abs(particle1lcms_pz - particle2lcms_pz);
                    double mDE = particle1lcms_e - particle2lcms_e;
                    QInv = std::sqrt(std::abs(QOut * QOut + QSide * QSide + QLong * QLong - mDE * mDE));
                }
                /**
                 * @brief Calculates the opening angle in deg between the two tracks (reimplemented from TLorenzVector)
                 * 
                 * @param part1 
                 * @param part2 
                 * @return float 
                 */
                float CalcOpeningAngle(const std::shared_ptr<TrackCandidate> &part1, const std::shared_ptr<TrackCandidate> &part2) const
                {
                    float ptot = std::sqrt((part1->Px*part1->Px + part1->Py*part1->Py + part1->Pz*part1->Pz) * (part2->Px*part2->Px + part2->Py*part2->Py + part2->Pz*part2->Pz));
                    if (ptot <= 0.)
                    {	
                        return 0.;
                    }
                    else
                    {
                        float arg = (part1->Px*part2->Px + part1->Py*part2->Py + part1->Pz*part2->Pz) / ptot;
                        if (arg > 1.) return 1.;
                        if (arg < -1.) return -1.;
                        return std::acos(arg);
                    }
                }
                /**
                 * @brief Gets rid of the angle wrap when calculating pair azimuthal angle. Used phi range is (-202.5,157.5] (I need this for asHBT, to have a in-plane and out-of-plane bin)
                 * 
                 * @param angle 
                 * @return float 
                 */
                float ConstrainAngle(float angle) const noexcept
                {
                    if (angle > 157.5)
                    {
                        return angle -360.;
                    }
                    else if (angle <= -202.5)
                    {
                        return angle + 360.;
                    }
                    else
                    {
                        return angle;
                    }
                }
                /**
                 * @brief Calculates the number of shared META cells between the two tracks
                 * 
                 * @param part1 Pointer to the first track
                 * @param part2 Pointer to the second tracks
                 * @return number of shared META cells
                 */
                unsigned CalcSharedMetaCells(const std::shared_ptr<TrackCandidate> &part1, const std::shared_ptr<TrackCandidate> &part2) const noexcept
                {
                    unsigned SMC = 0;
                    auto meta1 = part1->GetMetaHits();
                    auto meta2 = part2->GetMetaHits();
                    
                    for (const auto &hit1 : meta1)
                        for (const auto &hit2 : meta2)
                        {
                            if (hit1 == hit2) ++SMC;
                        }

                    return SMC;
                }
                /**
                 * @brief Reject pair if the fraction of distances between fired wires is grater than the set limit
                 * 
                 * @param fractionLimit in what fraction of the all hits the merging can occur
                 * @param cutoff how close the wires are allowed to be
                 * @return true if the calculated fraction exceeds the limit
                 */
                bool Reject(type<Behaviour::Uniform>,const float &fractionLimit, const unsigned &cutoff) const noexcept
                {
                    float closeLayers = 0;
                    float totalLayers = 0;
                    for (const auto &dist : wireDistances)
                    {
                        if (dist.has_value)
                        {
                            if (dist.value < cutoff)
                            {
                                closeLayers += 1.;
                            }
                            totalLayers += 1.;
                        }
                    }
                    
                    if (totalLayers > 0)
                        return ( (closeLayers/totalLayers) > fractionLimit) ? true : false;
                    else
                        return true;
                }
                /**
                 * @brief Reject pair if the fraction of distances between fired wires is grater than the set limit or if at any layer the distance beteen the fired wires is smaller than cutoff - 1
                 * 
                 * @param fractionLimit in what fraction of the all hits the merging can occur
                 * @param cutoff how close the wires are allowed to be
                 * @return true if the calculated fraction exceeds the limit
                 */
                bool Reject(type<Behaviour::OneUnder>,const float &fractionLimit, const unsigned &cutoff) const noexcept
                {
                    float closeLayers = 0;
                    float totalLayers = 0;
                    for (const auto &dist : wireDistances)
                    {
                        if (dist.has_value)
                        {
                            if (dist.value < cutoff)
                            {
                                if (std::abs(static_cast<long>(dist.value) - cutoff) > 1)
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
                        return ( (closeLayers/totalLayers) > fractionLimit) ? true : false;
                    else
                        return true;
                }
                /**
                 * @brief Reject pair if the fraction of distances between fired wires is grater than the set limit. The calculated fraction is weighted based on how much closer the fired wires are. Currently this fraction can be > 1, IDK what to do about it, but I don't use it so...
                 * 
                 * @param fractionLimit in what fraction of the all hits the merging can occur
                 * @param cutoff how close the wires are allowed to be
                 * @return true if the calculated fraction exceeds the limit
                 */
                bool Reject(type<Behaviour::Weighted>,const float &fractionLimit, const unsigned &cutoff) const noexcept
                {
                    float closeLayers = 0.;
                    float totalLayers = 0.;
                    for (const auto &dist : wireDistances)
                    {
                        if (dist.has_value)
                        {
                            if (dist.value < cutoff)
                            {
                                closeLayers += std::abs(static_cast<long>(dist.value) - cutoff);
                            }
                            totalLayers += 1.;
                        }
                    }

                    if (totalLayers > 0)
                        return ( (closeLayers/totalLayers) > fractionLimit) ? true : false; // this fraction could be above 1, not sure how to interpret this
                    else
                        return true;
                }
        };
    }

#endif
