#ifndef TrackCandidate_hxx
    #define TrackCandidate_hxx

#include "MdcWires.hxx"

#include "TLorentzVector.h"
#include "TCutG.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlemetamatcher.h"
#include "hgeantkine.h"

namespace Selection
{
    enum class Detector {RPC, ToF};

    class TrackCandidate
    {
        friend class PairCandidate; // this is here because I have a poorly structured code

        private:
            std::shared_ptr<TrackCandidate> GeantKineTrack;
            std::string TrackId;
            Detector System;
            bool isAtMdcEdge, isGoodMetaCell;
            short int PID, Charge, Sector;
            unsigned NBadLayers;
            float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass, Mass2, Beta, 
                PolarAngle, AzimuthalAngle, AzimuthalAngleWrtEP, ReactionPlaneAngle, innerSegChi2, outerSegChi2, metaMatchQuality, chi2;
            HADES::MDC::LayersTrack firedWiresCollection;
            std::array<short unsigned,HADES::MDC::WireInfo::numberOfPlanes> goodLayers;
            std::vector<unsigned> metaHits;

            static constexpr short badIndex = -1; // if there is no module or cell, the value is set to -1
            static constexpr short maxTofIndex = 64; // highest index of unique ToF cells
            static constexpr short maxRpcIndex = 250; // highest index of unique RPC cells
            static constexpr short tofCells = 8; // number of cells in each ToF module
            static constexpr short rpcCells = 31; // number of cells in each RPC collumn

            /**
             * @brief Create a new index for a given RPC hit, based on the module and cell numbers (values from maxTofIndex to maxRpcIndex - 1)
             * 
             * @param module ID of the module which was hit
             * @param cell ID of the cell which was hit
             * @return new index of the hit
             */
            [[nodiscard]] constexpr int MakeUniqueRpcIndex(int module, int cell) const noexcept
            {
                return maxTofIndex + module * rpcCells + cell;
            }
            /**
             * @brief Create a new index for a given ToF hit, based on the module and cell numbers (values from 0 to maxTofIndex - 1)
             * 
             * @param module ID of the module which was hit
             * @param cell ID of the cell which was hit
             * @return new index of the hit
             */
            [[nodiscard]] constexpr int MakeUniqueTofIndex(int module, int cell) const noexcept
            {
                return module * tofCells + cell;
            }
            /**
             * @brief Remove and mark MDC layers which are considered "bad", i.e. the track has fired too many wires in it
             * 
             * @param list collection of wires
             * @param cutoff maximal acceptable number of fired wires in each MDC layer 
             * @return number of bad layers
             */
            unsigned RemoveAndCountBadLayers(HADES::MDC::LayersTrack &list, std::size_t cutoff)
            {
                unsigned badLayers = 0;
                for (auto &layer : list)
                    if (layer.size() > cutoff)
                    {
                        ++badLayers;
                        layer.clear();
                        layer.resize(0);
                    }

                return badLayers;
            }
            /**
             * @brief calculate how many layers have been hit in each plane
             * 
             * @param list collection of wires
             * @return number of fired layers in each plane
             */
            [[nodiscard]] std::array<unsigned short,HADES::MDC::WireInfo::numberOfPlanes> CalculateLayersPerPlane(const HADES::MDC::LayersTrack &list) const
            {
                std::array<unsigned short,HADES::MDC::WireInfo::numberOfPlanes> tmpArray;

                for (const auto &plane : HADES::MDC::WireInfo::planeIndexing)
                {
                    short unsigned counter = 0;
                    for (const auto &layer : HADES::MDC::WireInfo::allLayerPerPlaneIndexing.at(plane))
                    {
                        if (list[layer].size() > 0)
                            ++counter;
                    }
                    tmpArray.at(plane) = counter;
                }

                return tmpArray;
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part HParticleCand pointer
             * @return vector of registered hits (each stored value represents a unique META cell; the vector::size() defines the number of hits for this track)
             */
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HParticleCand* part) const
            {
                std::vector<unsigned> output;
                int mod,cell,uniqueCell;

                if (System == Detector::RPC)
                {
                    for (const auto& hit : {0,1})
                    {
                        mod = part->getMetaModule(hit);
                        cell = part->getMetaCell(hit);
                        uniqueCell = MakeUniqueRpcIndex(mod,cell);
                        if (mod > badIndex && cell > badIndex && uniqueCell < maxRpcIndex)
                            output.push_back(static_cast<unsigned>(uniqueCell)); 
                    }
                }
                else if (System == Detector::ToF)
                {
                    for (const auto& hit : {0,1})
                    {
                        mod = part->getMetaModule(hit);
                        cell = part->getMetaCell(hit);
                        uniqueCell = MakeUniqueTofIndex(mod,cell);
                        if (mod > badIndex && cell > badIndex && uniqueCell < maxTofIndex)
                            output.push_back(static_cast<unsigned>(uniqueCell));
                    }
                }

                return output;
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part HParticleCandSim pointer
             * @return vector of registered hits (each stored value represents a unique META cell; the vector::size() defines the number of hits for this track)
             */
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HParticleCandSim* part) const
            {
                std::vector<unsigned> output;
                int mod,cell,uniqueCell;

                if (System == Detector::RPC)
                {
                    for (const auto& hit : {0,1})
                    {
                        mod = part->getMetaModule(hit);
                        cell = part->getMetaCell(hit);
                        uniqueCell = MakeUniqueRpcIndex(mod,cell);
                        if (mod > badIndex && cell > badIndex && uniqueCell < maxRpcIndex)
                            output.push_back(static_cast<unsigned>(uniqueCell)); 
                    }
                }
                else if (System == Detector::ToF)
                {
                    for (const auto& hit : {0,1})
                    {
                        mod = part->getMetaModule(hit);
                        cell = part->getMetaCell(hit);
                        uniqueCell = MakeUniqueTofIndex(mod,cell);
                        if (mod > badIndex && cell > badIndex && uniqueCell < maxTofIndex)
                            output.push_back(static_cast<unsigned>(uniqueCell));
                    }
                }

                return output;
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part HGeantKine pointer
             * @return vector of registered hits (each stored value represents a unique META cell; the vector::size() defines the number of hits for this track)
             */
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HGeantKine* part) const
            {
                std::vector<unsigned> output;
                int mod,cell,uniqueCell;

                if (System == Detector::RPC)
                {
                    std::vector<HGeantRpc*> rpcCellList;
                    if (part->getRpcHits(rpcCellList))
                    {
                        for (const auto elem : rpcCellList)
                        {
                            mod = elem->getColumn();
                            cell = elem->getCell();
                            uniqueCell = MakeUniqueRpcIndex(mod,cell);
                            if (mod > badIndex && cell > badIndex && uniqueCell < maxRpcIndex)
                                output.push_back(static_cast<unsigned>(uniqueCell));
                        }
                    }
                }
                else if (System == Detector::ToF)
                {
                    std::vector<HGeantTof*> tofCellList;
                    if (part->getTofHits(tofCellList))
                    {
                        for (const auto elem : tofCellList)
                        {
                            mod = elem->getModule();
                            cell = elem->getCell();
                            uniqueCell = MakeUniqueTofIndex(mod,cell);
                            if (mod > badIndex && cell > badIndex && uniqueCell < maxTofIndex)
                                output.push_back(static_cast<unsigned>(uniqueCell));
                        }
                    }
                }

                return output;
            }
            /**
             * @brief Gets rid of the angle wrap when calculating pair azimuthal angle. Used phi range is (-202.5,157.5] (I need this for asHBT, to have a in-plane and out-of-plane bin)
             * 
             * @param angle 
             * @return float 
             */
            [[nodiscard]] float ConstrainAngle(float angle) const noexcept
            {
                if (angle > 157.5)
                    return angle -360.;
                else if (angle <= -202.5)
                    return angle + 360.;
                else
                    return angle;
            }

        public:
            TrackCandidate(){}
            /**
             * @brief Construct a new Track Candidate object
             * 
             * @param particleCand HParticleCand object pointer
             * @param wires array of the fired wires
             * @param evtId unique ID of the underlying event
             * @param EP event plane angle of the underlying event
             * @param trackId unique ID of this track within the event (e.g. index of the track in the loop)
             * @param pid PID of the particle we want (when using DSTs put here whatever, just make sure the same PID is in the TrackCandidate::SelectTrack method)
             */
            TrackCandidate(HParticleCand* particleCand, const HADES::MDC::LayersTrack &wires, const std::string &evtId, float EP, std::size_t trackId, short pid) : 
                GeantKineTrack(nullptr), ReactionPlaneAngle(EP),NBadLayers(0),firedWiresCollection(wires),
                goodLayers(CalculateLayersPerPlane(firedWiresCollection)),metaHits(CalcMetaHits(particleCand))
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
                TLorentzVector vecTmp = *particleCand;
                NBadLayers = RemoveAndCountBadLayers(firedWiresCollection,2);

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhi();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhi() - ReactionPlaneAngle);
                Beta = particleCand->getBeta(); // or should I use this one Beta = 1 - (1/(1+(TotalMomentum*TotalMomentum/Mass2))); ?
                Charge = particleCand->getCharge();
                Sector = particleCand->getSector();
                Energy = vecTmp.E();
                isAtMdcEdge =particleCand->isAtAnyMdcEdge();
                isGoodMetaCell = HParticleTool::isGoodMetaCell(particleCand,4,kTRUE);
                Mass = particleCand->getMass();
                Mass2 = particleCand->getMass2();
                PID = pid;
                PolarAngle = particleCand->getTheta();
                Px = vecTmp.Px();
                Py = vecTmp.Py();
                Pz = vecTmp.Pz();
                Rapidity = vecTmp.Rapidity();
                if (particleCand->getSystem() == 0)
                    System = Detector::RPC;
                else if (particleCand->getSystem() == 1)
                    System = Detector::ToF;
                TotalMomentum = vecTmp.P();
                TransverseMomentum = vecTmp.Pt();

                innerSegChi2 = particleCand->getInnerSegmentChi2();
                outerSegChi2 = particleCand->getOuterSegmentChi2();
                metaMatchQuality = particleCand->getMetaMatchQuality();
                chi2 = particleCand->getChi2();
            }
            /**
             * @brief Construct a new Track Candidate object
             * 
             * @param particleCand HParticleCandSim object pointer (contains PID)
             * @param wires 
             * @param evtId unique ID of the underlying event
             * @param EP event plane angle of the underlying event
             * @param trackId unique ID of this track within the event (e.g. index of the track in the loop)
             * @param pid PID of the particle we want (put here whatever, we use HParticleCandSim for PID later)
             */
            TrackCandidate(HParticleCandSim* particleCand,HGeantKine* geantKine, const HADES::MDC::LayersTrack &wires, const std::string &evtId, float EP, std::size_t trackId, short pid) : 
                GeantKineTrack((geantKine == nullptr) ? nullptr : new TrackCandidate(geantKine,evtId,EP,trackId,pid)), 
                ReactionPlaneAngle(EP),NBadLayers(0),firedWiresCollection(wires),goodLayers(CalculateLayersPerPlane(firedWiresCollection)),
                metaHits(CalcMetaHits(particleCand))
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(particleCand->getGeantPID()));
                TLorentzVector vecTmp = *particleCand;
                NBadLayers = RemoveAndCountBadLayers(firedWiresCollection,2);

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhi();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhi() - ReactionPlaneAngle);
                Beta = particleCand->getBeta(); // or should I use this one Beta = 1 - (1/(1+(TotalMomentum*TotalMomentum/Mass2))); ?
                Charge = particleCand->getCharge();
                Sector = particleCand->getSector();
                Energy = vecTmp.E();
                isAtMdcEdge =particleCand->isAtAnyMdcEdge();
                isGoodMetaCell = HParticleTool::isGoodMetaCell(particleCand,4,kTRUE);
                Mass = particleCand->getMass();
                Mass2 = particleCand->getMass2();
                PID = particleCand->getGeantPID(); // only for simulations
                PolarAngle = particleCand->getTheta();
                Px = vecTmp.Px();
                Py = vecTmp.Py();
                Pz = vecTmp.Pz();
                Rapidity = vecTmp.Rapidity();
                if (particleCand->getSystem() == 0)
                    System = Detector::RPC;
                else if (particleCand->getSystem() == 1)
                    System = Detector::ToF;
                TotalMomentum = vecTmp.P();
                TransverseMomentum = vecTmp.Pt();
                
                innerSegChi2 = particleCand->getInnerSegmentChi2();
                outerSegChi2 = particleCand->getOuterSegmentChi2();
                metaMatchQuality = particleCand->getMetaMatchQuality();
                chi2 = particleCand->getChi2();
            }
            /**
             * @brief Construct a new Track Candidate object
             * 
             * @param particleCand HGeantKine object pointer (contains PID)
             * @param wires 
             * @param evtId unique ID of the underlying event
             * @param EP event plane angle of the underlying event
             * @param trackId unique ID of this track within the event (e.g. index of the track in the loop)
             * @param pid PID of the particle we want (put here whatever, we use HParticleCandSim for PID later)
             */
            TrackCandidate(HGeantKine* particleCand, const std::string &evtId, float EP, std::size_t trackId, short pid) :
                GeantKineTrack(nullptr), ReactionPlaneAngle(EP),innerSegChi2(std::numeric_limits<float>::max()), 
                outerSegChi2(std::numeric_limits<float>::max()), metaMatchQuality(std::numeric_limits<float>::max()), 
                chi2(std::numeric_limits<float>::max()),NBadLayers(0),firedWiresCollection({}),goodLayers({}),
                metaHits(CalcMetaHits(particleCand))
            {
                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhiDeg();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhiDeg() - ReactionPlaneAngle);
                Energy = particleCand->getE();
                isAtMdcEdge =particleCand->isAtAnyMdcEdge();
                Mass = particleCand->getM();
                Mass2 = Mass * Mass;
                PID = particleCand->getID(); // only for simulations
                Charge = HPhysicsConstants::charge(PID);
                Sector = particleCand->getSector();
                PolarAngle = particleCand->getThetaDeg();
                particleCand->getMomentum(Px,Py,Pz);
                Rapidity = particleCand->getRapidity();
                if (particleCand->getSystem() == 0)
                    System = Detector::RPC;
                else if (particleCand->getSystem() == 1)
                    System = Detector::ToF;
                TotalMomentum = particleCand->getTotalMomentum();
                TransverseMomentum = particleCand->getTransverseMomentum();
                Beta = 1 - (1/(1+(TotalMomentum*TotalMomentum/Mass2))); // beta = 1 - 1/(1+p^2/m_0^2) if my calculations are correct
            }
            /**
             * @brief Track selection method
             * 
             * @param rpcCut object pointer to the RPC bannana cut
             * @param tofCut object pointer to the ToF bannana cut
             * @param checkPID set a flag to additionally select tracks on PID
             * @return true if track is selected
             * @return false otherwise
             */
            bool SelectTrack(const TCutG *rpcCut, const TCutG *tofCut, bool checkPID = true) const
            {
                if (PID != 14 && checkPID)
                    return false;
                if (isAtMdcEdge)
                    return false;
                // if (!isGoodMetaCell)
                //     return false;
                // if (innerSegChi2 <= 0 || outerSegChi2 <= 0)
                //     return false;
                // if (metaMatchQuality >= 4)
                //     return false;
                // if (chi2 >= 400)
                //     return false;
                if (NBadLayers > 1)
                    return false;
                if (std::count_if(goodLayers.begin(),goodLayers.end(),[](unsigned i){return (i > 3);}) != 4)
                    return false;

                switch (System)
                {
                    case Detector::RPC:
                        if (rpcCut->IsInside(TotalMomentum*Charge,Beta))
                            return true;
                        break;

                    case Detector::ToF:
                        if (tofCut->IsInside(TotalMomentum*Charge,Beta))
                            return true;
                        break;
                }
                
                return false;
            }
            /**
             * @brief Get the indexes of wires at given layer
             * 
             * @param layer 
             * @return std::vector<unsigned> 
             * @throw std::out_of_range is thrown when layer value exceeds the HADES::MDC::WireInfo::numberOfAllLayers-1 value
             */
            [[nodiscard]] std::vector<unsigned> GetWires(std::size_t layer) const
            {
                return firedWiresCollection.at(layer);
            }
            /**
             * @brief Get the All Wires object
             * 
             * @return HADES::MDC::LayersTrack 
             */
            [[nodiscard]] HADES::MDC::LayersTrack GetAllWires() const noexcept
            {
                return firedWiresCollection;
            }
            /**
             * @brief Get the unique ID of this track
             * 
             * @return std::string 
             */
            [[nodiscard]] std::string GetID() const noexcept
            {
                return TrackId;
            }
            /**
             * @brief Get the detector which registered the track
             * 
             * @return Detector 
             */
            [[nodiscard]] Detector GetSystem() const noexcept
            {
                return System;
            }
            /**
             * @brief Get the Number Of Bad Layers object
             * 
             * @return unsigned
             */
            [[nodiscard]] unsigned GetNumberOfBadLayers() const noexcept
            {
                return NBadLayers;
            }
            /**
             * @brief Get the azimuthan angle of the track (in deg).
             * 
             * @return float 
             */
            [[nodiscard]] float GetPhi() const noexcept
            {
                return AzimuthalAngle;
            }
            /**
             * @brief Get the azimuthal angle of the track with respect to the Event Plane (in deg). Ranging from -202.5 to 157.5
             * 
             * @return float 
             */
            [[nodiscard]] float GetPhiWrtEP() const noexcept
            {
                return AzimuthalAngleWrtEP;
            }
            /**
             * @brief Get the polar angle of the track (in deg)
             * 
             * @return float 
             */
            [[nodiscard]] float GetTheta() const noexcept
            {
                return PolarAngle;
            }
            /**
             * @brief Get the Pt of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetPt() const noexcept
            {
                return TransverseMomentum;
            }
            /**
             * @brief Get the Rapidity of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetRapidity() const noexcept
            {
                return Rapidity;
            }
            /**
             * @brief Get the Charge of the track
             * 
             * @return short int 
             */
            [[nodiscard]] short int GetCharge() const noexcept
            {
                return Charge;
            }
            /**
             * @brief Get the assigned sector for this track 
             * 
             * @return short int 
             */
            [[nodiscard]] short int GetSector() const noexcept
            {
                return Sector;
            }
            /**
             * @brief Get momentum  of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetP() const noexcept
            {
                return TotalMomentum;
            }
            /**
             * @brief Get momentum  of the track along X-axis
             * 
             * @return float 
             */
            [[nodiscard]] float GetPx() const noexcept
            {
                return Px;
            }
            /**
             * @brief Get momentum  of the track along Y-axis
             * 
             * @return float 
             */
            [[nodiscard]] float GetPy() const noexcept
            {
                return Py;
            }
            /**
             * @brief Get momentum  of the track along Z-axis
             * 
             * @return float 
             */
            [[nodiscard]] float GetPz() const noexcept
            {
                return Pz;
            }
            /**
             * @brief Get the Beta of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetBeta() const noexcept
            {
                return Beta;
            }
            /**
             * @brief Get mass of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetM() const noexcept
            {
                return Mass;
            }
            /**
             * @brief Get squared mass of the track
             * 
             * @return float 
             */
            [[nodiscard]] float GetM2() const noexcept
            {
                return Mass2;
            }
            /**
             * @brief Get the Reaction Plane (in deg)
             * 
             * @return float 
             */
            [[nodiscard]] float GetReactionPlane() const noexcept
            {
                return ReactionPlaneAngle;
            }
            /**
             * @brief Get all META hits
             * 
             * @return std::vector<unsigned> 
             */
            [[nodiscard]] std::vector<unsigned> GetMetaHits() const noexcept
            {
                return metaHits;
            }
            /**
             * @brief Get the stored HGeantKine TrackCandidate
             * 
             * @return std::shared_ptr<TrackCandidate> 
             */
            [[nodiscard]] const std::shared_ptr<TrackCandidate>& GetGeantKine() const noexcept
            {
                return GeantKineTrack;
            }
            /**
             * @brief Get the Inner Seg Chi 2 object
             * 
             * @return float 
             */
            [[nodiscard]] float GetInnerSegChi2() const noexcept
            {
                return innerSegChi2;
            }
            /**
             * @brief Get the Outer Seg Chi 2 object
             * 
             * @return float 
             */
            [[nodiscard]] float GetOuterSegChi2() const noexcept
            {
                return outerSegChi2;
            }
            /**
             * @brief Get the Chi 2 object
             * 
             * @return float 
             */
            [[nodiscard]] float GetChi2() const noexcept
            {
                return chi2;
            }
            /**
             * @brief Get the Meta Match Quality object
             * 
             * @return float 
             */
            [[nodiscard]] float GetMetaMatchQuality() const noexcept
            {
                return metaMatchQuality;
            }
            /**
             * @brief Check if the track is at the edge of MDC
             * 
             * @return true if it is and false otherwise
             */
            [[nodiscard]] bool IsAtMdcEdge() const noexcept
            {
                return isAtMdcEdge;
            }
            /**
             * @brief Check if the track has a valid META cell
             * 
             * @return true if it has and false otherwise
             */
            [[nodiscard]] bool IsGoodMetaCell() const noexcept
            {
                return isGoodMetaCell;
            }
    };
} // namespace Selection

#endif
