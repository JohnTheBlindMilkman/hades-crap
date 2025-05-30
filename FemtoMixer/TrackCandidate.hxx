#ifndef TrackCandidate_hxx
    #define TrackCandidate_hxx

#include "../../HFiredWires/HFiredWires.hxx"

#include "TLorentzVector.h"
#include "TCutG.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hparticlemetamatcher.h"

namespace Selection
{
    enum class Detector {RPC, ToF};

    class TrackCandidate
    {
        friend class PairCandidate; // this is here because I have a badly structured code

        private:
            std::shared_ptr<TrackCandidate> GeantKineTrack;
            std::string TrackId;
            Detector System;
            bool IsAtMdcEdge;
            short int PID, Charge, Sector;
            unsigned NBadLayers;
            float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass, Mass2, Beta, PolarAngle, AzimuthalAngle, AzimuthalAngleWrtEP, ReactionPlaneAngle;
            std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> firedWiresCollection;
            std::array<short unsigned,HADES::MDC::WireInfo::numberOfPlanes> goodLayers;
            std::vector<unsigned> metaHits;

            /**
             * @brief Remove and mark MDC layers which are considered "bad", i.e. the track has fired too many wires in it
             * 
             * @param list collection of wires
             * @param cutoff lower bound of acceptable number of fired wires in each MDC layer 
             */
            void RemoveBadLayers(std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> &list, std::size_t cutoff)
            {
                for (auto &layer : firedWiresCollection)
                    if (layer.size() > cutoff)
                    {
                        ++NBadLayers;
                        layer.clear();
                        layer.resize(0);
                    }
            }
            /**
             * @brief calculate in how many layers a hit has been registered for each plane
             * 
             * @param list 
             */
            void CalculateLayersPerPlane(const std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> &list)
            {
                for (const auto &plane : HADES::MDC::WireInfo::planeIndexing)
                {
                    short unsigned counter = 0;
                    for (const auto &layer : HADES::MDC::WireInfo::allLayerPerPlaneIndexing.at(plane))
                    {
                        if (list[layer].size() > 0)
                            ++counter;
                    }
                    goodLayers.at(plane) = counter;
                }
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part 
             * @return std::vector<unsigned> 
             */
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HParticleCand* part) const
            {
                std::vector<unsigned> output;
                int mod,cell;

                switch (System)
                {
                    case Detector::RPC:
                        //mod = part->getMetaModule(0);
                        cell = part->getMetaCell(0);
                        if (cell > -1)
                            output.push_back(static_cast<unsigned>(64 + (31 - cell)));

                        //mod = part->getMetaModule(1);
                        cell = part->getMetaCell(1);
                        if (cell > -1)
                            output.push_back(static_cast<unsigned>(64 + (31 - cell)));
                        break;

                    case Detector::ToF:
                        mod = part->getMetaModule(0);
                        cell = part->getMetaCell(0);
                        if (mod > -1 && cell > -1)
                            output.push_back(static_cast<unsigned>(mod * 8 + cell));

                        mod = part->getMetaModule(1);
                        cell = part->getMetaCell(1);
                        if (mod > -1 && cell > -1)
                            output.push_back(static_cast<unsigned>(mod * 8 + cell));
                        break;
                    
                    default:
                        break;
                }

                return output;
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part 
             * @return std::vector<unsigned> 
             */
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HParticleCandSim* part) const
            {
                std::vector<unsigned> output;
                int mod,cell;

                switch (System)
                {
                    case Detector::RPC:
                        //mod = part->getMetaModule(0);
                        cell = part->getMetaCell(0);
                        if (cell > -1)
                            output.push_back(static_cast<unsigned>(64 + (31 - cell)));

                        //mod = part->getMetaModule(1);
                        cell = part->getMetaCell(1);
                        if (cell > -1)
                            output.push_back(static_cast<unsigned>(64 + (31 - cell)));
                        break;

                    case Detector::ToF:
                        mod = part->getMetaModule(0);
                        cell = part->getMetaCell(0);
                        if (mod > -1 && cell > -1)
                            output.push_back(static_cast<unsigned>(mod * 8 + cell));

                        mod = part->getMetaModule(1);
                        cell = part->getMetaCell(1);
                        if (mod > -1 && cell > -1)
                            output.push_back(static_cast<unsigned>(mod * 8 + cell));
                        break;
                    
                    default:
                        break;
                }
            }
            [[nodiscard]] std::vector<unsigned> CalcMetaHits(HGeantKine* part) const
            {
                std::vector<unsigned> output;
                std::vector<HGeantRpc*> rpcCells;
                std::vector<HGeantTof*> tofCells;
                int mod,cell;

                switch (System)
                {
                    case Detector::RPC:
                        if (part->getRpcHits(rpcCells))
                        {
                            for (const auto elem : rpcCells)
                            {
                                cell = elem->getCell();
                                if (cell > -1)
                                    output.push_back(static_cast<unsigned>(64 + (31 - cell)));
                            }
                        }
                        break;

                    case Detector::ToF:
                        if (part->getTofHits(tofCells))
                        {
                            for (const auto elem : tofCells)
                            {
                                mod = elem->getModule();
                                cell = elem->getCell();
                                if (mod > -1 && cell > -1)
                                    output.push_back(static_cast<unsigned>(mod * 8 + cell));
                            }
                        }
                        break;
                    
                    default:
                        break;
                }

                return output;
            }
            /**
             * @brief Gets rid of the angle wrap when calculating pair azimuthal angle. Used phi range is (-202.5,157.5] (I need this for asHBT, to have a in-plane and out-of-plane bin)
             * 
             * @param angle 
             * @return float 
             */
            [[nodiscard]] float ConstrainAngle(const float &angle)
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
             * @param wires array of the fired wires, obtained from HFiredWires class
             * @param evtId unique ID of the underlying event
             * @param EP event plane angle of the underlying event
             * @param trackId unique ID of this track within the event (e.g. index of the track in the loop)
             * @param pid PID of the particle we want (when using DSTs put here whatever, just make sure the same PID is in the TrackCandidate::SelectTrack method)
             */
            TrackCandidate(HParticleCand* particleCand, const std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> &wires, const std::string &evtId, float EP, std::size_t trackId, short pid)
            : GeantKineTrack(nullptr), ReactionPlaneAngle(EP),firedWiresCollection(wires),NBadLayers(0),goodLayers({}),metaHits({})
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
                TLorentzVector vecTmp = *particleCand;
                CalculateLayersPerPlane(firedWiresCollection);
                RemoveBadLayers(firedWiresCollection,2);

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhi();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhi() - ReactionPlaneAngle);
                Beta = particleCand->getBeta();
                Charge = particleCand->getCharge();
                Sector = particleCand->getSector();
                Energy = vecTmp.E();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
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
                metaHits = CalcMetaHits(particleCand);
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
            TrackCandidate(HParticleCandSim* particleCand,HGeantKine* geantKine, const std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> &wires, const std::string &evtId, float EP, std::size_t trackId, short pid) 
            : GeantKineTrack((geantKine == nullptr) ? nullptr : new TrackCandidate(geantKine,evtId,EP,trackId,pid)), 
            ReactionPlaneAngle(EP),firedWiresCollection(wires),NBadLayers(0),goodLayers({}),metaHits({})
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(particleCand->getGeantPID()));
                TLorentzVector vecTmp = *particleCand;
                CalculateLayersPerPlane(firedWiresCollection);
                RemoveBadLayers(firedWiresCollection,2);

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhi();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhi() - ReactionPlaneAngle);
                Beta = particleCand->getBeta();
                Charge = particleCand->getCharge();
                Sector = particleCand->getSector();
                Energy = vecTmp.E();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
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
                metaHits = CalcMetaHits(particleCand);
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
            TrackCandidate(HGeantKine* particleCand, const std::string &evtId, float EP, std::size_t trackId, short pid) 
            : GeantKineTrack(nullptr), ReactionPlaneAngle(EP),firedWiresCollection({}),NBadLayers(0),goodLayers({}),metaHits({})
            {
                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhiDeg();
                AzimuthalAngleWrtEP = ConstrainAngle(particleCand->getPhiDeg() - ReactionPlaneAngle);
                Energy = particleCand->getE();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
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
                metaHits = CalcMetaHits(particleCand);
            }
            /**
             * @brief Destroy the Track Candidate object
             * 
             */
            ~TrackCandidate(){}
            /**
             * @brief Track selection method
             * 
             * @param rpcCut object pointer to the RPC bannana cut
             * @param tofCut object pointer to the ToF bannana cut
             * @param checkPID set a flag to additionally select tracks on PID
             * @return true if track is selected
             * @return false otherwise
             */
            bool SelectTrack(const TCutG *rpcCut, const TCutG *tofCut, bool checkPID = true)
            {
                if (PID != 14 && checkPID)
                    return false;
                if (IsAtMdcEdge)
                    return false;
                /* if (TotalMomentum < 500. || TotalMomentum > 1500.)
                    return false; */
                /*if (TransverseMomentum < 300. || TransverseMomentum > 1000.)
                    return false;
                if (Rapidity < 0.14 || Rapidity > 1.34) // (-0.6,0.6) in c.m.
                    return false;*/
                /* if (Beta < 0.2)
                    return false; */
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
             * @brief Create a modernised wire collection object
             * 
             * @param wi HParticleWireInfo object, obtained from HParticleMetaMatcher
             * @return std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> 
             */
            [[nodiscard]] static std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> CreateWireArray(const HParticleWireInfo &wi)
            {
                std::array<std::vector<unsigned>,HADES::MDC::WireInfo::numberOfAllLayers> array;

                for (std::size_t mod = 0; mod < HADES::MDC::WireInfo::numberOfPlanes; ++mod)
                    for (std::size_t lay = 0; lay < HADES::MDC::WireInfo::numberOfLayersInPlane; ++lay)
                    {
                        std::vector<unsigned> tmpVec;

                        if (wi.ar[mod][lay][0] > -1)
                            tmpVec.push_back(static_cast<unsigned>(wi.ar[mod][lay][0]));
                        if (wi.ar[mod][lay][1] > -1)
                            tmpVec.push_back(static_cast<unsigned>(wi.ar[mod][lay][1]));

                        array[mod * HADES::MDC::WireInfo::numberOfLayersInPlane + lay] = tmpVec;
                    }

                return array;
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
             * @brief Set the 4-momentum
             * 
             * @param px x component
             * @param py y component
             * @param pz z component
             * @param ene energy
             */
            void SetMomentum(float px, float py, float pz, float ene) noexcept
            {
                Px = px;
                Py = py;
                Pz = pz;
                Energy = ene;
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
    };
} // namespace Selection

#endif
