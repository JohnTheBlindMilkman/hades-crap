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
            std::string TrackId;
            Detector System;
            bool IsAtMdcEdge;
            short int PID, Charge;
            unsigned int NBadLayers;
            float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass2, Beta, PolarAngle, AzimuthalAngle, ReactionPlaneAngle;
            std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> firedWiresCollection;
            std::vector<unsigned> goodLayers,metaCells;

            /**
             * @brief Remove and mark MDC layers which are considered "bad", i.e. the track has fired too many wires in it
             * 
             * @param list 
             * @param cutoff 
             */
            void RemoveBadLayers(std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> &list, short int cutoff)
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
             * @brief calculate in how many layers a hit has been registered for this track
             * 
             * @param list 
             */
            void CalculateLayersPerPlane(const std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> &list)
            {
                for (const auto &plane : HADES::MDC::WireInfo::allLayerPerPlaneIndexing)
                {
                    unsigned counter = 0;
                    for (const auto &layer : plane)
                    {
                        if (list[layer].size() > 0)
                            ++counter;
                    }
                    goodLayers.push_back(counter);
                }
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part 
             * @return std::vector<unsigned> 
             */
            std::vector<unsigned> GetMetaHits(HParticleCand* part)
            {
                std::vector<unsigned> output;
                unsigned hit;

                hit = part->getMetaModule(0);
                if (hit != -1)
                    output.push_back(hit);

                hit = part->getMetaModule(1);
                if (hit != -1)
                    output.push_back(hit);

                return output;
            }
            /**
             * @brief Get the Meta Hits information for this track
             * 
             * @param part 
             * @return std::vector<unsigned> 
             */
            std::vector<unsigned> GetMetaHits(HParticleCandSim* part)
            {
                std::vector<unsigned> output;
                unsigned hit;

                hit = part->getMetaModule(0);
                if (hit != -1)
                    output.push_back(hit);

                hit = part->getMetaModule(1);
                if (hit != -1)
                    output.push_back(hit);

                return output;
            }
            std::vector<unsigned> GetMetaHits(HGeantKine* part)
            {
                std::vector<unsigned> output;
                std::vector<HGeantTof*> tofCells;
                std::vector<HGeantRpc*> rpcCells;
                unsigned hit;

                if (part->getTofHits(tofCells))
                {
                    for (const auto elem : tofCells)
                    {
                        hit = static_cast<unsigned>(elem->getModule());
                        if (hit >= 0)
                            output.push_back(hit);
                    }
                }
                else if (part->getRpcHits(rpcCells))
                {
                    for (const auto elem: rpcCells)
                    {
                        hit = static_cast<unsigned>(elem->getCell());
                        if (hit >= 0)
                            output.push_back(hit);
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
            float ConstrainAngle(const float &angle)
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
            TrackCandidate(HParticleCand* particleCand, const std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> &wires, const std::string &evtId, const float &EP, const std::size_t &trackId, const short int &pid)
            : ReactionPlaneAngle(EP),firedWiresCollection(wires),NBadLayers(0),goodLayers({}),metaCells({})
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
                TLorentzVector vecTmp = *particleCand;
                CalculateLayersPerPlane(firedWiresCollection);
                RemoveBadLayers(firedWiresCollection,2);

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = ConstrainAngle(particleCand->getPhi() - TMath::RadToDeg()*ReactionPlaneAngle);
                Beta = particleCand->getBeta();
                Charge = particleCand->getCharge();
                Energy = vecTmp.E();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
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
                metaCells = GetMetaHits(particleCand);
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
            TrackCandidate(HParticleCandSim* particleCand, const std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> &wires, const std::string &evtId, const float &EP, const std::size_t &trackId, const short int &pid) : 
            ReactionPlaneAngle(EP),firedWiresCollection(wires),NBadLayers(0),goodLayers({}),metaCells({})
            {
                particleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
                TLorentzVector vecTmp = *particleCand;

                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhi();
                Beta = particleCand->getBeta();
                Charge = particleCand->getCharge();
                Energy = vecTmp.E();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
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
                metaCells = GetMetaHits(particleCand);
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
            TrackCandidate(HGeantKine* particleCand, const std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> &wires, const std::string &evtId, const float &EP, const std::size_t &trackId, const short int &pid) : 
            ReactionPlaneAngle(EP),firedWiresCollection(wires),NBadLayers(0),goodLayers({}),metaCells({})
            {
                TrackId = evtId + std::to_string(trackId);
                AzimuthalAngle = particleCand->getPhiDeg();
                Energy = particleCand->getE();
                IsAtMdcEdge =particleCand->isAtAnyMdcEdge();
                Mass2 = particleCand->getM()*particleCand->getM();
                PID = particleCand->getID(); // only for simulations
                Charge = HPhysicsConstants::charge(PID);
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
                metaCells = GetMetaHits(particleCand);
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
                /* if (TotalMomentum > 2000.)
                    return false; */
                /* if (TransverseMomentum > 1400.)
                    return false; */
                /* if (Beta < 0.2)
                    return false; */
                if (NBadLayers > 1)
                    return false;
                /* if (std::count_if(goodLayers.begin(),goodLayers.end(),[](unsigned i){return (i > 3);}) != 4)
                    return false; */

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
            static std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> CreateWireArray(const HParticleWireInfo &wi)
            {
                static constexpr std::size_t nLayers{6}, nModules{4};
                std::array<std::vector<int>,HADES::MDC::WireInfo::numberOfAllLayers> array;

                for (std::size_t mod = 0; mod < nModules; ++mod)
                    for (std::size_t lay = 0; lay < nLayers; ++lay)
                    {
                        std::vector<int> tmpVec;

                        if (wi.ar[mod][lay][0] != -1)
                            tmpVec.push_back(wi.ar[mod][lay][0]);
                        if (wi.ar[mod][lay][1] != -1)
                            tmpVec.push_back(wi.ar[mod][lay][1]);

                        array[mod*nLayers + lay] = tmpVec;
                    }

                return array;
            }
            /**
             * @brief Get the indexes of wires at given layer
             * 
             * @param layer 
             * @return std::vector<int> 
             * @throw std::out_of_range is thrown when layer value exceeds the HADES::MDC::WireInfo::numberOfAllLayers-1 value
             */
            std::vector<int> GetWires(std::size_t layer) const
            {
                return firedWiresCollection.at(layer);
            }
            std::string GetID() const
            {
                return TrackId;
            }
            /**
             * @brief Get the detector which registered the track
             * 
             * @return Detector 
             */
            Detector GetSystem() const
            {
                return System;
            }
            /**
             * @brief Get the Number Of Bad Layers object
             * 
             * @return unsigned int 
             */
            unsigned int GetNumberOfBadLayers() const
            {
                return NBadLayers;
            }
            /**
             * @brief Get the azimuthan angle of the track (in deg)
             * 
             * @return float 
             */
            float GetPhi() const
            {
                return AzimuthalAngle;
            }
            /**
             * @brief Get the polar angle of the track (in deg)
             * 
             * @return float 
             */
            float GetTheta() const
            {
                return PolarAngle;
            }
            /**
             * @brief Get the Pt of the track
             * 
             * @return float 
             */
            float GetPt() const
            {
                return TransverseMomentum;
            }
            /**
             * @brief Get the Rapidity of the track
             * 
             * @return float 
             */
            float GetRapidity() const 
            {
                return Rapidity;
            }
            /**
             * @brief Get the Charge of the track
             * 
             * @return short int 
             */
            short int GetCharge() const
            {
                return Charge;
            }
            /**
             * @brief Get momentum  of the track
             * 
             * @return float 
             */
            float GetP() const
            {
                return TotalMomentum;
            }
            /**
             * @brief Set the 4-momentum
             * 
             * @param px x component
             * @param py y component
             * @param pz z component
             * @param ene energy
             */
            void SetMomentum(float px, float py, float pz, float ene)
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
            float GetBeta() const
            {
                return Beta;
            }
            /**
             * @brief Get squared mass of the track
             * 
             * @return float 
             */
            float GetM2() const
            {
                return Mass2;
            }
            /**
             * @brief Get the Reaction Plane object
             * 
             * @return float 
             */
            float GetReactionPlane() const
            {
                return ReactionPlaneAngle;
            }
            /**
             * @brief Get the Meta Cells object
             * 
             * @return std::vector<unsigned> 
             */
            std::vector<unsigned> GetMetaCells() const
            {
                return metaCells;
            }
    };
} // namespace Selection

#endif
