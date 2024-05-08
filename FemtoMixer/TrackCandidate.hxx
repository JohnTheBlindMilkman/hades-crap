#ifndef TrackCandidate_hxx
    #define TrackCandidate_hxx

#include "TLorentzVector.h"
#include "TCutG.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"

namespace Selection
{
    enum class Detector {RPC, ToF};

    class TrackCandidate
    {
        friend class PairCandidate;

        private:
            std::string TrackId;
            Detector System;
            bool IsAtMdcEdge;
            short int PID, Charge;
            float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass2, Beta, PolarAngle, AzimuthalAngle;
            std::array<std::vector<int>,24> FiredWireList; // this now has drastically increaded the memory usage of each track... - JJ

            /**
             * @brief Fill the list of fired cells in all MDC layers which create this track
             * 
             * @param innerSeg HMdcSeg object pointer to the inner MDC segments (before magnet)
             * @param outerSeg HMdcSeg object pointer to the outer MDC segments (after magnet)
             * @param wireList reference to the wire list
             */
            void SetWires(HMdcSeg* innerSeg, HMdcSeg* outerSeg, std::array<std::vector<int>,24> &wireList)
            {
                const int segLayers = 12;
                int ncells = 0;

                if (innerSeg == nullptr || outerSeg == nullptr)
                    return;

                for (int layer = 0; layer < segLayers; ++layer)
                {
                    ncells = innerSeg->getNCells(layer);
                    if (ncells < 1)
                        continue;

                    for (int cell = 0; cell < ncells; ++cell)
                        wireList.at(layer).push_back(innerSeg->getCell(layer,cell));
                }
                for (int layer = segLayers; layer < segLayers+segLayers; ++layer)
                {
                    ncells = outerSeg->getNCells(layer);
                    if (ncells < 1)
                        continue;

                    for (int cell = 0; cell < ncells; ++cell)
                        wireList.at(layer).push_back(outerSeg->getCell(layer,cell));
                }
            }

        public:
            TrackCandidate(){}
            /**
             * @brief Construct a new Track Candidate object
             * 
             * @param particleCand HParticleCand object pointer
             * @param innerSeg HMdcSeg object pointer to the inner MDC segments (before magnet), set to nullptr if MDC segment info is not needed
             * @param outerSeg HMdcSeg object pointer to the outer MDC segments (after magnet), set to nullptr if MDC segment info is not needed
             * @param evtId unique ID of the underlying event
             * @param trackId unique ID of this track within the event (e.g. for (int i = 0; i < nTracks; ++i) trackId = i;)
             * @param pid PID of the particle we want (when using DSTs put here whatever, just make sure the same PID is in the TrackCandidate::SelectTrack method)
             */
            TrackCandidate(HParticleCand* particleCand, HMdcSeg* innerSeg, HMdcSeg* outerSeg, const std::string &evtId, const std::size_t &trackId, const short int &pid)
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
                PID = pid; // only for simulations
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
                SetWires(innerSeg,outerSeg,FiredWireList);
            }
            /**
             * @brief Construct a new Track Candidate object
             * 
             * @param particleCand HParticleCandSim object pointer (contains PID)
             * @param innerSeg HMdcSeg object pointer to the inner MDC segments (before magnet)
             * @param outerSeg HMdcSeg object pointer to the outer MDC segments (after magnet)
             * @param evtId unique ID of the underlying event
             * @param trackId unique ID of this track within the event (e.g. for (int i = 0; i < nTracks; ++i) trackId = i;)
             */
            TrackCandidate(HParticleCandSim* particleCand, HMdcSeg* innerSeg, HMdcSeg* outerSeg, const std::string &evtId, const std::size_t &trackId)
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
                SetWires(innerSeg,outerSeg,FiredWireList);
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
             * @return true if track is selected
             * @return false otherwise
             */
            bool SelectTrack(TCutG *rpcCut, TCutG *tofCut)
            {
                if (PID != 14)
                    return false;
                if (IsAtMdcEdge)
                    return false;
                if (TotalMomentum > 2000.)
                    return false;
                if (TransverseMomentum > 1400.)
                    return false;
                if (Beta < 0.2)
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
             * @brief Prints on the standard output the indexes of all the wires that were hit
             * 
             */
            void PrintWires()
            {
                std::cout << "\n\n---=== Inner Segment ===---\n";
                for (int i = 0; i < 12; ++i)
                {
                    std::cout << "Layer " << i+1 << ":\t";
                    for (const auto &wire : FiredWireList.at(i))
                        std::cout << wire << "\t";

                    std::cout << "\n";
                }

                std::cout << "---=== Outer Segment ===---\n";
                for (int i = 12; i < 24; ++i)
                {
                    std::cout << "Layer " << i+1 << ":\t";
                    for (const auto &wire : FiredWireList.at(i))
                        std::cout << wire << "\t";

                    std::cout << "\n";
                }
                std::cout << "\n\n";
            }
            /**
             * @brief Returns combined hash for this track. Required when used inside an std::map 
             * 
             * @return std::size_t 
             */
            std::size_t GetHash() const
            {
                return std::hash<std::string>{}(TrackId);
            }
            std::string GetID() const
            {
                return TrackId;
            }
            /**
             * @brief Get the Wire List object
             * 
             * @return std::array<std::vector<int>,24> 
             */
            std::array<std::vector<int>,24> GetWireList() const
            {
                return FiredWireList;
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
             * @brief Check if two tracks are equal
             * 
             * @param other 
             * @return true 
             * @return false 
             */
            bool operator==(const TrackCandidate &other) const
            {
                return (TrackId == other.TrackId);
            }
    };
} // namespace Selection

    template<>
    struct std::hash<Selection::TrackCandidate>
    {
        std::size_t operator()(const Selection::TrackCandidate &trck) const
        {
            return trck.GetHash();
        }
    };

#endif
