/**
 * @file EventCandidate.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Custom event-like object. Stores information about the event and the list of tracks assigned to it.
 * @version 0.1.0
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef EventCandidate_hxx
    #define EventCandidate_hxx

#include "TrackCandidate.hxx"
#include "Target.hxx"

#include <vector>

#include "heventheader.h"
#include "hparticleevtinfo.h"
#include "hparticleevtchara.h"

namespace Selection
{
    class EventCandidate
    {
        private:
            std::string EventId;
            short int Centrality, TargetPlate, ChargedTracks;
            float ReactionPlaneAngle;
            float X, Y, Z;
            std::vector<std::shared_ptr<TrackCandidate> > trackList;

        public:
            /**
             * @brief Construct a new Event Candidate object
             * 
             */
            EventCandidate(){}
            /**
             * @brief Construct a new Event Candidate object
             * 
             * @param evtId unique ID of the event (use HEventHeader::getEventSeqNumber and HEventHeader::getEventRunNumber)
             * @param vertx x position of the event vertex
             * @param verty y position of the event vertex
             * @param vertz z position of the event vertex
             * @param cent event centrality class (use HParticleEvtChara::getCentralityClass)
             * @param EP reaction plane angle (use HParticleEvtChara::getEventPlane)
             */
            EventCandidate(const std::string &evtId, float vertx, float verty, float vertz, int cent, float EP)
            : EventId(evtId),Centrality(cent),TargetPlate(-1),ChargedTracks(0),ReactionPlaneAngle((EP < 0) ? 0 : TMath::RadToDeg() * EP),X(vertx),Y(verty),Z(vertz) {}
            /**
             * @brief Construct a new Event Candidate object
             * 
             * @param evtHeader pointer to HEventHeader object
             * @param evtInfo pointer to HParticleEvtInfo object
             * @param cent event centrality class (use HParticleEvtChara::getCentralityClass)
             * @param EP reaction plane angle (use HParticleEvtChara::getEventPlane)
             */
            EventCandidate(HEventHeader* evtHeader, const HParticleEvtInfo* evtInfo,  int cent, float EP)
                : EventId(std::to_string(evtHeader->getEventRunNumber()) + std::to_string(evtHeader->getEventSeqNumber())),
                Centrality(cent), TargetPlate(-1),
                ChargedTracks(evtInfo->getSumRpcMultHitCut() + evtInfo->getSumTofMultCut()),
                ReactionPlaneAngle((EP < 0) ? 0 : TMath::RadToDeg() * EP),
                X(evtHeader->getVertexReco().getPos().X()),
                Y(evtHeader->getVertexReco().getPos().Y()),
                Z(evtHeader->getVertexReco().getPos().Z()) {}
            /**
             * @brief Select event for given centrality and vertex position
             * 
             * @tparam T HADES target setup
             * @param centIndex desired centrality classes (same layout as from HParticleEvtChara)
             * @param nSigmaX how many sigmas from the mean plate position should be accepted for given plate in X direction
             * @param nSigmaY how many sigmas from the mean plate position should be accepted for given plate in Y direction
             * @param nSigmaZ how many sigmas from the mean plate position should be accepted for given plate in Z direction
             * @return true if event is accepted or false otherwise
             * @throws std::runtime_error if specified nSigmaZ is > 2 
             */
            template <HADES::Target::Setup T>
            bool SelectEvent(const std::vector<int> centIndex = {1}, const float nSigmaX = 1, const float nSigmaY = 1, const float nSigmaZ = 1)
            {
                if (nSigmaZ > 2)
                    throw std::runtime_error("nSigmaZ >2 overlaps between neighbouring plates, please choose smaller value.\n If the specified value was intentional, please evaluate youe life choices...");

                // if this event's centrality is not in the centrality index list, reject
                if (std::find(centIndex.begin(),centIndex.end(),Centrality) == centIndex.end())
                    return false;

                if ((X < (HADES::Target::GetXTargetPosition<T>().first - nSigmaX*HADES::Target::GetXTargetPosition<T>().second)) || 
                    (X > (HADES::Target::GetXTargetPosition<T>().first + nSigmaX*HADES::Target::GetXTargetPosition<T>().second)))
                    return false;
                if ((Y < (HADES::Target::GetYTargetPosition<T>().first - nSigmaY*HADES::Target::GetYTargetPosition<T>().second)) || 
                    (Y > (HADES::Target::GetYTargetPosition<T>().first + nSigmaY*HADES::Target::GetYTargetPosition<T>().second)))
                    return false;

                std::size_t plateNo = 0;
                for(const auto &plate : HADES::Target::GetZPlatePositions<T>())
                {
                    if ((Z >= (plate.first - nSigmaZ*plate.second)) && (Z <= (plate.first + nSigmaZ*plate.second)))
                    {
                        TargetPlate = plateNo;
                        return true;
                    }
                    plateNo++;
                }

                return false;
            }
            /**
             * @brief Returns the unique event ID
             * 
             * @return std::string 
             */
            [[nodiscard]] std::string GetID() const noexcept
            {
                return EventId;
            }
            /**
             * @brief Check if two events are not equal
             * 
             * @param other 
             * @return true 
             * @return false 
             */
            [[nodiscard]] bool operator!=(const EventCandidate &other) const noexcept
            {
                return (EventId != other.EventId);
            }
            /**
             * @brief Get the Centrality class (same layout as from HParticleEvtChara)
             * 
             * @return short int 
             */
            [[nodiscard]] short int GetCentrality() const noexcept
            {
                return Centrality;
            }
            /**
             * @brief Get the Plate number to which the event vertex is assigned to
             * 
             * @return short int 
             */
            [[nodiscard]] short int GetPlate() const noexcept
            {
                return TargetPlate;
            }
            /**
             * @brief Get the total number of chaged tracks (nToF + nRPC)
             * 
             * @return int 
             */
            [[nodiscard]] int GetNCharged() const noexcept
            {
                return ChargedTracks;
            }
            /**
             * @brief Get the list of tracks assigned th this EventCandidate
             * 
             * @return std::vector<TrackCandidate> 
             */
            [[nodiscard]] std::vector<std::shared_ptr<TrackCandidate> > GetTrackList() const noexcept
            {
                return trackList;
            }
            /**
             * @brief Get the size of the list of tracks assigned th this EventCandidate
             * 
             * @return std::size_t 
             */
            [[nodiscard]] std::size_t GetTrackListSize() const noexcept
            {
                return trackList.size();
            }
            /**
             * @brief Get the Reaction Plane angle (in deg)
             * 
             * @return float 
             */
            [[nodiscard]] float GetReactionPlane() const noexcept
            {
                return ReactionPlaneAngle;
            }
            /**
             * @brief Get the X component of the event vertex
             * 
             * @return float 
             */
            [[nodiscard]] float GetX() const noexcept
            {
                return X;
            }
            /**
             * @brief Get the Y component of the event vertex
             * 
             * @return float 
             */
            [[nodiscard]] float GetY() const noexcept
            {
                return Y;
            }
            /**
             * @brief Get the Z component of the event vertex
             * 
             * @return float 
             */
            [[nodiscard]] float GetZ() const noexcept
            {
                return Z;
            }
            /**
             * @brief Calcluate average pT of all tracks currently stored in this event
             * 
             * @return average pT or 0 if event is empty
             */
            [[nodiscard]] double GetAveragePt() const
            {
                if (trackList.empty())
                    return 0.;
                
                double avgPt = 0.;
                std::for_each(trackList.begin(),trackList.end(),[&avgPt](const std::shared_ptr<TrackCandidate> &track){avgPt += track->GetPt();});
                return  avgPt / trackList.size();
            }
            /**
             * @brief Add new TrackCandidate to this EventCandidate
             * 
             * @param trck 
             */
            void AddTrack(const std::shared_ptr<TrackCandidate> &trck)
            {
                trackList.push_back(trck);
            }
    };
} // namespace Selection

#endif