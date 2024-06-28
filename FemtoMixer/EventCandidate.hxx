#ifndef EventCandidate_hxx
    #define EventCandidate_hxx

#include "TrackCandidate.hxx"

#include <vector>

namespace Selection
{
    class EventCandidate
    {
        private:
            std::string EventId;
            short int Centrality,TargetPlate;
            float ReactionPlaneAngle;
            float X, Y, Z;
            std::vector<TrackCandidate> trackList;

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
            EventCandidate(const std::string &evtId, const float &vertx, const float &verty, const float &vertz, const int &cent, const float &EP)
            : EventId(evtId),Centrality(cent),ReactionPlaneAngle(EP),X(vertx),Y(verty),Z(vertz){}
            /**
             * @brief Destroy the Event Candidate object
             * 
             */
            ~EventCandidate(){}
            /**
             * @brief Select event for given centrality and vertex position
             * 
             * @param centIndex desired centrality class (same layout as from HParticleEvtChara)
             * @param nSigma how many sigmas from the mean plate position should be accepted for given plate in X,Y, and Z directions
             * @return true if event is accepted
             * @return false otherwise
             */
            bool SelectEvent(const int centIndex = 1, const int nSigma = 1)
            {
                if (nSigma >= 2)
                    throw std::runtime_error("nSigma >=2 overlaps between neighbouring plates, please choose smaller value.\n If the specified value was intentional, please evaluate youe life choices...");
                
                // copied from zVertexPeakAndSigma.txt file
                // first: mean, second: stdev
                static const std::vector<std::pair<float,float>> plateVector = {{-54.7598,0.755846},{-51.6971,0.783591},{-47.7996,0.763387},{-44.5473,0.769386},{-40.569,0.781312},{-37.2151,0.762538},{-33.2948,0.76901},{-30.3726,0.742618},{-26.648,0.748409},{-22.5492,0.738462},{-18.9649,0.747727},{-15.5259,0.749724},{-11.8726,0.740386},{-8.45083,0.742672},{-4.58076,0.712394}};
                static const std::pair<float,float> xPosition = {0.1951,0.619};
                static const std::pair<float,float> yPosition = {0.7045,0.6187};

                if (Centrality != centIndex)
                    return false;

                if ((X < (xPosition.first - nSigma*xPosition.second)) || (X > (xPosition.first + nSigma*xPosition.second)))
                    return false;
                if ((Y < (yPosition.first - nSigma*yPosition.second)) || (Y > (yPosition.first + nSigma*yPosition.second)))
                    return false;

                short int plateNo = 0;
                for(const auto &plate : plateVector)
                {
                    if ((Z >= (plate.first - nSigma*plate.second)) && (Z <= (plate.first + nSigma*plate.second)))
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
            std::string GetID() const
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
            bool operator!=(const EventCandidate &other) const
            {
                return (EventId != other.EventId);
            }
            /**
             * @brief Get the Centrality class (same layout as from HParticleEvtChara)
             * 
             * @return short int 
             */
            short int GetCentrality() const
            {
                return Centrality;
            }
            /**
             * @brief Get the Plate number to which the event vertex is assigned to
             * 
             * @return short int 
             */
            short int GetPlate() const
            {
                return TargetPlate;
            }
            /**
             * @brief Get the list of tracks assigned th this EventCandidate
             * 
             * @return std::vector<TrackCandidate> 
             */
            std::vector<TrackCandidate> GetTrackList() const
            {
                return trackList;
            }
            /**
             * @brief Get the list of tracks assigned th this EventCandidate
             * 
             * @return std::vector<TrackCandidate> 
             */
            std::vector<TrackCandidate> GetTrackList()
            {
                return trackList;
            }
            /**
             * @brief Get the size of the list of tracks assigned th this EventCandidate
             * 
             * @return std::size_t 
             */
            std::size_t GetTrackListSize() const
            {
                return trackList.size();
            }
            /**
             * @brief Get the Reaction Plane angle (in deg)
             * 
             * @return float 
             */
            float GetReactionPlane() const
            {
                return ReactionPlaneAngle;
            }
            /**
             * @brief Add new TrackCandidate to this EventCandidate
             * 
             * @param trck 
             */
            void AddTrack(const TrackCandidate &trck)
            {
                trackList.push_back(trck);
            }
    };
} // namespace Selection

#endif