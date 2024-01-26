#ifndef EventCandidate_hxx
    #define EventCandidate_hxx

#include "TrackCandidate.hxx"

#include <vector>

namespace Selection
{
    struct EventCandidate
    {
        int Centrality;
        short int TargetPlate;
        float ReactionPlaneAngle;
        float X, Y, Z;
        std::vector<TrackCandidate> trackList;

        bool operator==(const EventCandidate &other) const
        {
            return (Centrality == other.Centrality && 
            TargetPlate == other.TargetPlate && 
            ReactionPlaneAngle == other.ReactionPlaneAngle &&
            abs(X - other.X) < std::numeric_limits<float>::epsilon() && 
            abs(Y - other.Y) < std::numeric_limits<float>::epsilon() &&
            abs(Z - other.Z) < std::numeric_limits<float>::epsilon() &&
            trackList == other.trackList);
        }

        bool operator!=(const EventCandidate &other)
        {
            return !(*this == other);
        }

    };

    void CreateEvent(EventCandidate &eventCand, const float &vertx, const float &verty, const float &vertz, const int &cent, const float &EP)
    {
        eventCand.Centrality = cent;
        eventCand.trackList.clear();
        eventCand.trackList.resize(0);
        eventCand.X = vertx;
        eventCand.Y = verty;
        eventCand.Z = vertz;
        eventCand.ReactionPlaneAngle = EP;
    }

    bool SelectEvent(EventCandidate &eventCand, const int centIndex = 1, const int nSigma = 1)
    {
        // copied from zVertexPeakAndSigma.txt file
        // first: mean, second: stdev
        static const std::vector<std::pair<float,float>> plateVector = {{-54.7598,0.755846},{-51.6971,0.783591},{-47.7996,0.763387},{-44.5473,0.769386},{-40.569,0.781312},{-37.2151,0.762538},{-33.2948,0.76901},{-30.3726,0.742618},{-26.648,0.748409},{-22.5492,0.738462},{-18.9649,0.747727},{-15.5259,0.749724},{-11.8726,0.740386},{-8.45083,0.742672},{-4.58076,0.712394}};
        static const std::pair<float,float> xPosition = {0.1951,0.619};
        static const std::pair<float,float> yPosition = {0.7045,0.6187};

        if (eventCand.Centrality != centIndex)
            return false;

        if ((eventCand.X < (xPosition.first - nSigma*xPosition.second)) || (eventCand.X > (xPosition.first + nSigma*xPosition.second)))
            return false;
        if ((eventCand.Y < (yPosition.first - nSigma*yPosition.second)) || (eventCand.Y > (yPosition.first + nSigma*yPosition.second)))
            return false;

        short int plateNo = 0;
        for(const auto &plate : plateVector)
        {
            if ((eventCand.Z < (plate.first - nSigma*plate.second)) || (eventCand.Z > (plate.first + nSigma*plate.second)))
            {
                eventCand.TargetPlate = plateNo;
                return true;
            }
            plateNo++;
        }

        return false;
    }
} // namespace Selection

    template<>
    struct std::hash<Selection::EventCandidate>
    {
        std::size_t operator()(const Selection::EventCandidate &evt) const
        {
            std::size_t res = 0;
            boost::hash_combine(res,evt.Centrality);
            boost::hash_combine(res,evt.ReactionPlaneAngle);
            boost::hash_combine(res,evt.TargetPlate);
            boost::hash_combine(res,evt.X);
            boost::hash_combine(res,evt.Y);
            boost::hash_combine(res,evt.Z);
            boost::hash_combine(res,evt.trackList);

            return res;
        }
    };

    template<> 
    std::size_t boost::hash_value(std::vector<Selection::TrackCandidate> const &trckVec)
    {
        std::size_t res = 17;
        for (const auto &track : trckVec)
            res = res *31 + std::hash<Selection::TrackCandidate>()(track);

        return res;
    }

#endif