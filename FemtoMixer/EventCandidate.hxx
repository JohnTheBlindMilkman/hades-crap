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
    };
} // namespace Selection
