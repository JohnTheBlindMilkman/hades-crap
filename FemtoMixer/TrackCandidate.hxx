#ifndef TrackCandidate_hxx
    #define TrackCandidate_hxx

#include "TLorentzVector.h"
#include "TCutG.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"

#include "boost/functional/hash.hpp"

namespace Selection
{
    enum class Detector {RPC, ToF};

    struct TrackCandidate
    {
        Detector System;
        bool IsAtMdcEdge;
        short int PID, Charge;
        float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass2, Beta, PolarAngle, AzimimuthalAngle;

        bool operator==(const TrackCandidate &other) const
        {
            return (System == other.System &&
            IsAtMdcEdge == other.IsAtMdcEdge &&
            PID == other.PID &&
            Charge == other.Charge &&
            abs(Rapidity - other.Rapidity) < std::numeric_limits<float>::epsilon() &&
            abs(TotalMomentum - other.TotalMomentum) < std::numeric_limits<float>::epsilon() &&
            abs(TransverseMomentum - other.TransverseMomentum) < std::numeric_limits<float>::epsilon() &&
            abs(Px - other.Px) < std::numeric_limits<float>::epsilon() &&
            abs(Py - other.Py) < std::numeric_limits<float>::epsilon() &&
            abs(Pz - other.Pz) < std::numeric_limits<float>::epsilon() &&
            abs(Energy - other.Energy) < std::numeric_limits<float>::epsilon() &&
            abs(Mass2 - other.Mass2) < std::numeric_limits<float>::epsilon() &&
            abs(Beta - other.Beta) < std::numeric_limits<float>::epsilon() &&
            abs(PolarAngle - other.PolarAngle) < std::numeric_limits<float>::epsilon() &&
            abs(AzimimuthalAngle - other.AzimimuthalAngle) < std::numeric_limits<float>::epsilon());
        }
    };

    void CreateTrack(TrackCandidate &trackCand, HParticleCand* particle_cand)
    {
        particle_cand->calc4vectorProperties(HPhysicsConstants::mass(14));
        TLorentzVector vecTmp = *particle_cand;

        trackCand.AzimimuthalAngle = particle_cand->getPhi();
        trackCand.Beta = particle_cand->getBeta();
        trackCand.Charge = particle_cand->getCharge();
        trackCand.Energy = vecTmp.E();
        trackCand.IsAtMdcEdge =particle_cand->isAtAnyMdcEdge();
        trackCand.Mass2 = particle_cand->getMass2();
        //trackCand.PID = particle_cand->getGeantPID(); // only for simulations
        trackCand.PolarAngle = particle_cand->getTheta();
        trackCand.Px = vecTmp.Px();
        trackCand.Py = vecTmp.Py();
        trackCand.Pz = vecTmp.Pz();
        trackCand.Rapidity = vecTmp.Rapidity();
        if (particle_cand->getSystem() == 0)
            trackCand.System = Detector::RPC;
        else if (particle_cand->getSystem() == 1)
            trackCand.System = Detector::ToF;
        trackCand.TotalMomentum = vecTmp.P();
        trackCand.TransverseMomentum = vecTmp.Pt();
    }

    bool SelectTrack(TrackCandidate &trackCand, TCutG *rpcCut, TCutG *tofCut)
    {
        //if (trackCand.PID != 14)
            //return false;
        if (trackCand.IsAtMdcEdge)
            return false;
        if (trackCand.TotalMomentum > 2000.)
            return false;
        if (trackCand.TransverseMomentum > 1400.)
            return false;
        if (trackCand.Beta < 0.2)
            return false;

        switch (trackCand.System)
        {
            case Detector::RPC:
                if (rpcCut->IsInside(trackCand.TotalMomentum*trackCand.Charge,trackCand.Beta))
                    return true;
                break;

            case Detector::ToF:
                if (tofCut->IsInside(trackCand.TotalMomentum*trackCand.Charge,trackCand.Beta))
                    return true;
                break;
        }
        
        return false;
    }

} // namespace Selection

    template<>
    struct std::hash<Selection::TrackCandidate>
    {
        std::size_t operator()(const Selection::TrackCandidate &trck) const
        {
            std::size_t res = 0;
            boost::hash_combine(res,trck.System);
            boost::hash_combine(res,trck.IsAtMdcEdge);
            boost::hash_combine(res,trck.PID);
            boost::hash_combine(res,trck.Charge);
            boost::hash_combine(res,trck.Rapidity);
            boost::hash_combine(res,trck.TotalMomentum);
            boost::hash_combine(res,trck.TransverseMomentum);
            boost::hash_combine(res,trck.Px);
            boost::hash_combine(res,trck.Py);
            boost::hash_combine(res,trck.Pz);
            boost::hash_combine(res,trck.Energy);
            boost::hash_combine(res,trck.Mass2);
            boost::hash_combine(res,trck.Beta);
            boost::hash_combine(res,trck.PolarAngle);
            boost::hash_combine(res,trck.AzimimuthalAngle);

            return res;
        }
    };

#endif