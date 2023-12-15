#include "TLorentzVector.h"
#include "TCutG.h"

namespace Selection
{
    enum class Detector {RPC, ToF};

    struct TrackCandidate
    {
        Detector System;
        bool IsAtMdcEdge;
        short int PID, Charge;
        float Rapidity, TotalMomentum, TransverseMomentum, Px, Py, Pz, Energy, Mass2, Beta, PolarAngle, AzimimuthalAngle;
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
