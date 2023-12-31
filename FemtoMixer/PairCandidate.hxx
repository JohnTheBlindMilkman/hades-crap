#ifndef PairCandidate_hxx
    #define PairCandidate_hxx

#include "EventCandidate.hxx"

#include <deque>
#include <random>

namespace Selection
{
    struct PairCandidate
    {
        TrackCandidate Particle1, Particle2;
        float QInv, QOut, QSide, QLong, Kt, Rapidity, AzimuthalAngle, OpeningAngle, DeltaPhi, DeltaTheta;
    };
        
    PairCandidate CFKinematics(const TrackCandidate &part1, const TrackCandidate &part2)
    {
        PairCandidate pair;
        pair.Particle1 = part1;
        pair.Particle2 = part2;

        // adapted from https://github.com/DanielWielanek/HAL/blob/main/analysis/femto/base/FemtoPairKinematics.cxx
        Double_t tPx = part1.Px + part2.Px;
        Double_t tPy = part1.Py + part2.Py;
        Double_t tPz = part1.Pz + part2.Pz;
        Double_t tE = part1.Energy + part2.Energy;
        Double_t tPt = tPx * tPx + tPy * tPy;
        Double_t tMt = tE * tE - tPz * tPz;  // mCVK;
        tMt = TMath::Sqrt(tMt);
        pair.Kt = TMath::Sqrt(tPt);
        Double_t tBeta  = tPz / tE;
        Double_t tGamma = tE / tMt;

        // Transform to LCMS

        Double_t particle1lcms_pz = tGamma * (part1.Pz - tBeta * part1.Energy);
        Double_t particle1lcms_e  = tGamma * (part1.Energy - tBeta * part1.Pz);
        Double_t particle2lcms_pz = tGamma * (part2.Pz - tBeta * part2.Energy);
        Double_t particle2lcms_e  = tGamma * (part2.Energy - tBeta * part2.Pz);

        // Rotate in transverse plane

        Double_t particle1lcms_px = (part1.Px * tPx + part1.Py * tPy) / pair.Kt;
        Double_t particle1lcms_py = (-part1.Px * tPy + part1.Py * tPx) / pair.Kt;

        Double_t particle2lcms_px = (part2.Px * tPx + part2.Py * tPy) / pair.Kt;
        Double_t particle2lcms_py = (-part2.Px * tPy + part2.Py * tPx) / pair.Kt;

        pair.QOut = abs(particle1lcms_px - particle2lcms_px);
        pair.QSide = abs(particle1lcms_py - particle2lcms_py);
        pair.QLong = abs(particle1lcms_pz - particle2lcms_pz);
        Double_t mDE = particle1lcms_e - particle2lcms_e;
        pair.QInv = TMath::Sqrt(TMath::Abs(pair.QOut * pair.QOut + pair.QSide * pair.QSide + pair.QLong * pair.QLong - mDE * mDE));

        return pair;
    }

    float CalcOpeningAngle(const TrackCandidate &part1, const TrackCandidate &part2)
    {
        float ptot = sqrt((part1.Px*part1.Px + part1.Py*part1.Py + part1.Pz*part1.Pz) * (part2.Px*part2.Px + part2.Py*part2.Py + part2.Pz*part2.Pz));
        if (ptot <= 0.)
        {	
            return 0.;
        }
        else
        {
            float arg = (part1.Px*part2.Px + part1.Py*part2.Py + part1.Pz*part2.Pz) / ptot;
            if (arg > 1.) return 1.;
            if (arg < -1.) return -1.;
            return acos(arg);
        }
    }

    PairCandidate CreatePair(const TrackCandidate &part1, const TrackCandidate &part2)
    {
        PairCandidate pairCand;
        pairCand = CFKinematics(part1,part2);
        pairCand.OpeningAngle = CalcOpeningAngle(part1,part2);
        pairCand.Rapidity = (part1.Rapidity + part2.Rapidity) / 2.;
        pairCand.AzimuthalAngle = (part1.AzimimuthalAngle + part2.AzimimuthalAngle) / 2.;
        pairCand.DeltaPhi = part1.AzimimuthalAngle - part2.AzimimuthalAngle;
        pairCand.DeltaTheta = part1.PolarAngle - part2.PolarAngle;

        return pairCand;
    }

    bool SelectPair(PairCandidate &pairCand, const float openAngle = 0.0)
    {
        std::vector<PairCandidate> acceptedPairs;

        if (pairCand.OpeningAngle <= openAngle)
            return false;

        return true;
    }

    void CalcSignal(std::vector<PairCandidate> &pairVec, const EventCandidate &eventCand, bool reorder, float openAngle = 0.0)
    {
        PairCandidate pairCand;
        for (size_t iter1 = 0; iter1 < eventCand.trackList.size(); iter1++)
            for (size_t iter2 = iter1+1; iter2 < eventCand.trackList.size(); iter2++)
            {
                if (reorder)
                    pairCand = CreatePair(eventCand.trackList.at(iter2),eventCand.trackList.at(iter1));
                else
                    pairCand = CreatePair(eventCand.trackList.at(iter1),eventCand.trackList.at(iter2));
                if (SelectPair(pairCand,openAngle))
                    pairVec.push_back(pairCand);
            }		
    }

    void CalcBackground(std::vector<PairCandidate> &pairVec, const EventCandidate &eventCand, const std::deque<TrackCandidate> &bckgVec, float openAngle = 0.0)
    {
        PairCandidate pairCand;
        for (auto &track : eventCand.trackList)
            for (auto &bckg : bckgVec)
            {
                pairCand = CreatePair(track,bckg);
                if (SelectPair(pairCand,openAngle))
                    pairVec.push_back(pairCand);
            }
                    
    }

    void MixBckg(const EventCandidate &eventCand, std::deque<TrackCandidate> &bckgVec, const int &maxSize, std::mt19937 &generator)
    {
        std::uniform_int_distribution<int> dist(0,eventCand.trackList.size()-1);
        bckgVec.push_back(eventCand.trackList.at(dist(generator)));
        if (bckgVec.size() > maxSize)
            bckgVec.pop_front();
    }

    int AssignKt(float kt)
    {
        static std::vector<std::pair<float,float> > ktBinVec = {{150,450},{450,750},{750,1050},{1050,1350},{1350,1650}};
        for (std::size_t i = 0; i < ktBinVec.size(); i++)
            if (kt >= ktBinVec.at(i).first && kt <= ktBinVec.at(i).second)
                return i;

        return -999;
    }

    float DegToRad(const int &angle)
    {
        return (TMath::Pi()/180)*angle;
    }

    bool IntToBool(const int &integ)
    {
        return integ;
    }

}

#endif