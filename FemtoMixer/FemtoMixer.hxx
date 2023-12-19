#ifndef FemtoMixer_hxx
    #define FemtoMixer_hxx

#include "Options.hxx"
#include "PairCandidate.hxx"
#include "boost/multi_array.hpp"
#include <iostream>
#include <queue>
#include <algorithm>
#include <random>

namespace FemtoCorrelation
{
    class FemtoMixer
    {
        using Settings = std::tuple<option::AzimuthallyDifferential,option::CentralityDifferential,option::KtDifferential,option::RapidityDifferential,option::EventsToMix>;
        
        public:
            enum class Distribution{Signal,Background};

            template <typename... Args, typename std::enable_if<details::are_settings_from_tuple<Settings, typename std::decay<Args>::type...>::value,void *>::type = nullptr>
            explicit FemtoMixer(Args &&... args) : 
                fSettings(details::get<details::MixingOption::azimuthally_differential>(option::AzimuthallyDifferential{std::vector<std::pair<float,float> >{}},std::forward<Args>(args)...),
                details::get<details::MixingOption::centrality_differential>(option::CentralityDifferential{std::vector<std::pair<int,int> >{}},std::forward<Args>(args)...),
                details::get<details::MixingOption::kT_differential>(option::KtDifferential{std::vector<std::pair<float,float> >{}},std::forward<Args>(args)...),
                details::get<details::MixingOption::rapidity_differential>(option::RapidityDifferential{std::vector<std::pair<float,float>>{}},std::forward<Args>(args)...),
                details::get<details::MixingOption::events_to_mix>(option::EventsToMix{50},std::forward<Args>(args)...)),
                fTargetPlates(15),
                fCentralityBins(get_value<details::MixingOption::centrality_differential>().size()),
                fRapidityBins(get_value<details::MixingOption::rapidity_differential>().size()),
                fKtBins(get_value<details::MixingOption::kT_differential>().size()),
                fAzimuthBins(get_value<details::MixingOption::azimuthally_differential>().size()),
                fIsCentralityDependent((get_value<details::MixingOption::centrality_differential>().size() > 0) ? true : false),
                fIsRapidityDependent((get_value<details::MixingOption::rapidity_differential>().size() > 0) ? true : false),
                fIsKtDependent((get_value<details::MixingOption::kT_differential>().size() > 0) ? true : false),
                fIsAzimuthallyDependent((get_value<details::MixingOption::azimuthally_differential>().size() > 0) ? true : false),
                fSignal(boost::extents[fTargetPlates][fCentralityBins][fRapidityBins][fKtBins][fAzimuthBins]),
                fBackground(boost::extents[fTargetPlates][fCentralityBins][fRapidityBins][fKtBins][fAzimuthBins])
                {}

            void PrintSettings()
            {
                std::cout << "\n";
                std::cout << "--------------------------------------------------------\n";
                std::cout << "FemtoMixer::PrintSettings()\n";
                std::cout << "--------------------------------------------------------\n";
                std::cout << "\n";

                std::cout << "Azimuthal bins:\n";
                for (const auto var : get_value<details::MixingOption::azimuthally_differential>())
                    std::cout << "[" << var.first <<", " << var.second <<")  ";
                
                std::cout << "\n Centrality bins:\n";
                for (const auto var : get_value<details::MixingOption::centrality_differential>())
                    std::cout << "[" << var.first <<", " << var.second <<")  ";

                std::cout << "\n kT bins:\n";
                for (const auto var : get_value<details::MixingOption::kT_differential>())
                    std::cout << "[" << var.first <<", " << var.second <<")  ";

                std::cout << "\n Rapidity bins:\n";
                for (const auto var : get_value<details::MixingOption::rapidity_differential>())
                    std::cout << "[" << var.first <<", " << var.second <<")  ";

                std::cout << "\n Events to mix:\t" << get_value<details::MixingOption::events_to_mix>() << "\n";

                std::cout << "--------------------------------------------------------\n";
                std::cout << "\n";
            }
            void MixAndDivide(Selection::EventCandidate &event)
            {
                std::size_t centIter = 0, plateIter = event.TargetPlate;

                if (fIsCentralityDependent && 
                (get_value<details::MixingOption::centrality_differential>().front().first <= event.Centrality &&
                get_value<details::MixingOption::centrality_differential>().back().second > event.Centrality))
                {
                    for(std::size_t ii = 0; ii < get_value<details::MixingOption::centrality_differential>().size(); ++ii)
                    {
                        if (get_value<details::MixingOption::centrality_differential>().at(ii).first <= event.Centrality 
                        && get_value<details::MixingOption::centrality_differential>().at(ii).second > event.Centrality)
                            centIter = ii;
                    }
                }

                ClearAndResizeVectors(fSignal);
                ClearAndResizeVectors(fBackground);

                std::mt19937 gen(fRandDevice());
                std::uniform_int_distribution<> randTrack(0,event.trackList.size()-1);

                std::shuffle(event.trackList.begin(),event.trackList.end(),gen); // tracks need to be randomised, becasue the sorter biased them

                std::vector<Selection::PairCandidate> signPairs = CalculateCombinations(event.trackList);
                if (signPairs.size() > 0)
                {
                    signPairs.erase(std::remove_if(signPairs.begin(),signPairs.end(),[](const Selection::PairCandidate &pair){return (pair.OpeningAngle > 0.0);}),signPairs.end());
                    DivideTracks(fSignal,plateIter,centIter,signPairs);
                }

                std::vector<Selection::PairCandidate> bckgPairs = CalculateCombinations(event.trackList,fBckgCocktail);
                if (bckgPairs.size() > 0)
                {
                    bckgPairs.erase(std::remove_if(bckgPairs.begin(),bckgPairs.end(),[](const Selection::PairCandidate &pair){return (pair.OpeningAngle > 0.0);}),bckgPairs.end());
                    DivideTracks(fBackground,plateIter,centIter,bckgPairs);
                }
                
                fBckgCocktail.push_back(event.trackList.at(randTrack(gen))); //bckg cocktail should have centrality and plate dependence
                if (fBckgCocktail.size() > get_value<details::MixingOption::events_to_mix>())
                    fBckgCocktail.pop_front();
            }

            std::vector<Selection::PairCandidate> GetPairsAtIndicies(std::size_t plate,std::size_t centrality,std::size_t rapidity, std::size_t kt, std::size_t azimuth,Distribution dist)
            {
                switch (dist)
                {
                    case Distribution::Signal:
                        return fSignal[plate][centrality][rapidity][kt][azimuth];
                    case Distribution::Background:
                        return fBackground[plate][centrality][rapidity][kt][azimuth];
                    default:
                        throw std::runtime_error("FemtoMixer::GetPairsAtIndicies - unknow Distribution option");
                }
            }

        private:
            template <details::MixingOption id>
            auto get_value()
                -> decltype((details::get_value<id>(std::declval<Settings &>()).value)) 
                {
                    return details::get_value<id>(fSettings).value;
                }

            template <details::MixingOption id>
            auto get_value() const 
                -> decltype((details::get_value<id>(std::declval<const Settings &>()).value)) 
                {
                    return details::get_value<id>(fSettings).value;
                }

            std::vector<Selection::PairCandidate> CalculateCombinations(const std::vector<Selection::TrackCandidate> &trackVec)
            {
                std::vector<Selection::PairCandidate> otpVec;

                for(std::size_t iter1 = 0; iter1 < trackVec.size(); ++iter1)
                    for(std::size_t iter2 = iter1+1; iter2 < trackVec.size(); ++iter2)
                        otpVec.push_back(Selection::CreatePair(trackVec.at(iter1),trackVec.at(iter2)));

                return otpVec;
            }

            std::vector<Selection::PairCandidate> CalculateCombinations(const std::vector<Selection::TrackCandidate> &trackVec,const std::deque<Selection::TrackCandidate> &backgVec)
            {
                std::vector<Selection::PairCandidate> otpVec;

                for(const auto &elem1 : trackVec)
                    for(const auto &elem2 : backgVec)
                        otpVec.push_back(Selection::CreatePair(elem1,elem2));

                return otpVec;
            }

            void DivideTracks(boost::multi_array<std::vector<Selection::PairCandidate>,5> &pairMultiArray,std::size_t plateIter,std::size_t centralityIter, const std::vector<Selection::PairCandidate> &pairVec)
            {
                std::size_t rapIter = 0,ktIter = 0,azimuthIter = 0;
                for (const auto &pair : pairVec)
                {
                    rapIter = 0;
                    ktIter = 0;
                    azimuthIter = 0;

                    if (fIsRapidityDependent && 
                    (get_value<details::MixingOption::rapidity_differential>().front().first <= pair.Rapidity 
                    && get_value<details::MixingOption::rapidity_differential>().back().second > pair.Rapidity))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::rapidity_differential>().size(); ++ii)
                        {
                            if (get_value<details::MixingOption::rapidity_differential>().at(ii).first <= pair.Rapidity 
                            && get_value<details::MixingOption::rapidity_differential>().at(ii).second > pair.Rapidity)
                                rapIter = ii;
                        }
                    }

                    if (fIsKtDependent && 
                    (get_value<details::MixingOption::kT_differential>().front().first <= pair.Kt 
                    && get_value<details::MixingOption::kT_differential>().back().second > pair.Kt))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::kT_differential>().size(); ++ii)
                        {
                            if (get_value<details::MixingOption::kT_differential>().at(ii).first <= pair.Kt 
                            && get_value<details::MixingOption::kT_differential>().at(ii).second > pair.Kt)
                                ktIter = ii;
                        }
                    }

                    if (fIsAzimuthallyDependent && 
                    (get_value<details::MixingOption::azimuthally_differential>().front().first <= pair.AzimuthalAngle 
                    && get_value<details::MixingOption::azimuthally_differential>().back().second > pair.AzimuthalAngle))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::azimuthally_differential>().size(); ++ii)
                        {
                            if (get_value<details::MixingOption::azimuthally_differential>().at(ii).first <= pair.AzimuthalAngle 
                            && get_value<details::MixingOption::azimuthally_differential>().at(ii).second > pair.AzimuthalAngle)
                                azimuthIter = ii;
                        }
                    }

                    pairMultiArray[plateIter][centralityIter][rapIter][ktIter][azimuthIter].push_back(pair);
                }
            }

            void ClearAndResizeVectors(boost::multi_array<std::vector<Selection::PairCandidate>,5> &pairMultiArray)
            {
                for (std::size_t plateIter = 0; plateIter < fTargetPlates; ++plateIter)
                    for (std::size_t centralityIter = 0; centralityIter < fCentralityBins; ++centralityIter)
                        for (std::size_t rapIter = 0; rapIter < fRapidityBins; ++rapIter)
                            for (std::size_t ktIter = 0; ktIter < fKtBins; ++ktIter)
                                for (std::size_t azimuthIter = 0; azimuthIter < fAzimuthBins; ++azimuthIter)
                                {
                                    pairMultiArray[plateIter][centralityIter][rapIter][ktIter][azimuthIter].clear();
                                    pairMultiArray[plateIter][centralityIter][rapIter][ktIter][azimuthIter].resize(0);
                                }
            }

            Settings fSettings;
            const std::size_t fTargetPlates,fCentralityBins,fRapidityBins,fKtBins,fAzimuthBins;
            bool fIsCentralityDependent,fIsRapidityDependent,fIsKtDependent,fIsAzimuthallyDependent;
            boost::multi_array<std::vector<Selection::PairCandidate>,5> fSignal,fBackground;
            std::deque<Selection::TrackCandidate> fBckgCocktail;
            std::random_device fRandDevice;

    };
} // namespace FemtoCorrelation

#endif
