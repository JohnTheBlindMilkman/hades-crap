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
                fSettings(details::get<details::MixingOption::azimuthally_differential>(option::AzimuthallyDifferential{std::vector<float>{0.}},std::forward<Args>(args)...),
                details::get<details::MixingOption::centrality_differential>(option::CentralityDifferential{std::vector<int>{1}},std::forward<Args>(args)...),
                details::get<details::MixingOption::kT_differential>(option::KtDifferential{std::vector<float>{0.}},std::forward<Args>(args)...),
                details::get<details::MixingOption::rapidity_differential>(option::RapidityDifferential{std::vector<float>{0.}},std::forward<Args>(args)...),
                details::get<details::MixingOption::events_to_mix>(option::EventsToMix{50},std::forward<Args>(args)...)),
                fTargetPlates(15),
                fCentralityBins(get_value<details::MixingOption::centrality_differential>().size() - 1),
                fRapidityBins(get_value<details::MixingOption::rapidity_differential>().size() - 1),
                fKtBins(get_value<details::MixingOption::kT_differential>().size() - 1),
                fAzimuthBins(get_value<details::MixingOption::azimuthally_differential>().size() - 1),
                fIsCentralityDependent((get_value<details::MixingOption::centrality_differential>().size() > 1) ? true : false),
                fIsRapidityDependent((get_value<details::MixingOption::rapidity_differential>().size() > 1) ? true : false),
                fIsKtDependent((get_value<details::MixingOption::kT_differential>().size() > 1) ? true : false),
                fIsAzimuthallyDependent((get_value<details::MixingOption::azimuthally_differential>().size() > 1) ? true : false),
                InitMatrix(fTargetPlates,fCentralityBins,fRapidityBins,fKtBins,fAzimuthBins){}

            void PrintSettings()
            {
                std::cout << "Azimuthal bins:\n";
                for (const auto var : get_value<details::MixingOption::azimuthally_differential>())
                    std::cout << var <<"\t";
                
                std::cout << "\n Centrality bins:\n";
                for (const auto var : get_value<details::MixingOption::centrality_differential>())
                    std::cout << var <<"\t";

                std::cout << "\n kT bins:\n";
                for (const auto var : get_value<details::MixingOption::kT_differential>())
                    std::cout << var <<"\t";

                std::cout << "\n Rapidity bins:\n";
                for (const auto var : get_value<details::MixingOption::rapidity_differential>())
                    std::cout << var <<"\t";

                std::cout << "\n Events to mix:\t" << get_value<details::MixingOption::events_to_mix>() << "\n";
            }
            void MixAndDivide(const Selection::EventCandidate &event)
            {
                //get track vector
                //transform to pairs
                //remove unwanted
                //divide
                //return what???

                std::size_t centIter = 0, plateIter = 0;

                if (fIsCentralityDependent && 
                (get_value<details::MixingOption::centrality_differential>().front() <= event.Centrality &&
                get_value<details::MixingOption::centrality_differential>().back() > event.Centrality))
                {
                    for(std::size_t ii = 0; ii < get_value<details::MixingOption::centrality_differential>().size() - 1; ++ii)
                    {
                        if (get_value<details::MixingOption::centrality_differential>().at(ii) <= event.Centrality 
                        && get_value<details::MixingOption::centrality_differential>().at(ii+1) > event.Centrality)
                            centIter = ii;
                    }
                }

                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> randTrack(0,event.trackList.size()-1);

                std::shuffle(std::begin(event.trackList),std::end(event.trackList),gen); // tracks need to be randomised, becasue the sorter biased them
                std::vector<Selection::PairCandidate> signPairs = CalculateCombinations(event.trackList);
                signPairs.erase(std::remove_if(signPairs.begin(),signPairs.end(),[](const Selection::PairCandidate &pair){return pair.OpeningAngle > 0.0;}),signPairs.end());

                fSignal = boost::multi_array<std::vector<Selection::PairCandidate>,5>(boost::extents[fTargetPlates][fCentralityBins][fRapidityBins][fKtBins][fAzimuthBins]);
                DivideTracks(fSignal,plateIter,centIter,signPairs);

                fBckgCocktail.push_back(event.trackList.at(randTrack(gen)));
                if (fBckgCocktail.size() > get_value<details::MixingOption::events_to_mix>())
                    fBckgCocktail.pop_front();

                std::vector<Selection::PairCandidate> bckgPairs = CalculateCombinations(event.trackList,fBckgCocktail);

                fBackground = boost::multi_array<std::vector<Selection::PairCandidate>,5>(boost::extents[fTargetPlates][fCentralityBins][fRapidityBins][fKtBins][fAzimuthBins]);
                DivideTracks(fBackground,plateIter,centIter,bckgPairs);
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

            void InitMatrix(const std::size_t plateSize, const std::size_t centSize, const std::size_t rapSize, const std::size_t ktSize, const std::size_t azimuthSize)
            {
                fSignal = boost::multi_array<std::vector<Selection::PairCandidate>,5>(boost::extents[plateSize][centSize][rapSize][ktSize][azimuthSize]);
                fBackground = boost::multi_array<std::deque<Selection::PairCandidate>,5>(boost::extents[plateSize][centSize][rapSize][ktSize][azimuthSize]);
            }

            std::vector<Selection::PairCandidate> CalculateCombinations(const std::vector<Selection::TrackCandidate> &trackVec)
            {
                std::vector<Selection::PairCandidate> otpVec;

                for(std::size_t iter1 = 0; iter1 < trackVec.size(); ++iter1)
                    for(std::size_t iter2 = iter1+1; iter2 < trackVec.size(); ++iter2)
                        otpVec.push_back(Selection::CreatePair(trackVec.at(iter1),trackVec.at(iter2)));

            }

            std::vector<Selection::PairCandidate> CalculateCombinations(const std::vector<Selection::TrackCandidate> &trackVec,const std::deque<Selection::TrackCandidate> &backgVec)
            {
                std::vector<Selection::PairCandidate> otpVec;

                for(const auto &elem1 : trackVec)
                    for(const auto &elem2 : backgVec)
                        otpVec.push_back(Selection::CreatePair(elem1,elem2));

            }

            void DivideTracks(boost::multi_array<std::vector<Selection::PairCandidate>,5> &pairMultiArray,std::size_t plateIter,std::size_t centralityIter, const std::vector<Selection::PairCandidate> &pairVec)
            {
                std::size_t rapIter = 0,ktIter = 0,azimuthIter = 0;
                for (auto &pair : pairVec)
                {
                    rapIter = 0;
                    ktIter = 0;
                    azimuthIter = 0;

                    if (fIsRapidityDependent && 
                    (get_value<details::MixingOption::rapidity_differential>().front() <= pair.Rapidity 
                    && get_value<details::MixingOption::rapidity_differential>().back() > pair.Rapidity))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::rapidity_differential>().size() - 1; ++ii)
                        {
                            if (get_value<details::MixingOption::rapidity_differential>().at(ii) <= pair.Rapidity 
                            && get_value<details::MixingOption::rapidity_differential>().at(ii+1) > pair.Rapidity)
                                rapIter = ii;
                        }
                    }

                    if (fIsKtDependent && 
                    (get_value<details::MixingOption::kT_differential>().front() <= pair.Kt 
                    && get_value<details::MixingOption::kT_differential>().back() > pair.Kt))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::kT_differential>().size() - 1; ++ii)
                        {
                            if (get_value<details::MixingOption::kT_differential>().at(ii) <= pair.Kt 
                            && get_value<details::MixingOption::kT_differential>().at(ii+1) > pair.Kt)
                                ktIter = ii;
                        }
                    }

                    if (fIsAzimuthallyDependent && 
                    (get_value<details::MixingOption::azimuthally_differential>().front() <= pair.AzimuthalAngle 
                    && get_value<details::MixingOption::azimuthally_differential>().back() > pair.AzimuthalAngle))
                    {
                        for(std::size_t ii = 0; ii < get_value<details::MixingOption::azimuthally_differential>().size() - 1; ++ii)
                        {
                            if (get_value<details::MixingOption::azimuthally_differential>().at(ii) <= pair.AzimuthalAngle 
                            && get_value<details::MixingOption::azimuthally_differential>().at(ii+1) > pair.AzimuthalAngle)
                                azimuthIter = ii;
                        }
                    }

                    pairMultiArray[plateIter][centralityIter][rapIter][ktIter][azimuthIter].push_back(pair);
                }
            }

            Settings fSettings;
            const std::size_t fTargetPlates,fCentralityBins,fRapidityBins,fKtBins,fAzimuthBins;
            bool fIsCentralityDependent,fIsRapidityDependent,fIsKtDependent,fIsAzimuthallyDependent;
            boost::multi_array<std::vector<Selection::PairCandidate>,5> fSignal,fBackground;
            std::deque<Selection::TrackCandidate> fBckgCocktail;

    };
} // namespace FemtoCorrelation


