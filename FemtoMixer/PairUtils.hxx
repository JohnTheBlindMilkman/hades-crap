/**
 * @file PairUtils.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Collection of classes for consistent use of pair grouping and cuts
 * @version 0.1.0
 * @date 2025-10-30
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef PairUtils_hxx
    #define PairUtils_hxx

    #include "JJUtils.hxx"
    #include "PairCandidate.hxx"

    #include <array>

    namespace Mixing
    {
        /**
         * @brief Class storing intervals and observables according to which the pairs are grouped
         * 
         */
        class PairGrouping
        {
            private:
                /**
                 * @brief Creates an index sequence starting from 1
                 * 
                 * @tparam Is Input index sequence of 0, 1, 2, ...
                 * @return Index sequence of 1, 2, 3, ...
                 */
                template <std::size_t ... Is>
                [[nodiscard]] constexpr auto MakeIndexSequence(std::index_sequence<Is...>) const noexcept 
                {
                    return std::array<std::size_t,sizeof...(Is)>{(Is + 1)...};
                }
                /**
                 * @brief Removes all zeros created when converting floating point variable to a string. This function will also remove dot if it is the last character of the string
                 * 
                 * @param val String value of a unmoidified number
                 * @return The same number, but without zeros at the end
                 */
                TString RemoveTrailingZeros(TString val) const noexcept
                {
                    if (val.Contains('.'))
                    {
                        val.Remove(TString::EStripType::kTrailing,'0');
                        if (val.EndsWith("."))
                        {
                            val.Resize(val.Length() - 1);
                        }
                    }
                    
                    return val;
                }

                static constexpr std::size_t m_ktIntervals1D = 12, m_rapIntervals1D = 9;
                static constexpr std::array<float, m_ktIntervals1D + 1> m_ktArr1D = {300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,2100};
                static constexpr std::array<float, m_rapIntervals1D + 1> m_rapArr1D = {0.09,0.19,0.29,0.39,0.49,0.59,0.69,0.79,0.89,0.99};

                static constexpr std::size_t m_ktIntervals3D = 12, m_rapIntervals3D = 9, m_psiIntervals3D = 8;
                static constexpr std::array<float, m_ktIntervals3D + 1> m_ktArr3D = {300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,2100};
                static constexpr std::array<float, m_rapIntervals3D + 1> m_rapArr3D = {0.09,0.39,0.49,0.59,0.69,0.79,0.89,0.99};
                static constexpr std::array<float,m_psiIntervals3D + 1> m_epArr3D = {-202.5,-157.5,-112.5,-67.5,-22.5,22.5,67.5,112.5,157.5};

            public:
                /**
                 * @brief Calculates and returns group ID to which the PairCandidate object is assigned to. Used for 1D analysis
                 * 
                 * @param pair pointer to the PairCandidate
                 * @param ktArr array of kT intervals
                 * @param rapArr array of rapidity intervals
                 * @return group ID
                 */
                [[nodiscard]] std::string GetPairIndex1D(const std::shared_ptr<Selection::PairCandidate> &pair, const std::array<float, m_ktIntervals1D + 1> &ktArr, const std::array<float, m_rapIntervals1D + 1> &rapArr) const
                {
                    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair->GetKt()) - ktArr.begin();
                    std::size_t yCut = std::lower_bound(rapArr.begin(),rapArr.end(),pair->GetRapidity()) - rapArr.begin();

                    // reject if value is below first slice or above the last
                    if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > rapArr.size()-1)
                        return "0";
                    else
                        return JJUtils::to_fixed_size_string((ktCut),2) + JJUtils::to_fixed_size_string(yCut,2);
                }
                /**
                 * @brief Calculates and returns group ID to which the PairCandidate object is assigned to. Used for 1D analysis
                 * 
                 * @param pair pointer to the PairCandidate
                 * @param ktArr array of kT intervals
                 * @param rapArr array of rapidity intervals
                 * @param psiArr array of azimuthal angle intervals
                 * @return group ID
                 */
                [[nodiscard]] std::string GetPairIndex3D(const std::shared_ptr<Selection::PairCandidate> &pair, const std::array<float, m_ktIntervals3D + 1> &ktArr, const std::array<float, m_rapIntervals3D + 1> &rapArr, const std::array<float,m_psiIntervals3D + 1> &psiArr) const
                {
                    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair->GetKt()) - ktArr.begin();
                    std::size_t yCut = std::lower_bound(rapArr.begin(),rapArr.end(),pair->GetRapidity()) - rapArr.begin();
                    std::size_t EpCut = std::lower_bound(psiArr.begin(),psiArr.end(),pair->GetPhi()) - psiArr.begin();

                    // reject if value is below first slice or above the last
                    if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > rapArr.size()-1 || EpCut == 0 || EpCut > psiArr.size() - 1)
                        return "0";
                    else
                        return JJUtils::to_fixed_size_string(ktCut,2) + JJUtils::to_fixed_size_string(yCut,2) + JJUtils::to_fixed_size_string(EpCut,2);
                }
                /**
                 * @brief Creates a wrapper function for GetPairIndex1D
                 * 
                 * @return std::function
                 */
                [[nodiscard]] std::function<std::string (const std::shared_ptr<Selection::PairCandidate> &)> MakePairGroupingFunction1D() const noexcept
                {
                    auto newKtArr = m_ktArr1D; // this is a workaround, because I have the arrays marked as static and the std::function (i think) tries to move them (which is a big no-no according to the compiler)
                    auto newRapArr = m_rapArr1D;
                    return [this,newKtArr,newRapArr](const std::shared_ptr<Selection::PairCandidate> &pair) -> std::string {return this->GetPairIndex1D(pair,newKtArr,newRapArr);};
                }
                /**
                 * @brief Creates a wrapper function for GetPairIndex3D
                 * 
                 * @return std::function
                 */
                [[nodiscard]] std::function<std::string (const std::shared_ptr<Selection::PairCandidate> &)> MakePairGroupingFunction3D() const noexcept
                {
                    auto newKtArr = m_ktArr3D; // this is a workaround, because I have the arrays marked as static and the std::function (i think) tries to move them (which is a big no-no according to the compiler)
                    auto newRapArr = m_rapArr3D;
                    auto newPsiArr = m_epArr3D;
                    return [this,newKtArr,newRapArr,newPsiArr](const std::shared_ptr<Selection::PairCandidate> &pair) -> std::string {return this->GetPairIndex3D(pair,newKtArr,newRapArr,newPsiArr);};
                }
                /**
                 * @brief Get the kT index sequence for 1D anaysis
                 * 
                 * @return collection of kT indexes of 1, 2, 3, etc.
                 */
                [[nodiscard]] constexpr auto GetKtIndexSequence1D() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_ktIntervals1D>{});
                }
                /**
                 * @brief Get the rapidity index sequence for 1D anaysis
                 * 
                 * @return collection of rapidity indexes of 1, 2, 3, etc.
                 */
                [[nodiscard]] constexpr auto GetRapIndexSequence1D() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_rapIntervals1D>{});
                }
                /**
                 * @brief Get the kT index sequence for 3D anaysis
                 * 
                 * @return collection of kT indexes of 1, 2, 3, etc.
                 */
                [[nodiscard]] constexpr auto GetKtIndexSequence3D() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_ktIntervals3D>{});
                }
                /**
                 * @brief Get the rapidity index sequence for 3D anaysis
                 * 
                 * @return collection of rapidity indexes of 1, 2, 3, etc.
                 */
                [[nodiscard]] constexpr auto GetRapIndexSequence3D() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_rapIntervals3D>{});
                }
                /**
                 * @brief Get the azimuthal angle index sequence for 3D anaysis
                 * 
                 * @return collection of azimuthal angle indexes of 1, 2, 3, etc.
                 */
                [[nodiscard]] constexpr auto GetPsiIndexSequence3D() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_psiIntervals3D>{});
                }
                /**
                 * @brief Get the kT index sequence together with string containing the interval range used for 1D analysis. Useful for making legends
                 * 
                 * @return collection of pairs of indexes with strings
                 */
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_ktIntervals1D> GetKtIndexIntervalPairs1D() const noexcept 
                {
                    auto ktIndices = GetKtIndexSequence1D();
                    std::array<std::pair<std::size_t,TString>,m_ktIntervals1D> ktIndecesAndIntervals;
                    std::transform(ktIndices.begin(),ktIndices.end(),ktIndecesAndIntervals.begin(),
                        [this,newKtArr = m_ktArr1D](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i))));
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return ktIndecesAndIntervals;
                }
                /**
                 * @brief Get the rapidity index sequence together with string containing the interval range used for 1D analysis. Useful for making legends
                 * 
                 * @return collection of pairs of indexes with strings
                 */
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_rapIntervals1D> GetRapIndexIntervalPairs1D() const noexcept 
                {
                    auto rapIndices = GetRapIndexSequence1D();
                    std::array<std::pair<std::size_t,TString>,m_rapIntervals1D> rapIndecesAndIntervals;
                    std::transform(rapIndices.begin(),rapIndices.end(),rapIndecesAndIntervals.begin(),
                        [this,newRapArr = m_rapArr1D](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i))));
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return rapIndecesAndIntervals;
                }
                /**
                 * @brief Get the kT index sequence together with string containing the interval range used for 3D analysis. Useful for making legends
                 * 
                 * @return collection of pairs of indexes with strings
                 */
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_ktIntervals3D> GetKtIndexIntervalPairs3D() const noexcept 
                {
                    auto ktIndices = GetKtIndexSequence1D();
                    std::array<std::pair<std::size_t,TString>,m_ktIntervals3D> ktIndecesAndIntervals;
                    std::transform(ktIndices.begin(),ktIndices.end(),ktIndecesAndIntervals.begin(),
                        [this,newKtArr = m_ktArr3D](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i))));
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return ktIndecesAndIntervals;
                }
                /**
                 * @brief Get the rapidity index sequence together with string containing the interval range used for 3D analysis. Useful for making legends
                 * 
                 * @return collection of pairs of indexes with strings
                 */
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_rapIntervals3D> GetRapIndexIntervalPairs3D() const noexcept 
                {
                    auto rapIndices = GetRapIndexSequence1D();
                    std::array<std::pair<std::size_t,TString>,m_rapIntervals1D> rapIndecesAndIntervals;
                    std::transform(rapIndices.begin(),rapIndices.end(),rapIndecesAndIntervals.begin(),
                        [this,newRapArr = m_rapArr3D](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i))));
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return rapIndecesAndIntervals;
                }
                /**
                 * @brief Get the azimuthal angle index sequence together with string containing the interval range used for 3D analysis. Useful for making legends
                 * 
                 * @return collection of pairs of indexes with strings
                 */
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_psiIntervals3D> GetPsiIndexIntervalPairs3D() const noexcept 
                {
                    auto psiIndices = GetPsiIndexSequence3D();
                    std::array<std::pair<std::size_t,TString>,m_psiIntervals3D> psiIndecesAndIntervals;
                    std::transform(psiIndices.begin(),psiIndices.end(),psiIndecesAndIntervals.begin(),
                        [this,newPsiArr = m_epArr3D](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newPsiArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newPsiArr.at(i))));
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return psiIndecesAndIntervals;
                }
        };

        /**
         * @brief Class storing the pair rejection function used in my analysis
         * 
         */
        class PairRejection
        {
            private:
            
            public:
                /**
                 * @brief Pair rejection function. Evaluates wheather a given PairCandidate should be removed from the analysis
                 * 
                 * @param pair pointer to PairCandidate
                 * @return true if pair should be removed
                 * @return false otherwise
                 */
                [[nodiscard]] bool Reject(const std::shared_ptr<Selection::PairCandidate> &pair) const noexcept
                {
                    using Behaviour = Selection::PairCandidate::Behaviour;

                    if (pair->AreTracksFromTheSameSector())
                    {
                        return pair->RejectPairByCloseHits<Behaviour::OneUnder>(0.75,3) ||
                            pair->GetBothLayers() < 20 ||
                            pair->GetSharedMetaCells() > 0;
                    }
                    else
                    {
                        return false;
                    }
                }
                /**
                 * @brief Creates a wrapper function for Reject
                 * 
                 * @return std::function
                 */
                [[nodiscard]] std::function<bool (const std::shared_ptr<Selection::PairCandidate> &)> MakePairRejectionFunction() const noexcept
                {
                    return [this](const std::shared_ptr<Selection::PairCandidate> &pair){return this->Reject(pair);};
                }
        };
    }

#endif