#ifndef PairUtils_hxx
    #define PairUtils_hxx

    #include "JJUtils.hxx"
    #include "PairCandidate.hxx"

    #include <array>

    namespace Mixing
    {
        class PairGrouping
        {
            private:
                template <std::size_t ... Is>
                [[nodiscard]] constexpr auto MakeIndexSequence(std::index_sequence<Is...>) const noexcept 
                {
                    return std::array<std::size_t,sizeof...(Is)>{(Is + 1)...};
                }

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

                static constexpr std::size_t m_ktIntervals = 12, m_rapIntervals = 9, m_psiIntervals = 8;
                static constexpr std::array<float, m_ktIntervals + 1> m_ktArr = {300,450,600,750,900,1050,1200,1350,1500,1650,1800,1950,2100};
                static constexpr std::array<float, m_rapIntervals + 1> m_rapArr = {0.09,0.19,0.29,0.39,0.49,0.59,0.69,0.79,0.89,0.99};
                static constexpr std::array<float,m_psiIntervals + 1> m_epArr{-202.5,-157.5,-112.5,-67.5,-22.5,22.5,67.5,112.5,157.5};

            public:
                [[nodiscard]] std::string GetPairIndex1D(const std::shared_ptr<Selection::PairCandidate> &pair, const std::array<float, m_ktIntervals + 1> &ktArr, const std::array<float, m_rapIntervals + 1> &rapArr) const
                {
                    std::size_t ktCut = std::lower_bound(ktArr.begin(),ktArr.end(),pair->GetKt()) - ktArr.begin();
                    std::size_t yCut = std::lower_bound(rapArr.begin(),rapArr.end(),pair->GetRapidity()) - rapArr.begin();

                    // reject if value is below first slice or above the last
                    if (ktCut == 0 || ktCut > ktArr.size()-1 || yCut == 0 || yCut > rapArr.size()-1)
                        return "0";
                    else
                        return JJUtils::to_fixed_size_string((ktCut),2) + JJUtils::to_fixed_size_string(yCut,2);
                }
                [[nodiscard]] std::string GetPairIndex3D(const std::shared_ptr<Selection::PairCandidate> &pair) const
                {
                    std::size_t ktCut = std::lower_bound(m_ktArr.begin(),m_ktArr.end(),pair->GetKt()) - m_ktArr.begin();
                    std::size_t yCut = std::lower_bound(m_rapArr.begin(),m_rapArr.end(),pair->GetRapidity()) - m_rapArr.begin();
                    std::size_t EpCut = std::lower_bound(m_epArr.begin(),m_epArr.end(),pair->GetPhi()) - m_epArr.begin();

                    // reject if value is below first slice or above the last
                    if (ktCut == 0 || ktCut > m_ktArr.size()-1 || yCut == 0 || yCut > m_rapArr.size()-1 || EpCut == 0 || EpCut > m_epArr.size() - 1)
                        return "0";
                    else
                        return JJUtils::to_fixed_size_string((ktCut),2) + JJUtils::to_fixed_size_string(yCut,2) + JJUtils::to_fixed_size_string(EpCut,2);
                }
                [[nodiscard]] std::function<std::string (const std::shared_ptr<Selection::PairCandidate> &)> MakePairGroupingFunction1D() const noexcept
                {
                    std::array<float, m_ktIntervals + 1> newKtArr = m_ktArr; // this is a workaround, because I have the arrays marked as static and the std::function (i think) tries to move them (which is a big no-no according to the compiler)
                    std::array<float, m_rapIntervals + 1> newRapArr = m_rapArr;
                    return [this,newKtArr,newRapArr](const std::shared_ptr<Selection::PairCandidate> &pair) -> std::string {return this->GetPairIndex1D(pair,newKtArr,newRapArr);};
                }
                [[nodiscard]] constexpr auto GetKtIndexSequence() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_ktIntervals>{});
                }
                [[nodiscard]] constexpr auto GetRapIndexSequence() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_rapIntervals>{});
                }
                [[nodiscard]] constexpr auto GetPsiIndexSequence() const noexcept 
                {
                    return MakeIndexSequence(std::make_index_sequence<m_psiIntervals>{});
                }
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_ktIntervals> GetKtIndexIntervalPairs() const noexcept 
                {
                    std::array<float, m_ktIntervals + 1> newKtArr = m_ktArr;
                    auto ktIndices = GetKtIndexSequence();
                    std::array<std::pair<std::size_t,TString>,m_ktIntervals> ktIndecesAndIntervals;
                    std::transform(ktIndices.begin(),ktIndices.end(),ktIndecesAndIntervals.begin(),
                        [this,newKtArr](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newKtArr.at(i))));
                            if (highEdge == RemoveTrailingZeros(TString(std::to_string(std::numeric_limits<float>::max()))))
                                highEdge = "#infty";
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return ktIndecesAndIntervals;
                }
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_rapIntervals> GetRapIndexIntervalPairs() const noexcept 
                {
                    std::array<float, m_rapIntervals + 1> newRapArr = m_rapArr;
                    auto rapIndices = GetRapIndexSequence();
                    std::array<std::pair<std::size_t,TString>,m_rapIntervals> rapIndecesAndIntervals;
                    std::transform(rapIndices.begin(),rapIndices.end(),rapIndecesAndIntervals.begin(),
                        [this,newRapArr](std::size_t i)
                        {
                            TString lowEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i - 1))));
                            TString highEdge = RemoveTrailingZeros(TString(std::to_string(newRapArr.at(i))));
                            if (highEdge == RemoveTrailingZeros(TString(std::to_string(std::numeric_limits<float>::max()))))
                                highEdge = "#infty";
                            return std::make_pair(i,TString::Format("(%s,%s)", lowEdge.Data(), highEdge.Data()));
                        }
                    );

                    return rapIndecesAndIntervals;
                }
                [[nodiscard]] std::array<std::pair<std::size_t,TString>,m_psiIntervals> GetPsiIndexIntervalPairs() const noexcept 
                {
                    std::array<float, m_psiIntervals + 1> newPsiArr = m_epArr;
                    auto psiIndices = GetPsiIndexSequence();
                    std::array<std::pair<std::size_t,TString>,m_psiIntervals> psiIndecesAndIntervals;
                    std::transform(psiIndices.begin(),psiIndices.end(),psiIndecesAndIntervals.begin(),
                        [newPsiArr](std::size_t i){return std::make_pair(i,TString::Format("(%f,%f)",newPsiArr.at(i - 1), newPsiArr.at(i)));}
                    );

                    return psiIndecesAndIntervals;
                }
        };

        class PairRejection
        {
            private:
            
            public:
                [[nodiscard]] bool Reject(const std::shared_ptr<Selection::PairCandidate> &pair) const noexcept
                {
                    using Behaviour = Selection::PairCandidate::Behaviour;

                    if (pair->AreTracksFromTheSameSector())
                    {
                        return pair->RejectPairByCloseHits<Behaviour::OneUnder>(0.75,3) ||
                            pair->GetBothLayers() < 21 ||
                            pair->GetSharedMetaCells() > 0;
                    }
                    else
                    {
                        return false;
                    }
                }
                [[nodiscard]] std::function<bool (const std::shared_ptr<Selection::PairCandidate> &)> MakePairRejectionFunction() const noexcept
                {
                    return [this](const std::shared_ptr<Selection::PairCandidate> &pair){return this->Reject(pair);};
                }
        };
    }

#endif