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

                static constexpr std::size_t m_ktIntervals = 10, m_rapIntervals = 13, m_psiIntervals = 8;
                static constexpr std::array<float, m_ktIntervals + 1> m_ktArr = {0,200,400,600,800,1000,1200,1400,1600,1800,2000};
                static constexpr std::array<float, m_rapIntervals + 1> m_rapArr = {0.09,0.19,0.29,0.39,0.49,0.59,0.69,0.79,0.89,0.99,1.09,1.19,1.29,1.39};
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
                    std::array<float, m_ktIntervals + 1> newKtArr = m_ktArr;
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
                        return pair->RejectPairByCloseHits<Behaviour::OneUnder>(0.5,3) ||
                            pair->GetBothLayers() < 20 ||
                            pair->GetSharedMetaCells() > 0;
                    }
                    else
                    {
                        return false;
                    }
                }
                [[nodiscard]] std::function<bool (const std::shared_ptr<Selection::PairCandidate> &)> MakePairRejectionFunction() noexcept
                {
                    return [this](const std::shared_ptr<Selection::PairCandidate> &pair){return this->Reject(pair);};
                }
        };
    }

#endif