#ifndef MdcWires_hxx
    #define MdcWires_hxx

    #include <array>
    #include <vector>
    #include <numeric>

    #include "hparticlemetamatcher.h"

    namespace HADES
    {
        namespace MDC
        {
            // some global constants which might be helpful with the wire list
            namespace WireInfo
            {
                constexpr short unsigned numberOfInnerLayers = 12;
                constexpr short unsigned numberOfOuterLayers = 12;
                constexpr short unsigned numberOfLayersInPlane = 6;
                constexpr short unsigned numberOfPlanes = 4;
                constexpr short unsigned numberOfAllLayers = numberOfInnerLayers + numberOfOuterLayers;
                constexpr short noWire = -1;
                constexpr std::array<short unsigned,numberOfPlanes> planeIndexing{0,1,2,3};
                constexpr std::array<short unsigned,numberOfInnerLayers> halfLayerIndexing{0,1,2,3,4,5,6,7,8,9,10,11};
                constexpr std::array<short unsigned,numberOfAllLayers> allLayerIndexing{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
                constexpr std::array<std::array<short unsigned,numberOfLayersInPlane>,numberOfPlanes> allLayerPerPlaneIndexing
                {
                    std::array<short unsigned,numberOfLayersInPlane>{0,1,2,3,4,5},
                    std::array<short unsigned,numberOfLayersInPlane>{6,7,8,9,10,11},
                    std::array<short unsigned,numberOfLayersInPlane>{12,13,14,15,16,17},
                    std::array<short unsigned,numberOfLayersInPlane>{18,19,20,21,22,23}
                };
            }

            /**
             * @brief Class representing all MDC layers. Holds an array of objects of type T for each MDC layer.
             * 
             * @tparam T object type which held at each layer
             */
            template <typename T>
            class MDCLayers
            {
                private:
                    std::array<T,WireInfo::numberOfAllLayers> m_layers;

                public:
                    using value_type = T;

                    constexpr auto begin() {return m_layers.begin();}
                    constexpr auto end() {return m_layers.end();}
                    constexpr const auto begin() const {return m_layers.begin();}
                    constexpr const auto end() const {return m_layers.end();}
                    [[nodiscard]] constexpr T& at(std::size_t index) {return m_layers.at(index);}
                    [[nodiscard]] constexpr const T& at(std::size_t index) const {return m_layers.at(index);}
                    [[nodiscard]] constexpr T& operator[](std::size_t index) noexcept {return m_layers[index];}
                    [[nodiscard]] constexpr const T& operator[](std::size_t index) const noexcept {return m_layers[index];}
            };
            /**
             * @brief Simple struct mimicing the behaviour of std::optional. If the default constructor is invoked, the struct is "empty", but if the main constructor is used it has some value.
             * 
             * @tparam T underlying object
             */
            template <typename T>
            struct OptionalDistance
            {
                T value;
                bool has_value;
                OptionalDistance() : value(T()), has_value(false) {}
                OptionalDistance(T val) : value(val), has_value(true) {}
                bool operator<(const OptionalDistance<T> &other) const noexcept 
                {
                    return (has_value && other.has_value) ? value < other.value : false;
                }
            };
            
            using LayersTrack = MDCLayers<std::vector<unsigned> >;
            using LayersPair = MDCLayers<std::pair<std::vector<unsigned>, std::vector<unsigned> > >;
            using WireDistances = MDCLayers<OptionalDistance<unsigned> >;

            /**
             * @brief Create a Track Layers object
             * 
             * @param wi HParticleWireInfo object, obtained from HParticleMetaMatcher
             * @return collection of wires fired by the track
             */
            [[nodiscard]] LayersTrack CreateTrackLayers(const HParticleWireInfo &wi)
            {
                LayersTrack array;

                for (std::size_t mod = 0; mod < HADES::MDC::WireInfo::numberOfPlanes; ++mod)
                    for (std::size_t lay = 0; lay < HADES::MDC::WireInfo::numberOfLayersInPlane; ++lay)
                    {
                        std::vector<unsigned> tmpVec;

                        if (wi.ar[mod][lay][0] > WireInfo::noWire) // safer than !=, because IDK if it won't be e.g. -2, -3, etc.
                            tmpVec.push_back(static_cast<unsigned>(wi.ar[mod][lay][0]));
                        if (wi.ar[mod][lay][1] > WireInfo::noWire)
                            tmpVec.push_back(static_cast<unsigned>(wi.ar[mod][lay][1]));

                        array[mod * HADES::MDC::WireInfo::numberOfLayersInPlane + lay] = tmpVec;
                    }

                return array;
            }
            /**
             * @brief Create a Pair Layers object
             * 
             * @param track1 collection of fired wires of the first track
             * @param track2 collection of fired wires of the second track
             * @return collection of pairs of fired wires
             */
            [[nodiscard]] LayersPair CreatePairLayers(const LayersTrack &track1, const LayersTrack &track2)
            {
                LayersPair pairLayers;
                std::transform(track1.begin(), track1.end(), track2.begin(), pairLayers.begin(),
                    [](const auto &layer1,const auto &layer2){return std::make_pair(layer1,layer2);});

                return pairLayers;
            }
            /**
             * @brief Calculate the number of MDC layers where both tracks had fired at least one wire
             * 
             * @param pairLayers collection of pairs of fired wires
             * @return  number of layers where both traks had fired a wire 
             */
            [[nodiscard]] unsigned CalculateBothLayers(const LayersPair &pairLayers)
            {
                MDCLayers<bool> bothLayers;
                std::transform(pairLayers.begin(), pairLayers.end(), bothLayers.begin(),
                    [](const auto &layer){return static_cast<bool>(!layer.first.empty() && !layer.second.empty());});

                return std::count_if(bothLayers.begin(),bothLayers.end(),[](bool areBoth){return areBoth;});
            }
            /**
             * @brief Calculate smallest distance between all the wires for each MDC layer
             * 
             * @param pairLayers collection of pairs of fired wires
             * @return collection of smallest distances 
             */
            [[nodiscard]] WireDistances CalculateWireDistances(const LayersPair &pairLayers)
            {
                WireDistances wireDistsances;

                auto uniaryOp = [](const std::pair<std::vector<unsigned>, std::vector<unsigned> > &layer)
                {
                    if (!layer.first.empty() && !layer.second.empty())
                    {
                        long minDist = std::numeric_limits<unsigned>::max();

                        for(const auto &wire1 : layer.first)
                            for (const auto &wire2 : layer.second)
                            {
                                // I cast to long instead of int to avoid narrowing conversion warning, even though I know it won't narrow it down
                                minDist = std::min(minDist, std::abs(static_cast<long>(wire1) - static_cast<long>(wire2)));
                            }

                        return OptionalDistance<unsigned>(minDist);
                    }
                    else
                    {
                        return OptionalDistance<unsigned>();
                    }
                };

                std::transform(pairLayers.begin(),pairLayers.end(),wireDistsances.begin(),uniaryOp);

                return wireDistsances;
            }
            /**
             * @brief Calculate number of shared wires within the pair of tracks
             * 
             * @param pairLayers collection of pairs of fired wires
             * @return number of shared wires
             */
            [[nodiscard]] unsigned CalculateSharedWires(const LayersPair &pairLayers)
            {
                MDCLayers<unsigned> sharedWires;
                auto uniaryOp = [](const std::pair<std::vector<unsigned>, std::vector<unsigned> > &layer)
                {
                    unsigned sw = 0;
                    for(const auto &wire1 : layer.first)
                        for (const auto &wire2 : layer.second)
                        {
                            if (wire1 == wire2)
                            {
                                ++sw;
                            }
                        }

                    return sw;
                };

                std::transform(pairLayers.begin(), pairLayers.end(), sharedWires.begin(), uniaryOp);

                return std::accumulate(sharedWires.begin(),sharedWires.end(),0);
            }
            /**
             * @brief Calculate pair splitting level according to Adams J., et al., Phys. Rev. C 71.4 (2005): 044906
             * 
             * @param pairLayers collection of pairs of fired wires
             * @return splitting level (e.g. SL = -0.5 no splitting, SL = 1 possible split tracks)
             */
            [[nodiscard]] double CalcluateSplittingLevel(const LayersPair &pairLayers)
            {
                unsigned nHits1 = 0, nHits2 = 0;
                MDCLayers<int> splittingLevels;
                auto uniaryOp = [&nHits1,&nHits2](const std::pair<std::vector<unsigned>, std::vector<unsigned> > &layer)
                {
                    int splittingLvl = 0;

                    if (!layer.first.empty() && !layer.second.empty())
                    {
                        --splittingLvl;
                        ++nHits1;
                        ++nHits2;
                    }
                    else if (!layer.first.empty())
                    {
                        ++splittingLvl;
                        ++nHits1;
                    }
                    else if (!layer.second.empty())
                    {
                        ++splittingLvl;
                        ++nHits2;
                    }

                    return splittingLvl;
                };

                std::transform(pairLayers.begin(), pairLayers.end(), splittingLevels.begin(), uniaryOp);

                return (nHits1 > 0 || nHits2 > 0) ? std::accumulate(splittingLevels.begin(),splittingLevels.end(),0.) / (nHits1 + nHits2) : 0.;
            }
        } // namespace MDC
        
    } // namespace HADES
    
#endif