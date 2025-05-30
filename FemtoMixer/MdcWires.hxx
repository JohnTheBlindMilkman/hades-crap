#ifndef MdcWires_hxx
    #define MdcWires_hxx

    #include <array>
    #include <vector>
    #include <numeric>

    #include "hparticlemetamatcher.h"

    #include "JJUtils.hxx"

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
             * @brief Struct representing an MDC Wire. Either stores a value (if the main constructor is used) or doesn't (if is default-constructed).  Behaves similarly as an std::optional<T>.
             * 
             * @tparam T stored type
             */
            template <typename T>
            struct Wire
            {
                T id;
                bool was_fired;
                Wire() : id(T()), was_fired(false) {}
                explicit Wire (T val) : id(val), was_fired(true) {}
                bool operator==(const Wire &other) const noexcept {return (id == other.id && was_fired && other.was_fired);}
            };

            /**
             * @brief Struct representing a pair of MDC Wires.
             * 
             * @tparam T stored type
             */
            template <typename T>
            struct WirePair
            {
                std::pair<Wire<T>, Wire<T> > pair;
                WirePair(Wire<T> &&wire1, Wire<T> &&wire2) : pair(std::make_pair(std::forward<Wire<T>>(wire1), std::forward<Wire<T>>(wire2))) {}
                T Distance(const WirePair<T> &other) const noexcept
                {
                    T dist1 = std::min(std::abs(pair.first.id - other.pair.first.id),std::abs(pair.second.id - other.pair.second.id));
                    T dist2 = std::min(std::abs(pair.first.id - other.pair.second.id),std::abs(pair.second.id - other.pair.first.id));
                }
                bool operator==(const WirePair &other) const noexcept
                {
                    if ((pair.first == other.pair.first && pair.second == other.pair.second) || 
                        (pair.first == other.pair.second && pair.second == other.pair.first));
                }
            };

            /**
             * @brief Struct representing an MDC layer. It either holds a value (if the main constructor is used) or is empty (if default constructor is used). Behaves similarly as std::optional<T>
             * 
             * @tparam T stored type
             */
            template <typename T>
            struct Layer
            {
                T value;
                bool is_empty;
                Layer() : value(T()), is_empty(true) {}
                explicit Layer(T val) : value(val), is_empty(false) {}
                bool operator==(const Layer &other) const noexcept {return (value == other.value && (!is_empty) && (!other.is_empty));}
            };

            /**
             * @brief Class representing all MDC layers. Holds an array of type Layer<T>.
             * 
             * @tparam T stored type
             */
            template <typename T>
            class Layers
            {
                private:
                    std::array<Layer<T>,WireInfo::numberOfAllLayers> m_layers;

                public:
                    Layer<T>& operator[](std::size_t index)
                    {
                        return m_layers[index];
                    }
                    const Layer<T>& operator[](std::size_t index) const
                    {
                        return m_layers[index];
                    }
                    Layer<T>& at(std::size_t index)
                    {
                        return m_layers.at(index);
                    }
                    const Layer<T>& at(std::size_t index) const
                    {
                        return m_layers.at(index);
                    }
                    auto begin()
                    {
                        return m_layers.begin();
                    }
                    const auto begin() const
                    {
                        return m_layers.begin();
                    }
                    auto end()
                    {
                        return m_layers.end();
                    }
                    const auto end() const
                    {
                        return m_layers.end();
                    }
            };

            // using TrackLayers = Layers<std::vector<unsigned> >;
            // using PairLayers = Layers<std::vector<std::pair<unsigned, unsigned> > >;
            // using PairDistances = Layers<unsigned>;

            class FiredWiresOps
            {
                /**
                 * @brief Create a modernised wire collection object
                 * 
                 * @param wi HParticleWireInfo object, obtained from HParticleMetaMatcher
                 * @return Collection of MDC layers
                 */
                [[nodiscard]] static Layers<WirePair<unsigned> > CreateWireArray(const HParticleWireInfo &wi)
                {
                    Layers<WirePair<unsigned> > array;

                    for (std::size_t mod = 0; mod < HADES::MDC::WireInfo::numberOfPlanes; ++mod)
                        for (std::size_t lay = 0; lay < HADES::MDC::WireInfo::numberOfLayersInPlane; ++lay)
                        {
                            int valWire1 = wi.ar[mod][lay][0];
                            int valWire2 = wi.ar[mod][lay][1];

                            array[mod * HADES::MDC::WireInfo::numberOfLayersInPlane + lay] = 
                                Layer<WirePair<unsigned> >(
                                    WirePair<unsigned>(
                                        ((valWire1 > -1) ? Wire<unsigned>(static_cast<unsigned>(valWire1)) : Wire<unsigned>()),
                                        ((valWire2 > -1) ? Wire<unsigned>(static_cast<unsigned>(valWire2)) : Wire<unsigned>())
                                    )
                                );
                        }

                    return array;
                }
                // [[nodiscard]] static PairLayers CreatePairLayers(const TrackLayers &layersTrack1, const TrackLayers &layersTrack2)
                // {
                //     PairLayers array;
                //     auto lambda = [](const Layer<std::vector<unsigned> > &layer1, const Layer<std::vector<unsigned> > &layer2)
                //     {
                //         if (layer1.has_value && layer2.has_value)
                //         {
                            

                //             return Layer<std::vector<std::pair<unsigned, unsigned> > >();
                //         }
                //         else
                //         {
                //             return Layer<std::vector<std::pair<unsigned, unsigned> > >();
                //         }
                //     };
                    
                // }
                /**
                 * @brief 
                 * 
                 * @param layersTrack1 
                 * @param layersTrack2 
                 * @return Layers<unsigned> 
                 */
                [[nodiscard]] static Layers<unsigned> CalculateWireDistances(const Layers<WirePair<unsigned> > &layersTrack1, const Layers<WirePair<unsigned> > &layersTrack2)
                {
                    Layers<unsigned> wireDistances;

                    std::transform(layersTrack1.begin(),layersTrack1.end(),layersTrack2.begin(),
                                    std::back_inserter(wireDistances),
                                    [](const auto &layer1, const auto &layer2)
                                    {
                                        if (layer1.has_value && layer2.has_value)
                                        {
                                            unsigned minDist = std::numeric_limits<unsigned>::max();

                                            for (const auto &bothWires : JJUtils::zip_to_bigger(layer1.value,layer2.value))
                                            {
                                                minDist = std::min(minDist,std::abs(bothWires.first - bothWires.second));
                                            }

                                            return Layer<unsgined>(minDist);
                                        }
                                        else
                                        {
                                            return Layer<unsigned>();
                                        }
                                    });

                    return wireDistances;
                }
                // [[nodiscard]] static double CalculateSplittingLevel(const Layers<std::vector<unsigned> > &layersTrack1, const Layers<std::vector<unsigned> > &layersTrack2)
                // {
                //     std::array<double,WireInfo::numberOfAllLayers> splittingVals;

                //     std::transform(layersTrack1.begin(),layersTrack1.end(),layersTrack2.begin(),
                //                     std::back_inserter(splittingVals),
                //                     [](const auto &layer1, const auto &layer2)
                //                     {
                //                         if (layer1.has_value && layer2.has_value)
                //                         {
                //                             double SL = 0;
                                            
                //                             // revisit the STAR paper

                //                             return SL;
                //                         }
                //                         else
                //                         {
                //                             return -1;
                //                         }
                //                     });

                //     //return std::accumulate(splittingVals.begin(),splittingVals.end(),0.) / ;
                // }
                // [[nodiscard]] static unsigned CalculateMinWireDistance()
                // {

                // }
            };

        } // namespace MDC
        
    } // namespace HADES
    
#endif