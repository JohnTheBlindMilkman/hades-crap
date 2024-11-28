/**
 * @file JJFemtoMixer.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Class responsible for mixing and storing pairs.
 * @version 1.0
 * @date 2024-01-15
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef JJFemtoMixer_hxx
    #define JJFemtoMixer_hxx

    #include <vector>
    #include <deque>
    #include <map>
    #include <functional>
    #include <random>
    #include <type_traits>

    namespace Mixing
    {
        template<typename Event, typename Track, typename Pair>
        class JJFemtoMixer
        {
            static_assert(std::is_class<Event>::value,"Provided event-type template parameter is not a class or a struct!");
            static_assert(std::is_class<Track>::value,"Provided track-type template parameter is not a class or a struct!");
            static_assert(std::is_class<Pair>::value,"Provided pair-type template parameter is not a class or a struct!");

            private:
                std::size_t fBufferSize;
                bool fWaitForBuffer,fEventHashingFunctionIsDefined,fPairHashingFunctionIsDefined,fPairCutFunctionIsDefined;
                std::size_t fDimensions;
                std::random_device fRandomDevice;
                std::mt19937 fGenerator;
                //map of similar events (of maximal size give by fBufferSize); each event has a single track (randomly chosen from the collection)
                std::map<std::size_t, std::deque<std::pair<Event, Track> > > fSimilarityMap;
                std::function<std::size_t(const Event &)> fEventHashingFunction;
                std::function<std::size_t(const Pair &)> fPairHashingFunction;
                std::function<bool(const Pair &)> fPairCutFunction;
                /**
                 * @brief Create pairs of identical particles from given tracks
                 * 
                 * @param tracks tracks vector
                 * @return std::vector<Pair> vector of pairs
                 */
                [[nodiscard]] std::vector<Pair> MakePairs(const std::vector<Track> &tracks);
                /**
                 * @brief Divide pairs into corresponding category (given by the pair hash)
                 * 
                 * @param pairs pairs vector
                 * @return std::map<std::size_t, std::vector<Pair> > map of sorted vectors (each "branch"/bucket is a single group of similar pairs)
                 */
                [[nodiscard]] std::map<std::size_t, std::vector<Pair> > SortPairs(const std::vector<Pair> &pairs);

            public:
                /**
                 * @brief Default constructor. Create mixer object with: buffer size 10, no wait for full buffer, no hashing in event and pair, and no pair cut
                 * 
                 */
                constexpr JJFemtoMixer() : fBufferSize(10), 
                                fWaitForBuffer(false), 
                                fEventHashingFunctionIsDefined(false),
                                fPairHashingFunctionIsDefined(false),
                                fPairCutFunctionIsDefined(false),
                                fGenerator(fRandomDevice()),
                                fEventHashingFunction([](const Event &){return 0;}),
                                fPairHashingFunction([](const Pair &){return 0;}),
                                fPairCutFunction([](const Pair &){return false;}) {}

                /**
                 * @brief Set the Event Hashing Function object.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetEventHashingFunction(const std::function<std::size_t(const Event &)> &func) noexcept {fEventHashingFunction = func; fEventHashingFunctionIsDefined = true;}
                /**
                 * @brief Set the Event Hashing Function object.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetEventHashingFunction(std::function<std::size_t(const Event &)> &&func) noexcept {fEventHashingFunction = std::move(func);  fEventHashingFunctionIsDefined = true;}
                /**
                 * @brief Get the corresponding hash for given Event object.
                 * 
                 * @param obj Event-type object
                 * @return std::size_t 
                 */
                [[nodiscard]] constexpr std::size_t GetEventHash(const Event &obj) const noexcept {return fEventHashingFunction(obj);}
                /**
                 * @brief Set the Pair Hashing Function object.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetPairHashingFunction(const std::function<std::size_t(const Pair &)> &func) noexcept {fPairHashingFunction = func; fPairHashingFunctionIsDefined = true;}
                /**
                 * @brief Set the Pair Hashing Function object.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetPairHashingFunction(std::function<std::size_t(const Pair &)> &&func) noexcept {fPairHashingFunction = std::move(func); fPairHashingFunctionIsDefined = true;}
                /**
                 * @brief Get the corresponding hash for given Pair object.
                 * 
                 * @param obj Pair-type object.
                 * @return std::size_t 
                 */
                [[nodiscard]] constexpr std::size_t GetPairHash(const Pair &obj) const noexcept {return fPairHashingFunction(obj);}
                /**
                 * @brief Set the Pair Cutting Function object. The function should return true if pair should be rejected and false if accepted.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetPairCuttingFunction(const std::function<bool(const Pair &)> &func) noexcept {fPairCutFunction = func; fPairCutFunctionIsDefined = true;}
                /**
                 * @brief Set the Pair Cutting Function object. The function should return true if pair should be rejected and false if accepted.
                 * 
                 * @param func Function object, can be lambda, standard function or std::function object.
                 */
                constexpr void SetPairCuttingFunction(std::function<bool(const Pair &)> &&func) noexcept {fPairCutFunction = std::move(func); fPairCutFunctionIsDefined = true;}
                /**
                 * @brief Get the result of pair selection function for a given Pair object.
                 * 
                 * @param obj Pair-type object.
                 * @return true Pair is rejected.
                 * @return false Pair is not rejected.
                 */
                [[nodiscard]] constexpr bool GetPairCutResult(const Pair &obj) const noexcept {return fPairCutFunction(obj);}
                /**
                 * @brief Set the max mixing buffer size for each "branch".
                 * 
                 * @param buffer Max buffer size.
                 */
                constexpr void SetMaxBufferSize(const std::size_t &buffer) noexcept {fBufferSize = buffer;}
                /**
                 * @brief Set the max mixing buffer size for each "branch".
                 * 
                 * @param buffer Max buffer size.
                 */
                constexpr void SetMaxBufferSize(std::size_t &&buffer) noexcept {fBufferSize = std::move(buffer);}
                /**
                 * @brief Get the max mixing buffer size.
                 * 
                 * @return std::size_t Max buffer size.
                 */
                [[nodiscard]] constexpr std::size_t GetMaxBufferSize() const noexcept {return fBufferSize;}
                /**
                 * @brief Define whether mixing should occur only for "branches" with size equal (if true) / equal or less than (if flase) the max buffer size.
                 * 
                 * @param isFixed Buffer size usage flag.
                 */
                constexpr void SetFixedBuffer(const bool &isFixed) noexcept {fWaitForBuffer = isFixed;}
                /**
                 * @brief Define whether mixing should occur only for "branches" with size equal (if true) / equal or less than (if flase) the max buffer size.
                 * 
                 * @param isFixed Buffer size usage flag.
                 */
                constexpr void SetFixedBuffer(bool &&isFixed) noexcept {fWaitForBuffer = std::move(isFixed);}
                /**
                 * @brief Get the mixing buffer flag.
                 * 
                 * @return true Mixing with fixed size.
                 * @return false Mixnig with not fixed size.
                 */
                [[nodiscard]] constexpr bool GetBufferState() const noexcept {return fWaitForBuffer;};
                /**
                 * @brief Prints to the standard output information about current setup of JJFemtoMixer.
                 * 
                 */
                void PrintSettings() const;
                /**
                 * @brief Prints to the standard output information about the amounts of tracks and events currently stored in JJFemtoMixer.
                 * 
                 */
                void PrintStatus() const;
                /**
                 * @brief Add currently processed event to the mixer.
                 * 
                 * @param event Current event.
                 * @param tracks Tracks from the current event.
                 * @return std::map<std::size_t, std::vector<Pair> > Sorted pairs from provided tracks for given event.
                 */
                [[nodiscard]] std::map<std::size_t, std::vector<Pair> > AddEvent(const Event &event, const std::vector<Track> &tracks);
                /**
                 * @brief Get the sorted pairs which come from similar events, but not from this event.
                 * 
                 * @param event Current event (the event from which we don't want to get tracks).
                 * @return std::map<std::size_t, std::vector<Pair> > Sorted pairs from stored tracks for similar events.
                 */
                [[nodiscard]] std::map<std::size_t, std::vector<Pair> > GetSimilarPairs(const Event &event);
        };

        template<typename Event, typename Track, typename Pair>
        std::vector<Pair> JJFemtoMixer<Event,Track,Pair>::MakePairs(const std::vector<Track> &tracks)
        {
            bool reverse = false;
            const std::size_t trackVecSize = tracks.size();

            std::vector<Pair> tmpVector;
            tmpVector.reserve((1 + trackVecSize) * trackVecSize / 2); // reserve the expected vector size (arithmetic sum, from 1 to N)

            for (std::size_t iter1 = 0; iter1 < trackVecSize; ++iter1)
                for (std::size_t iter2 = iter1 + 1; iter2 < trackVecSize; ++iter2)
                {
                    if (reverse)
                        tmpVector.push_back(Pair(tracks.at(iter2),tracks.at(iter1)));
                    else
                        tmpVector.push_back(Pair(tracks.at(iter1),tracks.at(iter2)));

                    reverse = !reverse; // reverse the order of tracks every other time (get rid of the bias from the track sorter)
                }

            return tmpVector;
        }

        template<typename Event, typename Track, typename Pair>
        std::map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::SortPairs(const std::vector<Pair> &pairs)
        {
            std::size_t currentPairHash = 0;
            std::map<std::size_t, std::vector<Pair> > pairMap;

            for (const auto &pair : pairs)
            {
                if (fPairCutFunction(pair))
                    continue;

                currentPairHash = fPairHashingFunction(pair);

                // an entry for given currentPairHash may already exist so we must check if that's the case
                if (pairMap.find(currentPairHash) != pairMap.end())
                {
                    pairMap.at(currentPairHash).push_back(pair);
                }
                else
                {
                    pairMap.emplace(currentPairHash,std::vector<Pair>(1,pair));
                }
            }

            return pairMap;
        }

        template<typename Event, typename Track, typename Pair>
        void JJFemtoMixer<Event,Track,Pair>::PrintSettings() const
        {
            std::cout << "\n------=========== JJFemtoMixer Settings ===========------\n";
            std::cout << "Max Background Mixing Buffer Size: " << fBufferSize << ((fWaitForBuffer) ? " (FIXED)\n" : " (LOOSE)\n");
            std::cout << "Event Hashing Function: " << ((fEventHashingFunctionIsDefined) ? " User-defined\n" : " Not set\n");
            std::cout << "Pair Hashing Function: " << ((fPairHashingFunctionIsDefined) ? " User-defined\n" : " Not set\n");
            std::cout << "Pair Rejection Function: " << ((fPairCutFunctionIsDefined) ? " User-defined\n" : " Not set\n");
            std::cout << "------==============================================------\n" << std::endl;
        }

        template<typename Event, typename Track, typename Pair>
        void JJFemtoMixer<Event,Track,Pair>::PrintStatus() const
        {
            std::cout << "\n------=========== JJFemtoMixer Status ===========------\n";
            std::cout << "Currently stored tracks:\n";
            for (const auto &[key,val] : fSimilarityMap)
                std::cout << val.size() << "/" << fBufferSize << "\t event hash: " << key << "\n";
            std::cout << "------=============================================------\n" << std::endl;
        }

        template<typename Event, typename Track, typename Pair>
        std::map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::AddEvent(const Event &event, const std::vector<Track> &tracks)
        {
            std::uniform_int_distribution<int> dist(0,tracks.size()-1);
            std::size_t evtHash = fEventHashingFunction(event);
            std::pair<Event, Track> trackPair{event,tracks.at(dist(fGenerator))};

            // an entry for given evtHash may not exist, so we must check if that's the case
            if (fSimilarityMap.find(evtHash) == fSimilarityMap.end())
            {
                fSimilarityMap.emplace(evtHash,std::deque<std::pair<Event, Track> >(1,trackPair));
            }
            else
            {
                fSimilarityMap.at(evtHash).push_back(trackPair);
                if (fSimilarityMap.at(evtHash).size() > fBufferSize)
                    fSimilarityMap.at(evtHash).pop_front(); 
            }

            return SortPairs(MakePairs(tracks));
        }

        template<typename Event, typename Track, typename Pair>
        std::map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::GetSimilarPairs(const Event &event)
        {
            std::vector<Track> outputVector;
            std::size_t evtHash = fEventHashingFunction(event);

            if (fSimilarityMap.at(evtHash).size() == fBufferSize || fWaitForBuffer == false)
            {
                for (const auto &[evt,trck] : fSimilarityMap.at(evtHash))
                    if (evt != event)
                        outputVector.push_back(trck);
            }

            return SortPairs(MakePairs(outputVector));
        }
    } // namespace Mixing
    
#endif