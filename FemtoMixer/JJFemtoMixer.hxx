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
    #include <unordered_map>
    #include <functional>
    #include <random>

    namespace Mixing
    {
        template<typename Event, typename Track, typename Pair>
        class JJFemtoMixer
        {
            private:
                std::size_t fBufferSize;
                bool fWaitForBuffer;
                std::size_t fDimensions;
                std::random_device fRandomDevice;
                std::mt19937 fGenerator;
                //map of similar events (of maximal size give by fBufferSize); each event has a single track (randomly chosen from the collection)
                std::unordered_map<std::size_t, std::deque<std::pair<Event, Track> > > fSimilarityMap;
                std::function<std::size_t(const Event &)> fEventHashingFunction;
                std::function<std::size_t(const Pair &)> fPairHashingFunction;
                std::function<bool(const Pair &)> fPairCutFunction;

                std::vector<Pair> MakePairs(const std::vector<Track> &tracks);
                std::unordered_map<std::size_t, std::vector<Pair> > SortPairs(const std::vector<Pair> &pairs);

            public:
                constexpr JJFemtoMixer() : fBufferSize(10), 
                                fWaitForBuffer(false), 
                                fGenerator(fRandomDevice()),
                                fEventHashingFunction([](const Event &){return 0;}),
                                fPairHashingFunction([](const Pair &){return 0;}),
                                fPairCutFunction([](const Pair &){return false;}) {}

                constexpr void SetEventHashingFunction(const std::function<std::size_t(const Event &)> &func) {fEventHashingFunction = func;}
                constexpr void SetEventHashingFunction(std::function<std::size_t(const Event &)> &&func) {fEventHashingFunction = std::move(func);}
                constexpr std::size_t GetEventHash(const Event &obj) const {return fEventHashingFunction(obj);}

                constexpr void SetPairHashingFunction(const std::function<std::size_t(const Pair &)> &func) {fPairHashingFunction = func;}
                constexpr void SetPairHashingFunction(std::function<std::size_t(const Pair &)> &&func) {fPairHashingFunction = std::move(func);}
                constexpr std::size_t GetPairHash(const Pair &obj) const {return fPairHashingFunction(obj);}

                constexpr void SetPairCuttingFunction(const std::function<bool(const Pair &)> &func) {fPairCutFunction = func;}
                constexpr void SetPairCuttingFunction(std::function<bool(const Pair &)> &&func) {fPairCutFunction = std::move(func);}
                constexpr bool GetPairCutResult(const Pair &obj) const {return fPairCutFunction(obj);}

                constexpr void SetMaxBufferSize(const std::size_t &buffer) {fBufferSize = buffer;}
                constexpr void SetMaxBufferSize(std::size_t &&buffer) {fBufferSize = std::move(buffer);}
                constexpr std::size_t GetMaxBufferSize() const {return fBufferSize;}

                constexpr void SetFixedBuffer(const bool &isFixed) {fWaitForBuffer = isFixed;}
                constexpr void SetFixedBuffer(bool &&isFixed) {fWaitForBuffer = std::move(isFixed);}
                constexpr bool GetBufferState() const {return fWaitForBuffer;};

                std::unordered_map<std::size_t, std::vector<Pair> > AddEvent(const Event &event, const std::vector<Track> &tracks);
                std::unordered_map<std::size_t, std::vector<Pair> > GetSimilarPairs(const Event &event);
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
        std::unordered_map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::SortPairs(const std::vector<Pair> &pairs)
        {
            std::size_t currentPairHash = 0;
            std::unordered_map<std::size_t, std::vector<Pair> > pairMap;

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
                    //std::vector<Pair> tmpVec{pair};
                    pairMap.emplace(currentPairHash,std::vector<Pair>(1,pair));
                }
            }

            return pairMap;
        }

        template<typename Event, typename Track, typename Pair>
        std::unordered_map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::AddEvent(const Event &event, const std::vector<Track> &tracks)
        {
            std::uniform_int_distribution<int> dist(0,tracks.size()-1);
            std::size_t evtHash = fEventHashingFunction(event);
            std::pair<Event, Track> trackPair{event,tracks.at(dist(fGenerator))};

            // an entry for given evtHash may not exist, so we must check if that's the case
            if (fSimilarityMap.find(evtHash) == fSimilarityMap.end())
            {
                //std::deque<std::pair<Event, Track> > tmpQueue;
                //tmpQueue.push_back(trackPair);
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
        std::unordered_map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::GetSimilarPairs(const Event &event)
        {
            std::vector<Track> outputVector;
            std::size_t evtHash = fEventHashingFunction(event);

            if (fSimilarityMap.at(evtHash).size() == fBufferSize || fWaitForBuffer == false)
            {
                for (std::size_t evtIter = 0; evtIter < fSimilarityMap.at(evtHash).size(); ++evtIter)
                    if (fSimilarityMap.at(evtHash).at(evtIter).first != event)
                        outputVector.push_back(fSimilarityMap.at(evtHash).at(evtIter).second);
            }

            return SortPairs(MakePairs(outputVector));
        }
    } // namespace Mixing
    

#endif