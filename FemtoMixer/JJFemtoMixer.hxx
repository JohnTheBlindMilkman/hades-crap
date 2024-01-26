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

                std::vector<Pair> MakePairs(const std::vector<Track> &tracks, const float &EPangle);
                std::vector<Pair> MakePairs(const std::vector<Track> &tracks);
                std::unordered_map<std::size_t, std::vector<Pair> > SortPairs(const std::vector<Pair> &pairs);

            public:
                JJFemtoMixer() : fBufferSize(10), 
                                fWaitForBuffer(false), 
                                fDimensions(0), 
                                fGenerator(fRandomDevice()),
                                fEventHashingFunction([](const Event &){return 0;}),
                                fPairHashingFunction([](const Pair &){return 0;}) {}

                void SetEventHashingFunction(std::function<std::size_t(const Event &)> func) {fEventHashingFunction = func;}
                std::size_t GetEventHash(const Event &obj) const {return fEventHashingFunction(obj);}

                void SetPairHashingFunction(std::function<std::size_t(const Pair &)> func) {fPairHashingFunction = func;}
                std::size_t GetPairHash(const Pair &obj) const {return fPairHashingFunction(obj);}

                void SetMaxBufferSize(std::size_t buffer) {fBufferSize = buffer;}
                std::size_t GetMaxBufferSize() const {return fBufferSize;}

                void SetFixedBuffer(bool isFixed) {fWaitForBuffer = isFixed;}
                bool GetBufferState() const {return fWaitForBuffer;};

                std::unordered_map<std::size_t, std::vector<Pair> > AddEvent(const Event &event, const std::vector<Track> &tracks);
                std::unordered_map<std::size_t, std::vector<Pair> > GetSimilarPairs(const Event &event);
                std::unordered_map<std::size_t, std::vector<Pair> > GetSimilarPairs(const Event &event, const float &EPangle);
        };

        template<typename Event, typename Track, typename Pair>
        std::vector<Pair> JJFemtoMixer<Event,Track,Pair>::MakePairs(const std::vector<Track> &tracks, const float &EPangle)
        {
            bool reverse = true;
            const std::size_t trackVecSize = tracks.size();

            std::vector<Pair> tmpVector;
            tmpVector.reserve((1 + trackVecSize) * trackVecSize / 2); // reserve the expected vector size (arithmetic sum, from 1 to N)

            for (std::size_t iter1 = 0; iter1 < trackVecSize; ++iter1)
                for (std::size_t iter2 = iter1 + 1; iter2 < trackVecSize; ++iter2)
                {
                    if (reverse)
                        tmpVector.push_back(Pair(tracks.at(iter2),tracks.at(iter1),EPangle));
                    else
                        tmpVector.push_back(Pair(tracks.at(iter1),tracks.at(iter2),EPangle));

                    reverse = !reverse; // reverse the order of tracks every other time (get rid of the bias from the track sorter)
                }

            return tmpVector;
        }

        template<typename Event, typename Track, typename Pair>
        std::vector<Pair> JJFemtoMixer<Event,Track,Pair>::MakePairs(const std::vector<Track> &tracks)
        {
            bool reverse = true;
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
                currentPairHash = fPairHashingFunction(pair);

                // an entry for given currentPairHash may already exist so we must check if that's the case
                if (pairMap.find(currentPairHash) != pairMap.end())
                {
                    pairMap.at(currentPairHash).push_back(pair);
                }
                else
                {
                    std::vector<Pair> tmpVec{pair};
                    pairMap.emplace(currentPairHash,tmpVec);
                }
            }

            return pairMap;
        }

        template<typename Event, typename Track, typename Pair>
        std::unordered_map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::AddEvent(const Event &event, const std::vector<Track> &tracks)
        {
            std::uniform_int_distribution<int> dist(0,tracks.size()-1);
            std::size_t evtHash = fEventHashingFunction(event);
            std::unordered_map<std::size_t, std::vector<Pair> > tmpPairMap = SortPairs(MakePairs(tracks));
            std::pair<Event, Track> trackPair{event,tracks.at(dist(fGenerator))};

            // an entry for given evtHash may not exist, so we must check if that's the case
            if (fSimilarityMap.find(evtHash) == fSimilarityMap.end())
            {
                std::deque<std::pair<Event, Track> > tmpQueue;
                tmpQueue.push_back(trackPair);
                fSimilarityMap.emplace(evtHash,tmpQueue);
            }
            else
            {
                fSimilarityMap.at(evtHash).push_back(trackPair);
                if (fSimilarityMap.at(evtHash).size() > fBufferSize)
                    fSimilarityMap.at(evtHash).pop_front(); 
            }

            return tmpPairMap;
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

        template<typename Event, typename Track, typename Pair>
        std::unordered_map<std::size_t, std::vector<Pair> > JJFemtoMixer<Event,Track,Pair>::GetSimilarPairs(const Event &event, const float &EPangle)
        {
            std::vector<Track> outputVector;
            std::size_t evtHash = fEventHashingFunction(event);

            if (fSimilarityMap.at(evtHash).size() == fBufferSize || fWaitForBuffer == false)
            {
                for (std::size_t evtIter = 0; evtIter < fSimilarityMap.at(evtHash).size(); ++evtIter)
                    if (fSimilarityMap.at(evtHash).at(evtIter).first != event)
                        outputVector.push_back(fSimilarityMap.at(evtHash).at(evtIter).second);
            }

            return SortPairs(MakePairs(outputVector,EPangle));
        }
    } // namespace Mixing
    

#endif