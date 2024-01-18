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

    namespace Mixing
    {
        template<typename Event, typename Track, typename Pair>
        class JJFemtoMixer
        {
            private:
                std::size_t fBufferSize;
                bool fWaitForBuffer;
                std::vector<long> fPIDList;
                std::size_t fDimensions;
                //map of similar events (of maximal size give by fBufferSize); each event has a map of similar tracks
                std::unordered_map<std::size_t, std::deque<std::unordered_map<Event, std::unordered_map<std::size_t, std::vector<Pair> > > > > fSimilarityMap;
                std::function<std::size_t(const Event &)> fEventHashingFunction;
                std::function<std::size_t(const Pair &)> fPairHashingFunction;

                std::vector<Pair> MakePairs(const std::vector<Track> &tracks);
                std::unordered_map<std::size_t, std::vector<Pair> > SortPairs(const std::vector<Pair> &tracks);

            public:
                JJFemtoMixer() : fBufferSize(10), 
                                    fWaitForBuffer(false), 
                                    fPIDList({}), 
                                    fDimensions(0), 
                                    fEventHashingFunction([](const Event &){return 0;}), 
                                    fPairHashingFunction([](const Pair &){return 0;}) {}

                void SetEventHashingFunction(std::function<std::size_t(const Event &)> func) {fEventHashingFunction = func;}
                std::size_t GetEventHash(const Event &obj) const {return fEventHashingFunction(obj);}

                void SetPairHashingFunction(std::function<std::size_t(const Pair &)> func) {fPairHashingFunction = func;}
                std::size_t GetPairHash(const Pair &obj) const {return fPairHashingFunction(obj);}

                // To-Do: implement PIDs
                void SetPIDs(const std::vector<long> &particleList) {fPIDList = particleList;}
                std::vector<long> GetPIDList() const {return fPIDList;}

                void SetMaxBufferSize(std::size_t buffer) {fBufferSize = buffer;}
                std::size_t GetMaxBufferSize() const {return fBufferSize;}

                void SetFixedBuffer(bool isFixed) {fWaitForBuffer = isFixed;}
                bool GetBufferState() const {return fWaitForBuffer;};

                std::unordered_map<std::size_t, std::vector<Pair> > AddEvent(const Event &event, const std::vector<Track> &tracks);
                std::vector<Pair> GetSimilarPairs(const Event &event, const Pair &pair) const;
        };

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
            // good idea would be to simplify the complex data structures:
            // using key = std::size_t;
            // using  value = std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > >;
            std::size_t evtHash = fEventHashingFunction(event);
            std::unordered_map<std::size_t, std::vector<Pair> > tmpPairMap = SortPairs(MakePairs(tracks));
            std::unordered_map<Event, std::unordered_map<std::size_t, std::vector<Pair> > > pairMap;

            pairMap.emplace(event,tmpPairMap);

            // an entry for given evtHash may not exist, so we must check if that's the case
            if (fSimilarityMap.find(evtHash) == fSimilarityMap.end())
            {
                std::deque<std::unordered_map<Event, std::unordered_map<std::size_t, std::vector<Pair> > > > tmpQueue;
                tmpQueue.push_back(pairMap);
                fSimilarityMap.emplace(evtHash,tmpQueue);
            }
            else
            {
                fSimilarityMap.at(evtHash).push_back(pairMap);
                if (fSimilarityMap.at(evtHash).size() > fBufferSize)
                    fSimilarityMap.at(evtHash).pop_front(); 
            }

            return tmpPairMap;
        }

        template<typename Event, typename Track, typename Pair>
        std::vector<Pair> JJFemtoMixer<Event,Track,Pair>::GetSimilarPairs(const Event &event, const Pair &pair) const
        {
            std::vector<Pair> outputVector = {};
            std::size_t evtHash = fEventHashingFunction(event);
            std::size_t trackHash = fPairHashingFunction(pair);

            if (fSimilarityMap.at(evtHash).size() == fBufferSize || fWaitForBuffer == false)
            {
                for (std::size_t evtIter = 0; evtIter < fSimilarityMap.at(evtHash).size(); ++evtIter)
                    if (fSimilarityMap.at(evtHash).at(evtIter).begin()->first != event)
                        if (fSimilarityMap.at(evtHash).at(evtIter).begin()->second.find(trackHash) != fSimilarityMap.at(evtHash).at(evtIter).begin()->second.end())
                            outputVector.push_back(fSimilarityMap.at(evtHash).at(evtIter).begin()->second.at(trackHash).front());
            }

            return outputVector;
        }

    } // namespace Mixing
    

#endif