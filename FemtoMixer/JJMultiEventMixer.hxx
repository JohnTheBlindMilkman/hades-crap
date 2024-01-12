#ifndef JJMultiEventMixer_hxx
    #define JJMultiEventMixer_hxx

    #include <vector>
    #include <deque>
    #include <unordered_map>
    #include <functional>

    namespace Mixing
    {
        template<typename T, typename U>
        class JJMultiEventMixer
        {
            private:
                std::size_t fBufferSize;
                bool fWaitForBuffer, fIsEventHashingSet, fIsTrackHashingSet;
                std::vector<long> fPIDList;
                std::size_t fDimensions;
                //map of similar events (of maximal size give by fBufferSize); each event has a map of similar tracks
                std::unordered_map<std::size_t, std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > > > fSimilarityMap;
                std::function<std::size_t(const T &)> fEventHashingFunction;
                std::function<std::size_t(const U &)> fTrackHashingFunction;

                std::unordered_map<std::size_t, std::vector<U> > SortTracks(const std::vector<U> &tracks);

            public:
                JJMultiEventMixer() : fBufferSize(10), 
                                    fWaitForBuffer(false), 
                                    fIsEventHashingSet(false),
                                    fIsTrackHashingSet(false),
                                    fPIDList({}), 
                                    fDimensions(0), 
                                    /* fSimilarityMap(), */
                                    fEventHashingFunction([](const T &){return 0;}), 
                                    fTrackHashingFunction([](const U &){return 0;}) {}

                void SetEventHashingFunction(std::function<std::size_t(const T &)> func) {fIsEventHashingSet = true; fEventHashingFunction = func;}
                std::size_t GetEventHash(const T &obj) const {return fEventHashingFunction(obj);}

                void SetTrackHashingFunction(std::function<std::size_t(const U &)> func) {fIsTrackHashingSet = true; fTrackHashingFunction = func;}
                std::size_t GetTrackHash(const U &obj) const {return fTrackHashingFunction(obj);}

                // To-Do: implement PIDs
                void SetPIDs(const std::vector<long> &particleList) {fPIDList = particleList;}
                std::vector<long> GetPIDList() const {return fPIDList;}

                void SetMaxBufferSize(std::size_t buffer) {fBufferSize = buffer;}
                std::size_t GetMaxBufferSize() const {return fBufferSize;}

                void SetFixedBuffer(bool isFixed) {fWaitForBuffer = isFixed;}
                bool GetBufferState() const {return fWaitForBuffer;};

                std::unordered_map<std::size_t, std::vector<U> > AddEvent(const T &event, const std::vector<U> &tracks);
                std::vector<U> GetSimilarTracks(const T &event, const U &track) const;
                std::unordered_map<std::size_t, std::vector<U> > GetSortedTracks(const T &event);
        };

        template<typename T, typename U>
        std::unordered_map<std::size_t, std::vector<U> > JJMultiEventMixer<T,U>::SortTracks(const std::vector<U> &tracks)
        {
            std::size_t currentTrackHash = 0;
            std::unordered_map<std::size_t, std::vector<U> > trackMap;

            for (const auto &track : tracks)
            {
                currentTrackHash = fTrackHashingFunction(track);

                // an entry for given currentTrackHash may already exist so we must check if that's the case
                if (trackMap.find(currentTrackHash) != trackMap.end())
                {
                    trackMap.at(currentTrackHash).push_back(track);
                }
                else
                {
                    std::vector<U> tmpVec{track};
                    //tmpVec.push_back(track);
                    trackMap.emplace(currentTrackHash,tmpVec);
                }
            }

            return trackMap;
        }

        template<typename T, typename U>
        std::unordered_map<std::size_t, std::vector<U> > JJMultiEventMixer<T,U>::AddEvent(const T &event, const std::vector<U> &tracks)
        {
            if (!fIsEventHashingSet)
            {
                throw std::runtime_error("Event hashnig is not set.");
            }
            if (!fIsTrackHashingSet)
            {
                throw std::runtime_error("Track hashing is not set.");
            }

            // good idea would be to simplify the complex data structures:
            // using key = std::size_t;
            // using  value = std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > >;
            std::size_t evtHash = fEventHashingFunction(event);
            std::unordered_map<std::size_t, std::vector<U> > tmpTrackMap = SortTracks(tracks);
            std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > trackMap;

            trackMap.emplace(event,tmpTrackMap);

            // an entry for given evtHash may not exist, so we must check if that's the case
            if (fSimilarityMap.find(evtHash) == fSimilarityMap.end())
            {
                std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > > tmpQueue;
                tmpQueue.push_back(trackMap);
                fSimilarityMap.emplace(evtHash,tmpQueue);
            }
            else
            {
                fSimilarityMap.at(evtHash).push_back(trackMap);
                if (fSimilarityMap.at(evtHash).size() > fBufferSize)
                    fSimilarityMap.at(evtHash).pop_front(); 
            }

            return tmpTrackMap;
        }

        template<typename T, typename U>
        std::vector<U> JJMultiEventMixer<T,U>::GetSimilarTracks(const T &event, const U &track) const
        {
            std::vector<U> outputVector = {};
            std::size_t evtHash = fEventHashingFunction(event);
            std::size_t trackHash = fTrackHashingFunction(track);

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