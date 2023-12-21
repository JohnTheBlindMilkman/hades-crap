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
                bool fWaitForBuffer;
                std::vector<std::size_t> fPIDList;
                std::size_t fDimensions;
                //map of similar events (of maximal size give by fBufferSize); each event has a map of similar tracks
                std::unordered_map<std::size_t, std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > > > fSimilarityMap;
                std::function<const std::size_t(const T &)> fEventHashingFunction;
                std::function<const std::size_t(const U &)> fTrackHashingFunction;

                std::unordered_map<std::size_t, std::vector<U> > SortTracks(const std::vector<U> &tracks);

            public:
                JJMultiEventMixer() : fBufferSize(0), 
                                    fWaitForBuffer(false), 
                                    fPIDList({}), 
                                    fDimensions(fPIDList.size()), 
                                    fSimilarityMap({}),
                                    fEventHashingFunction([](const T &){return 0;}), 
                                    fTrackHashingFunction([](const U &){return 0;}) {}

                void SetEventHashingFunction(std::function<const std::size_t(const T &)> func) {fEventHashingFunction = func;}
                const std::size_t GetEventHash(const T &obj) const {return fEventHashingFunction(obj);}

                void SetTrackHashingFunction(std::function<const std::size_t(const U &)> func) {fTrackHashingFunction = func;}
                const std::size_t GetTrackHash(const U &obj) const {return fTrackHashingFunction(obj);}

                // To-Do: implement PIDs
                void SetPIDs(const std::vector<std::size_t> &particleList) {fPIDList = particleList;}
                std::vector<std::size_t> GetPIDList() const {return fPIDList;}

                void SetMaxBufferSize(std::size_t buffer) {fBufferSize = buffer;}
                std::size_t GetMaxBufferSize() const {return fBufferSize;}

                void SetFixedBuffer(bool isFixed) {fWaitForBuffer = isFixed;}
                bool GetBufferState() const {return fWaitForBuffer};

                void AddEvent(const T &event, const std::vector<U> &tracks);
                std::vector<U> GetSimilarTracks(const T &event, const U &track) const;
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
                    std::vector<U> tmpVec;
                    tmpVec.push_back(track);
                    fSimilarityMap.emplace(currentTrackHash,tmpVec);
                }
            }

            return trackMap;
        }

        template<typename T, typename U>
        void JJMultiEventMixer<T,U>::AddEvent(const T &event, const std::vector<U> &tracks)
        {
            std::size_t evtHash = fEventHashingFunction(event);
            std::unordered_map<T, std::unordered_map<std::size_t,std::vector<U> > > trackMap;

            trackMap.emplace(event,SortTracks(tracks));

            // an entry for given evtHash may already exist so we must check if that's the case
            if (fSimilarityMap.find(evtHash) != fSimilarityMap.end())
            {
                fSimilarityMap.at(evtHash).push_back(trackMap);
            }
            else
            {
                std::deque<std::unordered_map<T, std::unordered_map<std::size_t, std::vector<U> > > > tmpDeque;
                tmpDeque.push_back(trackMap);
                fSimilarityMap.emplace(evtHash,tmpDeque);
            }

            if (fSimilarityMap[evtHash].size() > fBufferSize)
                fSimilarityMap[evtHash].pop_front(); 
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
                    if (fSimilarityMap.at(evtHash).at(evtIter).front() != fSimilarityMap.at(evtHash).at(evtIter).at(event))
                        if (fSimilarityMap.at(evtHash).at(evtIter).front().find(trackHash) != fSimilarityMap.at(evtHash).at(evtIter).front().end())
                            outputVector.push_back(fSimilarityMap.at(evtHash).at(evtIter).front().at(trackHash).front());
            }

            return outputVector;
        }

    } // namespace Mixing
    

#endif