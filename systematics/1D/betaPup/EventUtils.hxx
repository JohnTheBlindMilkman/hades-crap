#ifndef EventUtils_hxx
    #define EventUtils_hxx

    #include <string>
    #include <memory>

    #include "EventCandidate.hxx"
    #include "JJUtils.hxx"

    namespace Mixing
    {
        class EventGrouping
        {
            private:

            public:
                [[nodiscard]] std::string GetEventIndex(const std::shared_ptr<Selection::EventCandidate> &evt) const noexcept
                {
                    return JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetNCharged()/10),2) + 
                		JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
                }
                [[nodiscard]] std::function<std::string (const std::shared_ptr<Selection::EventCandidate> &)> MakeEventGroupingFunction() const noexcept
                {
                    return [this](const std::shared_ptr<Selection::EventCandidate> &evt){return this->GetEventIndex(evt);};
                }
        };
        
    } // namespace Mixing
    

#endif
