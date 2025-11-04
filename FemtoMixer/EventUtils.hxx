/**
 * @file EventUtils.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Collection of classes for consistent use of event grouping
 * @version 0.1.0
 * @date 2025-10-30
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef EventUtils_hxx
    #define EventUtils_hxx

    #include <string>
    #include <memory>

    #include "EventCandidate.hxx"
    #include "JJUtils.hxx"

    namespace Mixing
    {
        /**
         * @brief Class storing functions used for event grouping
         * 
         */
        class EventGrouping
        {
            private:

            public:
                /**
                 * @brief Calculates and returns group ID to which the EventCandidate object is assigned to
                 * 
                 * @param evt pointer to the EventCandidate
                 * @return group ID
                 */
                [[nodiscard]] std::string GetEventIndex(const std::shared_ptr<Selection::EventCandidate> &evt) const noexcept
                {
                    return JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetNCharged()/10),2) + 
                		JJUtils::to_fixed_size_string(static_cast<std::size_t>(evt->GetPlate()),2);
                }
                /**
                 * @brief Creates a wrapper for GetEventIndex
                 * 
                 * @return std::function
                 */
                [[nodiscard]] std::function<std::string (const std::shared_ptr<Selection::EventCandidate> &)> MakeEventGroupingFunction() const noexcept
                {
                    return [this](const std::shared_ptr<Selection::EventCandidate> &evt){return this->GetEventIndex(evt);};
                }
        };
        
    } // namespace Mixing
    

#endif