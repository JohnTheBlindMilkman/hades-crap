/**
 * @file JJUtils.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Class containing all bits and bobs used in my code
 * @version 0.1.0
 * @date 2024-11-29
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef JJUtils_hxx
    #define JJUtils_hxx

    #include  <random>
    #include  <iterator>

    namespace JJUtils
    {
        namespace Detail
        {
            template<typename Iter, typename RandomGenerator>
            Iter select_randomly(Iter start, Iter end, RandomGenerator& g) 
            {
                std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
                std::advance(start, dis(g));
                return start;
            }
        }

        /**
         * @brief Select random element from an STL container
         * 
         * @tparam Iter iterator type
         * @param start begin iterator
         * @param end end iterator
         * @return iterator pointing to the n-th element of the container
         */
        template<typename Iter>
        Iter select_randomly(Iter start, Iter end) 
        {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            return Detail::select_randomly(start, end, gen);
        }
    }

#endif