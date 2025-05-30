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

    #include <random>
    #include <iterator>
    #include <string>
    #include <sstream>

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
         * @return iterator pointing to a randomly selected element of the container
         */
        template<typename Iter>
        Iter select_randomly(Iter start, Iter end) 
        {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            return Detail::select_randomly(start, end, gen);
        }

        /**
         * @brief Turn value "val" into string with leading zeros. E.g. int val = 1 and width = 2, then the function returns "01"
         * 
         * @tparam T 
         * @param val any type representing a number (strings should work too?)
         * @param width length of the string (val and the leading zeros)
         * @return std::string of size equal to "width", containing the "val" with leading zeros
         */
        template <typename T>
        std::string to_fixed_size_string(T val, std::size_t width)
        {
            std::ostringstream oss;
            oss.width(width);
            oss.fill('0');
            oss << val;
            return oss.str();
        }
        
        /**
         * @brief Create new container of size std::min(size_of_container1,size_of_container2) which holds items from both containers at each index
         * 
         * @tparam T type of the containers
         * @param cont1 first container
         * @param cont2 second container
         * @return std::vector<std::pair<T,T> > output container
         */
        template <template<class> typename Container, typename T>
        Container<std::pair<T,T> > zip_to_smaller(const Container<T> &cont1, const Container<T> &cont2)
        {
            std::size_t maxIter = std::min(cont1.size(),cont2.size());
            Container<std::pair<T,T> > contOut(std::pair<T,T>(),maxIter);

            for (std::size_t i = 0; i < maxIter; ++i)
            {
                contOut[i] = std::make_pair(cont1[i],cont2[i]);
            }

            return contOut;
        }  
        
        /**
         * @brief Create new container of size std::max(size_of_container1,size_of_container2) which holds items from both containers at each index
         * 
         * @tparam T type of the containers
         * @param cont1 first container
         * @param cont2 second container
         * @return std::vector<std::pair<T,T> > output container
         */
        template <template<class> typename Container, typename T>
        Container<std::pair<T,T> > zip_to_bigger(const Container<T> &cont1, const Container<T> &cont2)
        {
            std::size_t maxIter = std::max(cont1.size(),cont2.size());
            Container<std::pair<T,T> > contOut(std::pair<T,T>(),maxIter);

            for (std::size_t i = 0; i < maxIter; ++i)
            {
                contOut[i] = std::make_pair(cont1[i],cont2[i]); // this is UB if maxIter > cont1.size() or cont2.size()!
            }

            return contOut;
        }  
    }

#endif