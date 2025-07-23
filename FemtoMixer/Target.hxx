/**
 * @file Target.hxx
 * @author Jędrzej Kołaś (jedrzej.kolas.dokt@pw.edu.pl)
 * @brief Namespace containing information about HADES target setup for each run
 * @version 0.1.0
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#ifndef Target_hxx
    #define Target_hxx
    
    #include <array>

    namespace HADES
    {
        namespace Target
        {
            enum class Setup{Apr12,Feb24Au};

            namespace Detail
            {
                template <Setup T> struct type{};

                static constexpr std::size_t NPlates(type<Setup::Apr12>) noexcept
                {
                    return 15;
                }
                static constexpr std::size_t NPlates(type<Setup::Feb24Au>) noexcept
                {
                    return 15;
                }

                static constexpr std::array<std::pair<double,double>,15> ZPlates(type<Setup::Apr12>) noexcept
                {
                    return {{{-54.7598,0.755846},{-51.6971,0.783591},{-47.7996,0.763387},{-44.5473,0.769386},{-40.569,0.781312},{-37.2151,0.762538},{-33.2948,0.76901},{-30.3726,0.742618},{-26.648,0.748409},{-22.5492,0.738462},{-18.9649,0.747727},{-15.5259,0.749724},{-11.8726,0.740386},{-8.45083,0.742672},{-4.58076,0.712394}}};
                }
                static constexpr std::array<std::pair<double,double>,15> ZPlates(type<Setup::Feb24Au>) noexcept
                {
                    return {{{-61.9605,1.37178},{-58.6207,1.52245},{-54.5652,1.39601},{-50.9504,1.73013},{-46.3113,1.66187},{-44.4911,1.30721},{-39.2215,1.54403},{-36.9806,1.43343},{-33.5816,1.34609},{-28.9379,1.32087},{-25.0631,1.67167},{-22.0411,1.38422},{-18.244,1.39024},{-15.115,1.47253},{-11.479,1.18757}}};
                }

                static constexpr std::pair<double,double> XPlates(type<Setup::Apr12>) noexcept
                {
                    return {0.1951,0.619};
                }
                static constexpr std::pair<double,double> XPlates(type<Setup::Feb24Au>) noexcept
                {
                    return {3.07761,0.936829};
                }

                static constexpr std::pair<double,double> YPlates(type<Setup::Apr12>) noexcept
                {
                    return {0.7045,0.6187};
                }
                static constexpr std::pair<double,double> YPlates(type<Setup::Feb24Au>) noexcept
                {
                    return {-2.17501,1.01304};
                }
            } // namespace Detail

            /**
             * @brief Get the number of plates for a given target setup
             * 
             * @tparam T HADES target setup
             * @return number of plates
             */
            template <Setup T>
            static constexpr std::size_t GetNumberOfPlates() noexcept
            {
                return Detail::NPlates(Detail::type<T>{});
            }
            /**
             * @brief Get the estimated values of mean position and std. dev. in the Z direction for each target plate
             * 
             * @tparam T HADES target setup
             * @return std::array of pairs (first value: mean, second vaule: std. dev.)
             */
            template <Setup T>
            static constexpr std::array<std::pair<double,double>,GetNumberOfPlates<T>() > GetZPlatePositions() noexcept
            {
                return Detail::ZPlates(Detail::type<T>{});
            }
            /**
             * @brief Get the estimated values of mean position and std. dev. of the target in the X direction
             * 
             * @tparam T HADES target setup
             * @return mean (first value) and std. dev. (secon value)
             */
            template <Setup T>
            static constexpr std::pair<double,double> GetXTargetPosition() noexcept
            {
                return Detail::XPlates(Detail::type<T>{});
            }
            /**
             * @brief Get the estimated values of mean position and std. dev. of the target in the Y direction
             * 
             * @tparam T HADES target setup
             * @return mean (first value) and std. dev. (secon value)
             */
            template <Setup T>
            static constexpr std::pair<double,double> GetYTargetPosition() noexcept
            {
                return Detail::YPlates(Detail::type<T>{});
            }

        } // namespace Target
    } // namespace HADES

#endif