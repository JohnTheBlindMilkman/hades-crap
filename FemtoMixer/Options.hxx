
/*
Activity Indicators for Modern C++
https://github.com/p-ranav/indicators

Licensed under the MIT License <http://opensource.org/licenses/MIT>.
SPDX-License-Identifier: MIT
Copyright (c) 2019 Dawid Pilarski <dawid.pilarski@panicsoftware.com>.

Permission is hereby  granted, free of charge, to any  person obtaining a copy
of this software and associated  documentation files (the "Software"), to deal
in the Software  without restriction, including without  limitation the rights
to  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell
copies  of  the Software,  and  to  permit persons  to  whom  the Software  is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR
IMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,
FITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE
AUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER
LIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/* 
    Adapted by Jedrzej Kolas <jedrzej.kolas.dokt@pw.edu.pl>
    Date: 08-12-2023
*/
#ifndef INDICATORS_SETTING
#define INDICATORS_SETTING

#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace FemtoCorrelation 
{

    namespace details 
    {

        template <bool condition> struct if_else;

        template <> struct if_else<true> { using type = std::true_type; };

        template <> struct if_else<false> { using type = std::false_type; };

        template <bool condition, typename True, typename False> struct if_else_type;

        template <typename True, typename False> struct if_else_type<true, True, False> { using type = True; };

        template <typename True, typename False> struct if_else_type<false, True, False> { using type = False; };

        template <typename... Ops> struct conjuction;

        template <> struct conjuction<> : std::true_type {};

        template <typename Op, typename... TailOps>
        struct conjuction<Op, TailOps...> : if_else_type<!Op::value, std::false_type, conjuction<TailOps...>>::type {};

        template <typename... Ops> struct disjunction;

        template <> struct disjunction<> : std::false_type {};

        template <typename Op, typename... TailOps>
        struct disjunction<Op, TailOps...> : if_else_type<Op::value, std::true_type, disjunction<TailOps...>>::type {};

        enum class MixingOption 
        {
            events_to_mix = 0,
            kT_differential,
            rapidity_differential,
            azimuthally_differential,
            centrality_differential
        };

        template <typename T, MixingOption Id> struct Setting 
        {
            template <typename... Args, typename = typename std::enable_if<std::is_constructible<T, Args...>::value>::type>
            explicit Setting(Args &&... args) : value(std::forward<Args>(args)...) {}
            Setting(const Setting &) = default;
            Setting(Setting &&) = default;

            static constexpr auto id = Id;
            using type = T;

            T value{};
        };

        template <typename T> struct is_setting : std::false_type {};

        template <MixingOption Id, typename T> struct is_setting<Setting<T, Id>> : std::true_type {};

        template <typename... Args>
        struct are_settings : if_else<conjuction<is_setting<Args>...>::value>::type {};

        template <> struct are_settings<> : std::true_type {};

        template <typename Setting, typename Tuple> struct is_setting_from_tuple;

        template <typename Setting> struct is_setting_from_tuple<Setting, std::tuple<>> : std::true_type {};

        template <typename Setting, typename... TupleTypes>
        struct is_setting_from_tuple<Setting, std::tuple<TupleTypes...>> : if_else<disjunction<std::is_same<Setting, TupleTypes>...>::value>::type {};

        template <typename Tuple, typename... Settings>
        struct are_settings_from_tuple : if_else<conjuction<is_setting_from_tuple<Settings, Tuple>...>::value>::type {};

        template <MixingOption Id> struct always_true { static constexpr auto value = true; };

        template <MixingOption Id, typename Default> Default &&get_impl(Default &&def) { return std::forward<Default>(def); }

        template <MixingOption Id, typename Default, typename T, typename... Args>
        auto get_impl(Default && /*def*/, T &&first, Args &&... /*tail*/) ->
            typename std::enable_if<(std::decay<T>::type::id == Id), decltype(std::forward<T>(first))>::type { return std::forward<T>(first); }

        template <MixingOption Id, typename Default, typename T, typename... Args>
        auto get_impl(Default &&def, T && /*first*/, Args &&... tail) ->
            typename std::enable_if<(std::decay<T>::type::id != Id), decltype(get_impl<Id>(std::forward<Default>(def), std::forward<Args>(tail)...))>::type 
            { 
                return get_impl<Id>(std::forward<Default>(def), std::forward<Args>(tail)...); 
            }

        template <MixingOption Id, typename Default, typename... Args, typename = typename std::enable_if<are_settings<Args...>::value, void>::type>
        auto get(Default &&def, Args &&... args) -> 
            decltype(details::get_impl<Id>(std::forward<Default>(def), std::forward<Args>(args)...)) 
            {
                return details::get_impl<Id>(std::forward<Default>(def), std::forward<Args>(args)...);
            }

        template <MixingOption Id> using StringSetting = Setting<std::string, Id>;

        template <MixingOption Id> using FloatSetting = Setting<float, Id>;

        template <MixingOption Id> using IntegerSetting = Setting<std::size_t, Id>;

        template <MixingOption Id> using BooleanSetting = Setting<bool, Id>;

        template <MixingOption Id, typename Tuple, std::size_t counter = 0> struct option_idx;

        template <MixingOption Id, typename T, typename... Settings, std::size_t counter>
        struct option_idx<Id, std::tuple<T, Settings...>, counter> : if_else_type<(Id == T::id), std::integral_constant<std::size_t, counter>, option_idx<Id, std::tuple<Settings...>, counter + 1>>::type {};

        template <MixingOption Id, std::size_t counter> struct option_idx<Id, std::tuple<>, counter> 
        {
            static_assert(always_true<(MixingOption)Id>::value, "No such option was found");
        };

        template <MixingOption Id, typename Settings>
        auto get_value(Settings &&settings) -> 
            decltype((std::get<option_idx<Id, typename std::decay<Settings>::type>::value>(std::declval<Settings &&>()))) 
            {
                return std::get<option_idx<Id, typename std::decay<Settings>::type>::value>(std::forward<Settings>(settings));
            }

    } // namespace details

    namespace option 
    {
        using EventsToMix = details::IntegerSetting<details::MixingOption::events_to_mix>;
        using KtDifferential = details::Setting<std::vector<std::pair<float,float> >,details::MixingOption::kT_differential>;
        using RapidityDifferential = details::Setting<std::vector<std::pair<float,float> >,details::MixingOption::rapidity_differential>;
        using AzimuthallyDifferential = details::Setting<std::vector<std::pair<float,float> >,details::MixingOption::azimuthally_differential>;
        using CentralityDifferential = details::Setting<std::vector<std::pair<int,int> >,details::MixingOption::centrality_differential>;
    } // namespace option
} // namespace indicators

#endif