//
// Created by evgen on 15.10.2024.
//

#ifndef EXTRACTIONS_HPP
#define EXTRACTIONS_HPP
#include <functional>
#include <tuple>

namespace Extractions {

template<typename F>
struct function_info{};

template <typename ReturnType, typename... Args> struct function_info<ReturnType(Args...)> {
    using return_type = ReturnType;
    using args_tuple = std::tuple<Args...>;
    static constexpr unsigned dimention = sizeof...(Args);
};

template <typename ReturnType, typename... Args> struct function_info<std::function<ReturnType(Args...)>> {
    using return_type = ReturnType;
    using args_tuple = std::tuple<Args...>;
    static constexpr unsigned dimention = sizeof...(Args);
};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_info<ReturnType (ClassType::*)(Args...) const> {
    using result_type = ReturnType;
    using arg_tuple = std::tuple<Args...>;
    static constexpr auto dimention = sizeof...(Args);
};
} // namespace Extractions

#endif // EXTRACTIONS_HPP
