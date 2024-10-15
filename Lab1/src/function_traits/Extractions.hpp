//
// Created by evgen on 15.10.2024.
//

#ifndef EXTRACTIONS_HPP
#define EXTRACTIONS_HPP
#include <tuple>

namespace Extractions {

template<typename Type>
struct function_info{};

template<typename ReturnType, typename... Args>
struct function_info<ReturnType(Args...)> {
  using return_type = ReturnType;
  using args_tuple = std::tuple<Args...>;
  static constexpr unsigned dimention = sizeof...(Args);
};

}

#endif //EXTRACTIONS_HPP
