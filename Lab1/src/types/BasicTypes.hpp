//
// Created by evgen on 15.10.2024.
//

#ifndef BASICTYPES_HPP
#define BASICTYPES_HPP

#include <complex>
#include <array>

namespace Types {
  using scalar = double;
  using complex = std::complex<double>;
  template<int N>
  using array_s = std::array<scalar, N>;
  template<typename t, int N>
  using array = std::array<t, N>;

}

#endif //BASICTYPES_HPP
