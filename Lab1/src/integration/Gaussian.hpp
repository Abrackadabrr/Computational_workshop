//
// Created by evgen on 15.10.2024.
//

#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include "integration/Node.hpp"
#include "types/BasicTypes.hpp"

namespace Math::Integration::QuadratureTypes {
template <unsigned int N> struct GaussianPoints {};

template <> struct GaussianPoints<1> {
  static constexpr int dimention = 1;
  static constexpr Types::array<Points::Node<dimention>, 1> nodes = {
    Points::Node<dimention>{{1. / 2}, 1}
  };
};

template<> struct GaussianPoints<2> {
  static constexpr int dimention = 1;
  static constexpr Types::array<Points::Node<dimention>, 2> nodes = {
    Points::Node<dimention>{{0}, 1./2},
    Points::Node<dimention>{{1}, 1./2}
  };
};

template<> struct GaussianPoints<3> {
  static constexpr int dimention = 1;
  static constexpr Types::array<Points::Node<dimention>, 3> nodes = {
    Points::Node<dimention>{{0.}, 1./6},
    Points::Node<dimention>{{1./2}, 2./3},
    Points::Node<dimention>{{1.}, 1./6}
  };
};

} // namespace Math::Integration::QuadratureTypes

#endif // GAUSSIAN_HPP
