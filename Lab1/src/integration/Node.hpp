//
// Created by evgen on 15.10.2024.
//

#ifndef NODE_HPP
#define NODE_HPP

#include "types/BasicTypes.hpp"

namespace Math::Integration::QuadratureTypes::Points {

/*
 * Тип узла интегрирования в барицентрических координатах
 */
template <int dim> struct Node {
  using point_t = Types::array_s<dim>;
  // барицентрические координаты
  point_t point;
  // веса квадратуры в этих координатах
  Types::scalar weight;

  constexpr Node(const Types::array_s<dim> &_point, const Types::scalar &_weight)
      : point(_point), weight(_weight) {}
};

} // namespace Math::Integration::QuadratureTypes::Points
#endif // NODE_HPP
