//
// Created by evgen on 15.10.2024.
//

#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "Node.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

#include <inmost.h>

namespace Math::Integration {
namespace detail {

/*
 * Возвращает значение точки на ячейке element в зависимости от переданных
 * барицентрических координат
 * @param dimention -- это последовательность int
 *                     от 0 до "размерность квадратуры - 1"
 */
template <std::size_t... dimention>
Types::array<Types::scalar, 3> parametrisation(
    const typename QuadratureTypes::Points::Node<sizeof...(dimention)>::point_t
        &point,
    const INMOST::Element &element,
    const std::index_sequence<dimention...>
        & /* param to deduce ...dimention */) {
  // количество узлов на каждом элементе завязано на размерность квадратуры
  const int n_nodes{sizeof...(dimention)};
  // вычисление каждой координаты в векторе-параметризации: r(l1, ..., ln)
  const Types::scalar x =
      ((element.getNodes()[dimention].Coords()[0] * point[dimention]) + ...) +
      element.getNodes()[n_nodes].Coords()[0] *
          (1 - ((point[dimention]) + ...));
  const Types::scalar y =
      ((element.getNodes()[dimention].Coords()[1] * point[dimention]) + ...) +
      element.getNodes()[n_nodes].Coords()[1] *
          (1 - ((point[dimention]) + ...));
  const Types::scalar z =
      ((element.getNodes()[dimention].Coords()[2] * point[dimention]) + ...) +
      element.getNodes()[n_nodes].Coords()[1] *
          (1 - ((point[dimention]) + ...));
  return {x, y, z};
}

/*
 * Расчитывает сумму по предоставленной квадратуре
 * function = f(r(барицентрические координаты от 1 до n-1))
 * до n - 1 потому что последняя вычисляется тривиально через сумму равной 1
 */
template <typename Quadrature, typename Callable, std::size_t... indexes>
Types::scalar sub_integrate(const Callable &function_r,
                            std::index_sequence<indexes...> /* no name */) {
  // fold expression instead of for
  return ((Quadrature::nodes[indexes].weight *
           function_r(Quadrature::nodes[indexes].point)) +
          ...);
};

/*
 * Рассчитывает интеграл по одной заданной параметризации
 * Эта параметризация в барицентрических координатах задается
 * при помощи функции r от n-1 барицентрических коорданат
 */
template <typename Quadrature, typename Parametrization, typename Callable>
Types::scalar integrate_over_cell(const Parametrization &r,
                                  const Callable &function) {
  return sub_integrate<Quadrature>(
      [&](const typename QuadratureTypes::Points::Node<
          Quadrature::dimention>::point_t &point) {
        return function(r(point));
      }, std::make_index_sequence<Quadrature::nodes.size()>{});
};
} // namespace detail

/*
 * Интегрирование скалярной функции по сетке (многообразие)
 */
template <typename Quadrature, typename Callable>
Types::scalar integrate(const Callable &function, INMOST::Mesh &mesh,
                        const Mesh::CellType &cellType) {
  // размерность квадратуры
  constexpr int dimention = Quadrature::dimention;
  Types::scalar result = 0;
  for (auto ielem = mesh.BeginElement(INMOST::CELL | INMOST::FACE);
       ielem != mesh.EndElement(); ++ielem) {
    if (cellType == ielem->Integer(mesh.GetTag("GMSH_TAGS"))) {
      // параметризация конкретной ячейки
      const auto r = [&](const typename QuadratureTypes::Points::Node<dimention>::point_t &point) {
        return detail::parametrisation(point, ielem->self(), std::make_index_sequence<dimention>{});
      };
      result += detail::integrate_over_cell<Quadrature>(r, function);
    }
  }
  return result;
}

}
#endif // QUADRATURE_HPP
