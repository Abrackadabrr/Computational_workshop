//
// Created by evgen on 15.10.2024.
//

#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "Node.hpp"
#include "function_traits/Extractions.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

#include <inmost.h>

namespace Math::Integration {
namespace detail {

/*
 * Возвращает значение точки на ячейке element в зависимости от переданных
 * барицентрических координат
 * @param point -- барицентрические координаты (их n штук, где n -- размерность
 * элемента) всего барицентрических координат будет n+1 штука (по количеству
 * вершин на элементе)
 *
 * @param element -- это элементарная ячейка, по которому происходит
 * интегрирование
 * @note dimention -- это последовательность int
 *                     от 0 до "размерность квадратуры"
 */
template <std::size_t... dimention>
Types::point_t parametrisation(const typename QuadratureTypes::Points::Node<sizeof...(dimention)>::point_t &point,
                               const INMOST::Element &element,
                               const std::index_sequence<dimention...> & /* param to deduce ...dimention */) {
    // количество узлов на каждом элементе завязано на размерность квадратуры
    const int n_nodes{sizeof...(dimention)};
    // вычисление каждой координаты в векторе-параметризации r(l1, ..., ln) через
    // барицентрически координаты (point)
    const Types::scalar x = ((element.getNodes()[dimention].Coords()[0] * point[dimention]) + ...) +
                            element.getNodes()[n_nodes].Coords()[0] * (1 - ((point[dimention]) + ...));
    const Types::scalar y = ((element.getNodes()[dimention].Coords()[1] * point[dimention]) + ...) +
                            element.getNodes()[n_nodes].Coords()[1] * (1 - ((point[dimention]) + ...));
    const Types::scalar z = ((element.getNodes()[dimention].Coords()[2] * point[dimention]) + ...) +
                            element.getNodes()[n_nodes].Coords()[1] * (1 - ((point[dimention]) + ...));
    return {x, y, z};
}

/*
 * Расчитывает сумму по предоставленной квадратуре
 * function = f(r(барицентрические координаты от 1 до n-1))
 * до n - 1 потому что последняя вычисляется тривиально через сумму равной 1
 */
template <typename Quadrature, typename Callable, std::size_t... indexes>
decltype(auto) sub_integrate(const Callable &function_r, std::index_sequence<indexes...> /* no name */) {
    // fold expression instead of for
    return ((Quadrature::nodes[indexes].weight * function_r(Quadrature::nodes[indexes].point)) + ...);
};

template <typename ReturnType> ReturnType measure(const INMOST::Element &element);

template <> inline Types::scalar measure<Types::scalar>(const INMOST::Element &element) {
    switch (element.GetElementType()) {
    case INMOST::EDGE:
        return element.getAsEdge().Length();
    case INMOST::CELL:
        return element.getAsCell().Volume();
    case INMOST::FACE:
        return element.getAsFace().Area();
    default:
        return 0;
    }
}

/*
 * Рассчитывает интеграл по одной заданной параметризации
 * Эта параметризация в барицентрических координатах задается
 * при помощи функции r от n-1 барицентрических коорданат
 */
template <typename Quadrature, typename Parametrization, typename Callable>
typename Extractions::function_info<Callable>::return_type integrate_over_cell(const Parametrization &r,
                                                                               const Callable &function) {
    return sub_integrate<Quadrature>(
        [&](const typename QuadratureTypes::Points::Node<Quadrature::Base::dim>::point_t &point) {
            return function(r(point));
        },
        std::make_index_sequence<Quadrature::nodes.size()>{});
}
} // namespace detail

/**
 * Более удобный вариант предыдущей функции
 */
template <typename Quadrature, typename Callable>
Types::scalar integrate_function_over_cell(const Callable &function, const INMOST::Element &cell) {
    constexpr int dimention = Quadrature::Base::dim;
    const auto r = [&](const typename QuadratureTypes::Points::Node<dimention>::point_t &point) {
        return detail::parametrisation(point, cell, std::make_index_sequence<dimention>{});
    };
    return detail::integrate_over_cell<Quadrature>(r, function) *
        detail::measure<typename Extractions::function_info<Callable>::return_type>(cell);
}

/*
 * Интегрирование скалярной функции по сетке (многообразие)
 */
template <typename Quadrature, typename Callable>
typename Extractions::function_info<Callable>::return_type integrate(const Callable &function, INMOST::Mesh &mesh,
                                                                     const Mesh::CellType &cellType) {
    // размерность квадратуры
    constexpr int dimention = Quadrature::Base::dim;
    typename Extractions::function_info<Callable>::return_type result{};  // zero initialization
  for (auto ielem = mesh.BeginElement(INMOST::EDGE | INMOST::CELL | INMOST::FACE);
       ielem != mesh.EndElement(); ++ielem) {
    if (cellType == ielem->Integer(mesh.GetTag("GMSH_TAGS")) || cellType == Mesh::notag) {
      // параметризация конкретной ячейки
      const auto r = [&](const typename QuadratureTypes::Points::Node<dimention>::point_t &point) {
        return detail::parametrisation(point, ielem->self(), std::make_index_sequence<dimention>{});
      };
      result = result + detail::integrate_over_cell<Quadrature>(r, function) *
        detail::measure<typename Extractions::function_info<Callable>::return_type>(ielem->self());
    }
  }
  return result;
}

}
#endif // QUADRATURE_HPP
