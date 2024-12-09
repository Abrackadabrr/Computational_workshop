//
// Created by evgen on 08.12.2024.
//

#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include "function_traits/Extractions.hpp"
#include "integration/Node.hpp"
#include "integration/Quadrature.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

namespace FEM::Integration {
namespace detail {
/*
 * Функция, интегрирующая с весом в барицентрических координатах
 * При этом функция веса считается заданной в барицентрических координатах
 * Это необходимо для того, чтобы реализовать так называемый "расчет на референсном элементе"
 */

template <typename Quadrature, typename Parametrization, typename Callable, typename WeightInL>
decltype(auto) integrate_over_cell_with_weight(const Parametrization &r, const Callable &function, const WeightInL &w) {
    // здесь параметризовывается функция function с помощью замены r: lambda -> x
    // при этом к весу w замена не применяется
    return Math::Integration::detail::sub_integrate<Quadrature>(
        [&](const typename Math::Integration::QuadratureTypes::Points::Node<Quadrature::Base::dim>::point_t &point) {
            const auto f_result = function(r(point));
            const auto w_result = w(point);
            return f_result * w_result;
        },
        std::make_index_sequence<Quadrature::nodes.size()>{});
}

/*
 * Функция, интегрирующая билинейную форму тензора D <w2, Dw1> в барицентрических координатах
 * При этом функция веса считается заданной в барицентрических координатах
 * Это необходимо для того, чтобы реализовать так называемый "расчет на референсном элементе"
 * @param w1 векторнозначная функция от барицентрических координат.
 * @param w2 векторнозначная функция от барицентрических координат.
 * @param D матричнозначная функия от декакртовых координат
 * @prama r парматеризация ячейки с помощью бирицентрических координат
 */
template <typename Quadrature, typename Parametrization, typename Callable, typename W1InL, typename W2InL>
decltype(auto) integrate_bilinear_form_over_cell(const Parametrization &r, const Callable &D, const W1InL &w1,
                                                 const W2InL &w2) {
    // здесь параметризовывается функция function с помощью замены r: lambda -> x
    // при этом к весу w замена не применяется
    return Math::Integration::detail::sub_integrate<Quadrature>(
        [&](const typename Math::Integration::QuadratureTypes::Points::Node<Quadrature::Base::dim>::point_t &point) {
            const auto tensor_result = D(r(point));
            const auto w1_result = w1(point);
            const auto w2_result = w2(point);
            return (tensor_result * w1_result).dot(w2_result);
        },
        std::make_index_sequence<Quadrature::nodes.size()>{});
}
} // namespace detail

/**
 * Интегрирование скалярных функиций с весом по ячейке, при этом вес задан с барицентрических координатах ячейки
 *
 * Эта функция позволяет расчитать локальный вектор правой части и локальную матрицу масс
 */
template <typename Quadrature, typename Callable, typename WeightInL>
Types::scalar integrate_over_cell_with_weight(const INMOST::Element &cell, const Callable &function,
                                              const WeightInL &w) {
    constexpr int dimention = Quadrature::Base::dim;
    // здесь создается функция для замены переменных x <-> lambda
    const auto r = [&](const typename Math::Integration::QuadratureTypes::Points::Node<dimention>::point_t &point) {
        return Math::Integration::detail::parametrisation(point, cell, std::make_index_sequence<dimention>{});
    };
    return detail::integrate_over_cell_with_weight<Quadrature>(r, function, w) *
           Math::Integration::detail::measure<Types::scalar>(cell);
}

/**
 * Интегрирование градиентной билинейной формы тензора D <w2, Dw1> по ячейке, при этом вес задан
 * в барицентрических координатах ячейки
 * @param w1 векторнозначная функция от барицентрических координат.
 * @param w2 векторнозначная функция от барицентрических координат.
 * @param D матричнозначная функия от декакртовых координат
 * Эта функция позволяет расчитать локальный вектор правой части и локальную матрицу масс
 */
template <typename Quadrature, typename Callable, typename W1InL, typename W2InL>
Types::scalar integrate_bilinear_form_over_cell(const INMOST::Element &cell, const Callable &function,
                                                          const W1InL &w1, const W2InL &w2) {
    constexpr int dimention = Quadrature::Base::dim;
    // здесь создается функция для замены переменных x <-> lambda
    const auto r = [&](const typename Math::Integration::QuadratureTypes::Points::Node<dimention>::point_t &point) {
        return  Math::Integration::detail::parametrisation(point, cell, std::make_index_sequence<dimention>{});
    };
    return detail::integrate_bilinear_form_over_cell<Quadrature>(r, function, w1, w2) *
            Math::Integration::detail::measure<Types::scalar>(cell);
}

}

#endif //INTEGRATION_HPP
