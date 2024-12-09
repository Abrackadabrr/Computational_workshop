//
// Created by evgen on 08.12.2024.
//

#ifndef MATRIXASSEMBLY_HPP
#define MATRIXASSEMBLY_HPP

#include "fem/Parametrisation.hpp"
#include "fem/assembly/Integration.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

namespace FEM {

namespace detail {

template <unsigned int N> using submatrix = Eigen::Matrix<Types::scalar, N, N>;
template <unsigned int N> using subvector = Eigen::Vector<Types::scalar, N>;

/*
 * @tparam BFT for basis function type
 */
template <typename Quadrature, typename BFT, typename Callable_D>
submatrix<BFT::n> getSubmatrix(const Types::cell_t &cell, const Callable_D &D) {
    const Types::index n_Dof = BFT::n;
    submatrix<n_Dof> submatrix{};
    const Types::Matrix3d &A_inverse = getLocalParametrizationMatrix(cell).inverse();

    // расчет интегралов от функций

    const auto tensor = [&A_inverse, &D](const Types::point_t &x) -> Types::Matrix3d {
        return A_inverse * D(x) * A_inverse.transpose();
    };

    for (Types::index i = 0; i < n_Dof; i++) {
        for (Types::index j = 0; j < n_Dof; j++) {
            submatrix(i, j) =
                Integration::integrate_bilinear_form_over_cell<Quadrature>(cell, tensor, BFT::getGradient(i), BFT::getGradient(j));
        }
    }
    return submatrix;
}

/*
 * @tparam BFT for basis function type
 */
template <typename Quadrature, typename BFT, typename Callable>
subvector<BFT::n> getSubrhs(const Types::cell_t &cell, const Callable &rhs) {
    const subvector<BFT::n> subrhs{};
    for (int i = 0; i < BFT::n; i++) {
        subrhs[i] = Integration::integrate_over_cell_with_weight<Quadrature>(cell, rhs, BFT::getFunction(i));
    }
    return subrhs;
} // namespace detail

} // namespace FEM

#endif //MATRIXASSEMBLY_HPP

