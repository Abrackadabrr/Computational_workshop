//
// Created by evgen on 08.12.2024.
//

#ifndef MATRIXASSEMBLY_HPP
#define MATRIXASSEMBLY_HPP

#include "fem/Parametrisation.hpp"
#include "fem/assembly/Integration.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

namespace FEM::Assembly {
namespace detail {

template <unsigned int N> using submatrix = Eigen::Matrix<Types::scalar, N, N>;
template <unsigned int N> using subvector = Eigen::Vector<Types::scalar, N>;

template <typename BFT> Types::index getNumberOfDOF(const Types::mesh_t &mesh);

template <> inline Types::index getNumberOfDOF<FiniteElements::P1>(const Types::mesh_t &mesh) {
    return mesh.NumberOfNodes();
}

/*
 * @tparam BFT for basis function type
 */
template <typename Quadrature, typename BFT, typename Callable_D, typename Callable_algebraic>
submatrix<BFT::n> getSubmatrix(const Types::cell_t &cell, const Callable_D &D, const Callable_algebraic &c) {
    const Types::index n_Dof = BFT::n;
    submatrix<n_Dof> submatrix{};
    const Types::Matrix3d &A_inverse = Parametrisation::getLocalParametrizationMatrix(cell).inverse();

    // расчет интегралов от функций

    const auto tensor = [&A_inverse, &D](const Types::point_t &x) -> Types::Matrix3d {
        return A_inverse * D(x) * A_inverse.transpose();
    };

    for (Types::index i = 0; i < n_Dof; i++) {
        for (Types::index j = 0; j < n_Dof; j++) {

            const auto scalar_weight = [&](const FiniteElements::Barycentric &b) {
                return BFT::getFunction(i)(b) * BFT::getFunction(j)(b);
            };

            submatrix(i, j) = Integration::integrate_bilinear_form_over_cell<Quadrature>(
                                  cell, tensor, BFT::getGradient(i), BFT::getGradient(j)) +
                              Integration::integrate_over_cell_with_weight<Quadrature>(cell, c, scalar_weight);
        }
    }
    return submatrix;
}

/*
 * @tparam BFT тип базисных функций
 */
template <typename Quadrature, typename BFT, typename Callable>
subvector<BFT::n> getSubrhs(const Types::cell_t &cell, const Callable &rhs) {
    subvector<BFT::n> subrhs{};
    for (int i = 0; i < BFT::n; i++) {
        subrhs[i] = Integration::integrate_over_cell_with_weight<Quadrature>(cell, rhs, BFT::getFunction(i));
    }
    return subrhs;
}
} // namespace detail

template <typename Quadrature, typename BFT, typename d_tensor_t, typename algebraic_t>
Types::SparseMatrixXd getMatrix(Types::mesh_t &mesh, const d_tensor_t &D, const algebraic_t &c) {
    // триплеты для задания разреженной матрицы
    Types::vector<Eigen::Triplet<Types::scalar>> triplets;
    triplets.reserve(mesh.NumberOfCells() * BFT::n * BFT::n);
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto &global_indexes = Parametrisation::getGlobalIndexesOfDOF<BFT>(ielem->self());
        const auto &A_sub = detail::getSubmatrix<Quadrature, BFT>(ielem->self(), D, c);
        for (int i = 0; i < BFT::n; i++) {
            for (int k = 0; k < BFT::n; k++) {
                triplets.push_back(Eigen::Triplet<Types::scalar>(global_indexes[i], global_indexes[k], A_sub(i, k)));
            }
        }
    }
    Types::SparseMatrixXd result{detail::getNumberOfDOF<BFT>(mesh), detail::getNumberOfDOF<BFT>(mesh)};
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

template <typename Quadrature, typename BFT, typename rsh_t>
Types::VectorXd getRhs(Types::mesh_t &mesh, const rsh_t &rhs) {
    Types::VectorXd result = Types::VectorXd::Zero(detail::getNumberOfDOF<BFT>(mesh));
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto &global_indexes = Parametrisation::getGlobalIndexesOfDOF<BFT>(ielem->self());
        const auto &rhs_sub = detail::getSubrhs<Quadrature, BFT>(ielem->self(), rhs);
        for (int i = 0, j = 0; i < BFT::n; i++) {
            result[global_indexes[i]] += rhs_sub(i);
        }
    }
    return result;
}

template <typename Quadrature, typename BFT, typename rsh_t, typename d_tensor_t, typename algebraic_t>
Types::SLAE getSLAE(Types::mesh_t &mesh, const rsh_t &rhs, const d_tensor_t &D, const algebraic_t &c) {
    const auto &matrix = getMatrix<Quadrature, BFT>(mesh, D, c);
    const auto &rhs_res = getRhs<Quadrature, BFT>(mesh, rhs);
    return {matrix, rhs_res};
}

} // namespace FEM::Assembly

#endif //MATRIXASSEMBLY_HPP
