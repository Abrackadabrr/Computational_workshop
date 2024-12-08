//
// Created by evgen on 07.11.2024.
//

#include "Scheme.hpp"
#include "flux_calculation/TPFA.hpp"
#include "integration/Quadrature.hpp"
#include "integration/Gaussian.hpp"

namespace FVM {
/**
 * Расчет вектора части правой части, отвечающей за неоднородность уравнения
 * @param mesh
 * @param rhs
 * @return
 */
Types::VectorXd getRHS(INMOST::Mesh &mesh, const std::function<Types::scalar(const Types::point_t &)> &rhs_function) {
    const int cells_n = mesh.NumberOfCells();
    Types::VectorXd rhs{cells_n};
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        // интеграл методом средней точки
        rhs[ielem->GlobalID()] = rhs_function(Mesh::Utils::Cell::Centroid(ielem->self())) * ielem->Volume();
        // Math::Integration::detail::integrate_function_over_cell<
        // Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(rhs_function, ielem->self());
    }
    return rhs;
}

Types::SLAE getSLAE(INMOST::Mesh &mesh,
                    const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient,
                    const std::function<Types::scalar(const Types::point_t &)> &rhs_function,
                    const std::function<Types::scalar(const Types::point_t &)> &dirichlet_boundary,
                    const std::function<Types::Vector3d(const Types::point_t &)> &gradient_for_boundary_neuman) {
    Types::VectorXd b = getRHS(mesh, rhs_function);
    Types::vector<Eigen::Triplet<Types::scalar>> triplets;
    triplets.reserve(mesh.NumberOfFaces());

    for (auto ielem = mesh.BeginFace(), end = mesh.EndFace(); ielem != end; ++ielem) {
        if (!ielem->Boundary()) {
            const auto flux_data = TPFA::interiorTPFA(mesh, ielem->self(), diffusion_coefficient);
            // std::cout << flux_data.coef[0] << ' ' << flux_data.coef[1] << std::endl;
            const int back_cell_index = flux_data.indexes[0];
            const int front_cell_index = flux_data.indexes[1];
            // как уравнение, записанное для back_cell ячейки
            triplets.emplace_back(back_cell_index, back_cell_index, flux_data.coef[0]);
            triplets.emplace_back(back_cell_index, front_cell_index, flux_data.coef[1]);
            // как уравнение, записанное для front_cell ячейки
            triplets.emplace_back(front_cell_index, back_cell_index, -flux_data.coef[0]);
            triplets.emplace_back(front_cell_index, front_cell_index, -flux_data.coef[1]);
        } else {
            // для граничных элементов
            // условия дирихле
            if (ielem->Integer(mesh.GetTag("Boundary_type")) == Mesh::BoundaryType::DIRICHLET) {
                const auto flux_data =
                    TPFA::boundaryDirichletTPFA(mesh, ielem->self(), diffusion_coefficient, dirichlet_boundary);

                b[flux_data.indexes[0]] -= flux_data.b;
                triplets.emplace_back(flux_data.indexes[0], flux_data.indexes[0], flux_data.coef[0]);
            }
            // условие Неймана
            if (ielem->Integer(mesh.GetTag("Boundary_type")) == Mesh::BoundaryType::NEUMANN) {
                const auto flux_data = TPFA::bounbaryNeumanTPFA(mesh, ielem->self(), gradient_for_boundary_neuman);
                b[flux_data.indexes[0]] -= flux_data.b;
            }
            // условие Робэна
        }
    }
    Types::SparseMatrixXd A{mesh.NumberOfCells(), mesh.NumberOfCells()};
    A.setFromTriplets(triplets.begin(), triplets.end());
    return {A, b};
}

Types::VectorXd getMeshProjection(INMOST::Mesh &mesh,
                                  const std::function<Types::scalar(const Types::point_t &)> &function) {
    Types::VectorXd res{mesh.NumberOfCells()};
    for (auto ci = mesh.BeginCell(), end = mesh.EndCell(); ci != end; ++ci) {
        res[ci->GlobalID()] = function(Mesh::Utils::Cell::Centroid(ci->self()));
    }
    return res;
}

}
