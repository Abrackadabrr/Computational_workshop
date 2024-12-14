//
// Created by evgen on 14.12.2024.
//

#ifndef IMPLICITEULER_HPP
#define IMPLICITEULER_HPP
#include <types/BasicTypes.hpp>
#include <types/MeshTypes.hpp>
#include "fem/assembly/MatrixAssembly.hpp"
#include "fem/time_discretisation/TimeDiscretisation.hpp"
#include "auxilary/SchemeParameters.hpp"

namespace FEM::Time::Schemes {

template<typename BFT, typename initial_t>
Types::VectorXd getInitialCondition(Types::mesh_t& mesh, const initial_t& initial_cond) {
    Types::VectorXd result; result.resize(FiniteElements::getNumberOfDOF<BFT>(mesh));

    for (auto i_cell = mesh.BeginCell(); i_cell != mesh.EndCell(); ++i_cell) {
        const auto local_dof = Assembly::Utils::getLocalDOF<BFT>(i_cell->self());
        const auto global_indexes_of_dof = Assembly::Utils::getGlobalIndexesOfDOF<BFT>(i_cell->self());

        for (int i = 0; i < local_dof.size(); ++i) {
            result[global_indexes_of_dof[i]] = initial_cond(Mesh::Utils::getPoint(local_dof[i]));
        }
    }
    return result;
}

template<typename MethodSLAE, typename Quadrature, typename BFT, typename initial_t, typename rsh_t, typename d_tensor_t, typename algebraic_t,
          typename neumann_t, typename dirichlet_t>
Types::vector<Types::VectorXd> solveDynamicProblem(Types::mesh_t& mesh,
                                                   const Aux::TimeIntegrationParameters integration_params,
                                                   const initial_t& initial_cond,
                                                   const rsh_t &rhs,
                                                   const d_tensor_t &D,
                                                   const algebraic_t &c,
                                                   const neumann_t &neumann,
                                                   const dirichlet_t &dirichlet) {
    Types::vector<Types::VectorXd> result; result.reserve(integration_params.N);
    const Types::SLAE slae = Assembly::getSLAE<Quadrature, BFT>(mesh, rhs, D, c, neumann, dirichlet);
    const Types::VectorXd f = slae.b * integration_params.timeStep;
    const Types::SparseMatrixXd time_matrix = Time::getMatrix<Quadrature, BFT>(mesh);
    const Types::SparseMatrixXd matrix = integration_params.timeStep * slae.A + time_matrix;

    Types::VectorXd current_state = getInitialCondition<BFT>(mesh, initial_cond);
    result.push_back(current_state);

    MethodSLAE method;
    method.compute(matrix);
    method.setTolerance(1e-10);
    std::cout << "Time Integration Starts" << std::endl;
    for (int i = 0; i < integration_params.N; i++) {
        const Types::VectorXd b = time_matrix * current_state + f;
        current_state = method.solve(b);
        result.push_back(current_state);
        if (i % 50 == 0) {
            std::cout << i << " done" << std::endl;
        }
    }
    return result;
}

}

#endif //IMPLICITEULER_HPP
