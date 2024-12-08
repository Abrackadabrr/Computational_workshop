//
// Created by evgen on 07.11.2024.
//

#ifndef SCHEME_HPP
#define SCHEME_HPP

#include "Scheme.hpp"
#include "types/BasicTypes.hpp"

#include <functional>
#include <inmost_mesh.h>

namespace FVM {
/**
 * Расчет вектора части правой части, отвечающей за неоднородность уравнения
 * @param mesh
 * @param rhs
 * @return
 */
Types::VectorXd getRHS(INMOST::Mesh &mesh, const std::function<Types::scalar(const Types::point_t &)> &rhs_function);

/**
 * Расчет матрицы системы с опорой на некоторый метод расчета потоков
 */
Types::SLAE getSLAE(INMOST::Mesh &mesh,
                    const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient,
                    const std::function<Types::scalar(const Types::point_t &)> &rhs_function,
                    const std::function<Types::scalar(const Types::point_t &)> &dirichlet_boundary,
                    const std::function<Types::Vector3d(const Types::point_t &)> &gradient_for_boundary_neuman);

Types::VectorXd getMeshProjection(INMOST::Mesh& mesh, const std::function<Types::scalar(const Types::point_t &)> &function);


}

#endif //SCHEME_HPP
