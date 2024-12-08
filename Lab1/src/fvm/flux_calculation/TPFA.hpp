//
// Created by evgen on 12.11.2024.
//

#ifndef TPFA_HPP
#define TPFA_HPP

#include "FluxCalculationMetod.hpp"
#include "types/BasicTypes.hpp"

#include "inmost.h"

namespace FVM::TPFA {

Flux::flux<1>
bounbaryNeumanTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                   const std::function<Types::Vector3d(const Types::point_t &)> &gradient_for_boundary_neuman);

Flux::flux<2> interiorTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                           const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient);

Flux::flux<1> boundaryDirichletTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                                    const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient,
                                    const std::function<Types::scalar(const Types::point_t &)> &dirichlet_function);

} // namespace FVM::TPFA

#endif // TPFA_HPP
