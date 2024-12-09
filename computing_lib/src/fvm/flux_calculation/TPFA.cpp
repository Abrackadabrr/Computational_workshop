//
// Created by evgen on 12.11.2024.
//

#include "TPFA.hpp"
#include "mesh/MeshUtils.hpp"

namespace FVM::TPFA {
Flux::flux<1>
bounbaryNeumanTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                   const std::function<Types::Vector3d(const Types::point_t &)> &gradient_for_boundary_neuman) {
    const auto back_cell = face.BackCell();
    return {{0},
            {back_cell->GlobalID()},
            gradient_for_boundary_neuman(Mesh::Utils::Face::Centroid(face)).dot(Mesh::Utils::Face::UnitNormal(face)) *
                face.Area()};
}

Flux::flux<1> boundaryDirichletTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                                    const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient,
                                    const std::function<Types::scalar(const Types::point_t &)> &dirichlet_function) {
    INMOST::Cell back_cell = face.BackCell();
    const auto centroid1 = Mesh::Utils::Cell::Centroid(back_cell);
    const auto face_centroid = Mesh::Utils::Face::Centroid(face);
    const Types::Vector3d conormal = diffusion_coefficient(face_centroid) * Mesh::Utils::Face::UnitNormal(face);
    const Types::Vector3d r1f = face_centroid - centroid1;
    const Types::scalar l1 = face.Area() * conormal.squaredNorm() / (r1f.dot(conormal));
    return {{l1}, {back_cell.GlobalID()}, -l1 * dirichlet_function(face_centroid)};
}

Flux::flux<2> interiorTPFA(INMOST::Mesh &mesh, INMOST::Face &face,
                           const std::function<Types::Matrix3d(const Types::point_t &)> &diffusion_coefficient) {
    INMOST::Cell back_cell = face.BackCell();
    INMOST::Cell front_cell = face.FrontCell();
    const auto centroid1 = Mesh::Utils::Cell::Centroid(back_cell);
    const auto centroid2 = Mesh::Utils::Cell::Centroid(front_cell);
    const auto face_centroid = Mesh::Utils::Face::Centroid(face);
    const Types::Vector3d conormal = diffusion_coefficient(face_centroid) * Mesh::Utils::Face::UnitNormal(face);
    const Types::Vector3d r1f = face_centroid - centroid1;
    const Types::Vector3d r2f = face_centroid - centroid2;
    const Types::scalar l1 = conormal.squaredNorm() / (r1f.dot(conormal));
    const Types::scalar l2 = conormal.squaredNorm() / (r2f.dot(-conormal));
    const Types::scalar tau = face.Area() * (l1 * l2 / (l1 + l2));
    // std::cout << back_cell.GlobalID() << ' ' << front_cell.GlobalID() <<std::endl;
    return Flux::flux<2>{Types::array_s<2>{tau, -tau},
                         Types::array<int, 2>{back_cell.GlobalID(), front_cell.GlobalID()}, 0};
}
} // namespace FVM::TPFA