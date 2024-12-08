//
// Created by evgen on 13.11.2024.
//

#include "MeshUtils.hpp"
#include "inmost.h"

namespace Mesh::Utils {
void fill(const Types::VectorXd &scalarField, INMOST::Mesh &mesh, const std::string &name) {
    INMOST::Tag u_tag = mesh.CreateTag(name, INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci)
        ci->Real(u_tag) = scalarField[ci->GlobalID()];
}

Types::scalar norm_over_mesh(const Types::VectorXd &scalarField, INMOST::Mesh &mesh) {
    Types::scalar result = 0;
    for (auto ci = mesh.BeginCell(), end = mesh.EndCell(); ci != end; ++ci) {
        const auto index = ci->GlobalID();
        result += scalarField(index) * scalarField(index) * ci->Volume();
    }
    return std::sqrt(result);
}

} // namespace Mesh::Utils
