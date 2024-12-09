#ifndef MESHUTILS_HPP
#define MESHUTILS_HPP

#include "inmost.h"
#include "types/BasicTypes.hpp"

namespace Mesh {

enum CellType {
    cut = 1,
    endocard = 2,
    epicard = 3,
    notag = -1,
};

enum BoundaryType { DIRICHLET = 1, NEUMANN = 2, ROBIN = 3 };

#define CREATE_EIGEN_COVER(Element, Fucntion)                                                                          \
    inline Types::Vector3d Fucntion(const INMOST::Element &el) {                                                       \
        auto *value = new Types::scalar[3];                                                                            \
        el->Fucntion(value);                                                                                           \
        return {value[0], value[1], value[2]};                                                                         \
    }

namespace Utils {
namespace Face {
CREATE_EIGEN_COVER(Face, UnitNormal)
CREATE_EIGEN_COVER(Face, Centroid)
} // namespace Face

namespace Cell {
CREATE_EIGEN_COVER(Cell, Centroid)
}

Types::point_t inline getPoint(const INMOST::Node &node) {
    return {node.Coords()[0], node.Coords()[1], node.Coords()[2]};
}

void assign_to_nodes(const Types::VectorXd &scalarField, INMOST::Mesh &mesh, const std::string &name);

void assign_to_cells(const Types::VectorXd &scalarField, INMOST::Mesh &mesh, const std::string &name);

Types::scalar norm_over_mesh(const Types::VectorXd &scalarField, INMOST::Mesh &mesh);

template <typename Callable>
void fill_scalar_with_function(const Callable &func, INMOST::Mesh &mesh, const std::string &name) {
    INMOST::Tag u_tag = mesh.CreateTag(name, INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci)
        ci->Real(u_tag) = func(Cell::Centroid(ci->self()));
}

template <typename Callable>
void fill_vector3d_with_function(const Callable &func, INMOST::Mesh &mesh, const std::string &name) {
    INMOST::Tag u_tag = mesh.CreateTag(name, INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 3);
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci)
        ci->RealArray(u_tag) = func(Cell::Centroid(ci->self())).data();
}

template <typename Callable> Types::VectorXd evaluate_function_over_cells(const Callable &func, INMOST::Mesh &mesh) {
    const int cells_n = mesh.NumberOfCells();
    Types::VectorXd result{cells_n};
    int i = 0;
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci)
        result[i++] = func(Cell::Centroid(ci->self()));
    return result;
}

} // namespace Utils
} // namespace Mesh

#endif // MESHUTILS_HPP
