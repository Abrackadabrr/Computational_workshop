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

enum BoundaryType { NOT_BOUNDARY = 0, DIRICHLET = 1, NEUMANN = 2, ROBIN = 3 };

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

template<typename Callable>
void assign_analytical_to_nodes(const Callable &f, INMOST::Mesh &mesh, const std::string &name) {
    INMOST::Tag u_tag = mesh.CreateTag(name, INMOST::DATA_REAL, INMOST::NODE, INMOST::NONE, 1);
    for (auto ci = mesh.BeginNode(); ci != mesh.EndNode(); ++ci)
        ci->Real(u_tag) = f(Mesh::Utils::getPoint(ci->self()));
}

template <typename Callable>
void assign_analytical_to_cells(const Callable &func, INMOST::Mesh &mesh, const std::string &name) {
    INMOST::Tag u_tag = mesh.CreateTag(name, INMOST::DATA_REAL, INMOST::CELL, INMOST::NONE, 1);
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci)
        ci->Real(u_tag) = func(Cell::Centroid(ci->self()));
}

Types::scalar norm_over_mesh(const Types::VectorXd &scalarField, INMOST::Mesh &mesh);

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

template <typename Callable> Types::VectorXd evaluate_over_nodes(const Callable &func, INMOST::Mesh &mesh) {
    const int n = mesh.NumberOfNodes();
    Types::VectorXd result{n};
    int i = 0;
    for (auto ci = mesh.BeginNode(); ci != mesh.EndNode(); ++ci)
        result[i++] = func(ci->self());
    return result;
}

template<typename ElementType, typename InmostAwfulArray>
Types::vector<ElementType> getArray(const InmostAwfulArray &array) {
    Types::vector<ElementType> array_data{};
    for (auto ci = array.begin(); ci != array.end(); ++ci) {
        array_data.push_back(ci->self());
    }
    return array_data;
}

} // namespace Utils
} // namespace Mesh

#endif // MESHUTILS_HPP
