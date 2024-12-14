//
// Created by evgen on 09.12.2024.
//

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCholesky"

// fem
#include "fem/assembly/MatrixAssembly.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "fem/time_discretisation/ImplicitEuler.hpp"
#include "integration/Gaussian.hpp"

// test
#include "../Lab2/TestFucntions.hpp"

// mesh
#include "mesh/MeshUtils.hpp"

#include <ranges>
#include <filesystem>


std::vector<unsigned> set_boundary_tag(INMOST::Mesh &mesh) {
    INMOST::Tag u_tag = mesh.CreateTag("Boundary_type", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::NONE, 1);
    std::vector<unsigned> boundary_nodes;
    const auto &gmsh_tag = mesh.GetTag("GMSH_TAGS");
    for (auto ci = mesh.BeginFace(); ci != mesh.EndFace(); ++ci) {
        auto face = ci->getAsFace();
        if (face.Boundary()) {
            const auto nodes = face.getNodes();
            for (int i = 0; i < nodes.size(); ++i) {
                boundary_nodes.push_back(nodes[i].GlobalID());
            }
            if (face.Integer(gmsh_tag) == Mesh::cut)
                face.Integer(u_tag) = Mesh::BoundaryType::DIRICHLET;
            else
                face.Integer(u_tag) = Mesh::BoundaryType::NEUMANN;
        } else {
            face.Integer(u_tag) = Mesh::BoundaryType::NOT_BOUNDARY;
        }
    }
    return boundary_nodes;
}

Types::scalar evaluate_mean_diameter(INMOST::Mesh &mesh) {
    Types::scalar result = 0;
    for (auto ci = mesh.BeginCell(); ci != mesh.EndCell(); ++ci) {
        const auto faces = ci->getFaces();
        for (auto fi = faces.begin(); fi != faces.end(); ++fi) {
            Types::scalar local_res = ci->Volume() * 3 / fi->Area();
            result = std::max(result, local_res);
        }
    }
    return result;
}

namespace Case = Test::Case1;

Types::scalar neumann(const Types::point_t &x, const Types::face_t &face) {
    if (x.norm() > 2) {
        return 1;
    }
    return 0;
}

Types::scalar initial_conditions(const Types::point_t &x) { return 1; }

Types::scalar dirichlet(const Types::point_t &x) { return 1; }

Types::scalar rhs(const Types::point_t &x) { return 0; }

Types::Matrix3d D(const Types::point_t &x) { return (1./ 10) * Case::diff(x); }

Types::scalar c(const Types::point_t &x) { return 0; }

int main() {

    using method = Eigen::ConjugateGradient<Types::SparseMatrixXd, Eigen::Lower | Eigen::Upper>;
    using quadrature = Math::Integration::QuadratureTypes::GaussianPoints<3, 3>;
    using finite_elements = FEM::FiniteElements::P1;

    Eigen::setNbThreads(14);

    INMOST::Mesh mesh{};
    mesh.LoadMSH(
        "/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
        "examples/Lab2/meshes/2.msh");
    mesh.AssignGlobalID(INMOST::NODE);
    const auto debug_info = set_boundary_tag(mesh);
    std::cout << mesh.NumberOfCells() << std::endl;
    std::cout << mesh.NumberOfNodes() << std::endl;

    // решение методом конечных элементов
    const auto result = FEM::Time::Schemes::solveDynamicProblem<method, quadrature, finite_elements>(
        mesh, {0.001, 200}, initial_conditions, rhs, D, c, neumann, dirichlet);

    Mesh::Utils::assign_analytical_to_nodes(Case::analytical_solution, mesh, "analytical_stationary_solution");

    for (int i = 0; i < result.size(); i++) {
        Mesh::Utils::assign_to_nodes(result[i], mesh, "numerical_solution");
        mesh.Save("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
                  "examples/"
                  "results/fem/fem_mesh_" + std::to_string(i) + ".vtk");
    }
}
