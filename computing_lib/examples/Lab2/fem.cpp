//
// Created by evgen on 09.12.2024.
//

#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCholesky"

// fem
#include "fem/assembly/MatrixAssembly.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "integration/Gaussian.hpp"

// test
#include "TestFucntions.hpp"

// mesh
#include "mesh/MeshUtils.hpp"

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
                face.Integer(u_tag) = Mesh::BoundaryType::NEUMANN;
            else
                face.Integer(u_tag) = Mesh::BoundaryType::DIRICHLET;
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

Types::VectorXd solve(const Types::SLAE &slae) {
    Eigen::setNbThreads(14);
    Eigen::ConjugateGradient<Types::SparseMatrixXd, Eigen::Lower | Eigen::Upper> method;
    method.compute(slae.A);
    method.setTolerance(1e-10);
    const Types::VectorXd res = method.solve(slae.b);
    std::cout << method.info() << std::endl;
    std::cout << method.iterations() << std::endl;
    std::cout << method.error() << std::endl;
    return res;
}

void check_solution(Types::mesh_t &mesh, const std::vector<unsigned> &debug_info, const Types::SLAE &slae) {
    Types::index i = 0;
    for (auto ielem = mesh.BeginNode(); ielem != mesh.EndNode(); ++ielem) {
        if (ielem->Boundary()) {
            const unsigned j = ielem->GlobalID();
            if (std::find(debug_info.begin(), debug_info.end(), j) == debug_info.end()) {
                std::cout << "Numbering changes unexpectedly" << std::endl;
            }
#if 0
            std::cout << j << std::endl;
            std::cout << slae.A.row(j) << std::endl << std::endl;
            std::cout << slae.b[j] << std::endl << std::endl;
            std::cout << "// ------------------- //" << std::endl;
            ++i;
#endif
        }
    }
    std::cout << i << std::endl;
}

namespace Case = Test::Case1;

Types::scalar neumann(const Types::point_t &x, const Types::face_t &face) {
    return Mesh::Utils::Face::UnitNormal(face).dot(Case::DgradU(x));
}

Types::scalar dirichlet(const Types::point_t &x) { return Case::analytical_solution(x); }

Types::scalar rhs(const Types::point_t &x) { return Case::rhs(x); }

Types::Matrix3d D(const Types::point_t &x) { return Case::diff(x); }

Types::scalar c(const Types::point_t &x) { return 0; }

int main() {
    INMOST::Mesh mesh{};
    mesh.LoadMSH(
        "/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
        "examples/Lab2/meshes/1.msh");
    mesh.AssignGlobalID(INMOST::NODE);
    const auto debug_info = set_boundary_tag(mesh);
    std::cout << mesh.NumberOfCells() << std::endl;
    std::cout << mesh.NumberOfNodes() << std::endl;

    // решение методом конечных элементов
    const auto slae =
        FEM::Assembly::getSLAE<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>, FEM::FiniteElements::P1>(
            mesh, rhs, D, c, neumann, dirichlet);
    const Types::SparseMatrixXd &matrix = slae.A;
    std::cout << "Assembled" << std::endl;

    const auto res = solve(slae);

    Mesh::Utils::assign_to_nodes(res, mesh, "numerical_solution");

    Mesh::Utils::assign_analytical_to_nodes(Case::analytical_solution, mesh, "analytical");

    const auto analytical_vector = Mesh::Utils::evaluate_over_nodes(
        [](const Types::node_t &node) { return Case::analytical_solution(Mesh::Utils::getPoint(node)); }, mesh);

    Mesh::Utils::assign_to_nodes((res - analytical_vector).cwiseAbs(), mesh, "difference");

    mesh.Save("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/examples/"
              "results/fem_mesh_2.vtk");
}
