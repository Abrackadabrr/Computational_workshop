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
            // if (face.Integer(gmsh_tag) == Mesh::cut)
            //     face.Integer(u_tag) = Mesh::BoundaryType::NEUMANN;
            // else
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
    method.setTolerance(1e-16);
    const Types::VectorXd res = method.solve(slae.b);
    std::cout << method.info() << std::endl;
    std::cout << method.iterations() << std::endl;
    std::cout << method.error() << std::endl;
    return res;
}

void check_solution(Types::mesh_t &mesh, const std::vector<unsigned>& debug_info,  const Types::SLAE &slae) {
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

Types::scalar neumann(const Types::point_t &x, const Types::face_t & /*no name*/) { return 1; }

Types::scalar dirichlet(const Types::point_t &x) { return 50; }

Types::scalar rhs(const Types::point_t &x) { return std::sin(x.norm()); }

Types::Matrix3d D(const Types::point_t &x) { return Types::Matrix3d::Identity(); }

Types::scalar c(const Types::point_t &x) { return 0; }

int main() {
    INMOST::Mesh mesh{};
    mesh.LoadMSH(
        "/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
        "examples/Lab2/meshes/2.msh");
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

    // std::cout << matrix.diagonal() << std::endl;

    // const Types::scalar det = slae.A.toDense().determinant();
    // std::cout << "Determinant: " << det << std::endl;
#if 1
    // проверка на самосопряженность
    for (int i = 0; i < slae.A.rows(); ++i) {
        for (int j = 0; j < slae.A.cols(); ++j) {
            if (std::abs(matrix.coeff(i, j) - matrix.coeff(j, i)) > 1e-15) {
                std::cout << "Not selfadjoint " << std::endl;
                std::cout << std::abs(matrix.coeff(i, j) - matrix.coeff(j, i)) << std::endl;
            }
        }
    }
    std::cout << (matrix * Types::VectorXd::Zero(slae.b.size()) - slae.b).norm() / slae.b.norm() << std::endl;
#endif
    const auto res = solve(slae);

    Mesh::Utils::assign_to_nodes(res, mesh, "numerical_solution");

    check_solution(mesh, debug_info, slae);

    mesh.Save("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/examples/"
              "results/fem_mesh_2.vtk");
}
