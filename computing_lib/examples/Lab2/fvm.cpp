//
// Created by evgen on 08.11.2024.
//

#include "../examples/Lab2/TestFucntions.hpp"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCholesky"
#include "fvm/Scheme.hpp"
#include "mesh/MeshUtils.hpp"

void set_boundary_tag(INMOST::Mesh &mesh) {
    INMOST::Tag u_tag = mesh.CreateTag("Boundary_type", INMOST::DATA_INTEGER, INMOST::FACE, INMOST::FACE, 1);
    for (auto ci = mesh.BeginFace(); ci != mesh.EndFace(); ++ci) {
        auto face = ci->getAsFace();
        if (face.Boundary()) {
            face.Integer(u_tag) = Mesh::BoundaryType::NEUMANN;
        }
    }
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
    // std::cout << method.info() << std::endl;
    // std::cout << method.iterations() << std::endl;
    // std::cout << method.error() << std::endl;
    return res;
}

namespace testCase = Test::Case2;

int main() {
    INMOST::Mesh mesh{};
    mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
                 "examples/Lab2/meshes/1.msh");
    mesh.AssignGlobalID(INMOST::CELL);
    set_boundary_tag(mesh);
    std::cout << mesh.NumberOfCells() << std::endl;

    const auto slae =
        FVM::getSLAE(mesh, testCase::diff, testCase::rhs, testCase::analytical_solution, testCase::DgradU);

    std::cout << evaluate_mean_diameter(mesh) << ' ';

    const auto res = solve(slae);
    const Types::VectorXd error = res - FVM::getMeshProjection(mesh, testCase::analytical_solution);
    std::cout << Mesh::Utils::norm_over_mesh(error, mesh) / Mesh::Utils::norm_over_mesh(res, mesh) << std::endl;
    Mesh::Utils::fill(res, mesh, "numerical_solution");
    Mesh::Utils::fill(error, mesh, "numerical_error");
    Mesh::Utils::fill_scalar_with_function(testCase::analytical_solution, mesh, "U_anal");
    // mesh.Save("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/examples/"
    //           "results/mesh_3.vtk");
}
// 4.83 0.0300383
// 2.36 0.0272742

// 0.589887 0.0249499
// 0.700703 0.0262324
// 0.955483 0.0251233
// 1.15115 0.0245842
// 1.50582 0.0256542
// 1.88466 0.0262035
// 2.80281 0.0300383
