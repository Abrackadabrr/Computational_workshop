//
// Created by evgen on 08.12.2024.
//

#include "fem/Parametrisation.hpp"
#include "fem/assembly/Integration.hpp"
#include "fem/finite elements/FiniteElements.hpp"
#include "integration/Gaussian.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

#include <gtest/gtest.h>

using barycentric_v = FEM::FiniteElements::Barycentric::barycentric_v;

Types::scalar weight(const FEM::FiniteElements::Barycentric::barycentric_v & /**/) { return 1; }

Types::scalar f(const Types::point_t &x) { return 1; }

Types::scalar first_bar_coord(const Types::point_t &x, const Types::cell_t &cell) {
    const auto &nodes = cell.getNodes();
    const Types::Matrix3d &matrix = FEM::getLocalParametrizationMatrix(cell).inverse();
    const Types::point_t lambda = matrix * (x - Mesh::Utils::getPoint(nodes[3]));
    return lambda.x();
}

class INTEGRATION : public ::testing::Test {
  protected:
    Types::mesh_t mesh;

    void SetUp() override {
        mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/Lab1/"
                     "examples/Lab2/meshes/2.msh");
        mesh.AssignGlobalID(INMOST::CELL);
        mesh.AssignGlobalID(INMOST::NODE);
    }
};

/*
 * Проверка коннеткронсти базисных функций Р1
 */
TEST_F(INTEGRATION, BASIC_TEST_P1) {
    ASSERT_NEAR(FEM::FiniteElements::P1::first(FEM::FiniteElements::Barycentric::barycentric_v{0, 1, 0}), 0, 1e-16);
    ASSERT_NEAR(FEM::FiniteElements::P1::second(FEM::FiniteElements::Barycentric::barycentric_v{0, 1, 0}), 1, 1e-16);
    ASSERT_NEAR(FEM::FiniteElements::P1::third(FEM::FiniteElements::Barycentric::barycentric_v{0, 1, 0}), 0, 1e-16);
    ASSERT_NEAR(FEM::FiniteElements::P1::fourth(FEM::FiniteElements::Barycentric::barycentric_v{0, 1, 0}), 0, 1e-16);
    ASSERT_NEAR(FEM::FiniteElements::P1::first(FEM::FiniteElements::Barycentric::barycentric_v{0, 0, 0}), 0, 1e-16);
    ASSERT_NEAR(FEM::FiniteElements::P1::first(FEM::FiniteElements::Barycentric::barycentric_v{0, 0, 1}), 0, 1e-16);
}

TEST_F(INTEGRATION, INTEGRATION_BASIC_TEST_WITH_SIMPLE_F_W) {
#if 1
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const Types::scalar res =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), f, weight);
        ASSERT_NEAR(res, ielem->Volume(), 1e-15) << "Constant f constant weight integration failed" << std::endl;
    }

    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const Types::scalar res =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), f, FEM::FiniteElements::P1::first);
        ASSERT_NEAR(res, 0.25 * ielem->Volume(), 1e-15) << "Constant f Linear weight integration failed" << std::endl;
    }
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const Types::scalar res =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), f, FEM::FiniteElements::P1::second);
        ASSERT_NEAR(res, 0.25 * ielem->Volume(), 1e-15) << "Constant f Linear weight integration failed" << std::endl;
    }
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const Types::scalar res =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), f, FEM::FiniteElements::P1::third);
        ASSERT_NEAR(res, 0.25 * ielem->Volume(), 1e-15) << "Constant f Linear weight integration failed" << std::endl;
    }
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const Types::scalar res =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), f, FEM::FiniteElements::P1::fourth);
        ASSERT_NEAR(res, 0.25 * ielem->Volume(), 1e-15) << "Constant f Linear weight integration failed" << std::endl;
    }
#endif
#if 1

    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {

        const auto l1 = [&](const Types::point_t &x) { return first_bar_coord(x, ielem->self()); };

        const auto l4 = [&](const Types::point_t &x) {
            const auto &nodes = ielem->getNodes();
            const Types::Matrix3d &matrix = FEM::getLocalParametrizationMatrix(ielem->self()).inverse();
            const Types::point_t lambda = matrix * (x - Mesh::Utils::getPoint(nodes[3]));
            return 1 - lambda.x() - lambda.y() - lambda.z();
        };

        const Types::scalar res1 =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), l1, weight);
        const Types::scalar res4 =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), l4, weight);
        ASSERT_NEAR(res1, 0.25 * ielem->Volume(), 1e-14) << "Linear f constant weight integration failed" << std::endl;
        ASSERT_NEAR(res4, 0.25 * ielem->Volume(), 1e-14) << "Linear f constant weight integration failed" << std::endl;
    }
#endif
#if 1

    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {

        const auto l1 = [&](const Types::point_t &x) { return first_bar_coord(x, ielem->self()); };

        const auto l4 = [&](const Types::point_t &x) {
            const auto &nodes = ielem->getNodes();
            const Types::Matrix3d &matrix = FEM::getLocalParametrizationMatrix(ielem->self()).inverse();
            const Types::point_t lambda = matrix * (x - Mesh::Utils::getPoint(nodes[3]));
            return 1 - lambda.x() - lambda.y() - lambda.z();
        };

        const Types::scalar res1 =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), l1, FEM::FiniteElements::P1::first);
        const Types::scalar res4 =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), l4, FEM::FiniteElements::P1::fourth);
        ASSERT_NEAR(res1, 0.1 * ielem->Volume(), 1e-14) << "Linear f linear weight integration failed" << std::endl;
        ASSERT_NEAR(res4, 0.1 * ielem->Volume(), 1e-14) << "Linear f linear weight integration failed" << std::endl;
    }
#endif

#if 1
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {

        const auto l4_q = [&](const Types::point_t &x) {
            const auto &nodes = ielem->getNodes();
            const Types::Matrix3d &matrix = FEM::getLocalParametrizationMatrix(ielem->self()).inverse();
            const Types::point_t lambda = matrix * (x - Mesh::Utils::getPoint(nodes[3]));
            return (1 - lambda.x() - lambda.y() - lambda.z()) * (1 - lambda.x() - lambda.y() - lambda.z());
        };

        const Types::scalar res1 =
            FEM::Integration::integrate_over_cell_with_weight<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>>(
                ielem->self(), l4_q, FEM::FiniteElements::P1::first);

        ASSERT_NEAR(res1, (1./60) * ielem->Volume(), 1e-14) << "Quardatic f linear weight integration failed" << std::endl;
    }
#endif
}
