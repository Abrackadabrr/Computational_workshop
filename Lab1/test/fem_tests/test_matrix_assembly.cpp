//
// Created by evgen on 09.12.2024.
//
#include "fem/Parametrisation.hpp"
#include "fem/assembly/MatrixAssembly.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "integration/Node.hpp"
#include "integration/Quadrature.hpp"
#include "integration/Gaussian.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

#include <gtest/gtest.h>

using barycentric_v = FEM::FiniteElements::Barycentric::barycentric_v;

class ASSEMBLY : public ::testing::Test {
  protected:
    Types::mesh_t mesh;

    void SetUp() override {
        mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/Lab1/"
                     "examples/Lab2/meshes/2.msh");
        mesh.AssignGlobalID(INMOST::CELL);
        mesh.AssignGlobalID(INMOST::NODE);
    }

    Types::scalar c(const Types::point_t& x) {
        return 0;
    }
};

/*
 * Проверка корректности вычисления локальной матрицы
 */
TEST_F(ASSEMBLY, BASIC_TEST_LOCAL_MATRIX) {
    bool a = true;
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto id = [&](const Types::point_t &x) -> Types::Matrix3d { return Types::Matrix3d::Identity(); };

        const auto &res = FEM::Assembly::detail::getSubmatrix<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>,
                                                    FEM::FiniteElements::P1>(ielem->self(), id, c);
        // расчет руками
        const Types::Matrix3d &A_inverse = FEM::Parametrisation::getLocalParametrizationMatrix(ielem->self()).inverse();
        const Types::Matrix3d A_inverse_A_inv_tr = A_inverse * A_inverse.transpose();
        FEM::Assembly::detail::submatrix<4> expected_res{};

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                expected_res(i, j) =
                    (A_inverse_A_inv_tr * FEM::FiniteElements::P1::getGradient(i)(barycentric_v{0, 0, 0}))
                        .dot(FEM::FiniteElements::P1::getGradient(j)(barycentric_v{0, 0, 0}));
            }
        }

        expected_res = ielem->Volume() * expected_res;

        if (a) {
            std::cout << res << std::endl;
            std::cout << std::endl;
            std::cout << expected_res << std::endl;
            a = !a;
        }

        ASSERT_NEAR((res - expected_res).norm(), 0, 5e-15);
    }
}
#if 1
/**
 * Корректность вычисления локальной правой части
 */
TEST_F(ASSEMBLY, BASIC_TEST_LOCAL_RHS) {
    bool a = true;
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto id = [&](const Types::point_t &x) -> Types::scalar { return 1; };

        const auto &res =
            FEM::Assembly::detail::getSubrhs<Math::Integration::QuadratureTypes::GaussianPoints<3, 3>, FEM::FiniteElements::P1>(
                ielem->self(), id);
        const FEM::Assembly::detail::subvector<4> expected_res =
            ielem->Volume() * (1.0 / 4.0) * FEM::Assembly::detail::subvector<4>{1, 1, 1, 1};

        if (a) {
            std::cout << res << std::endl;
            std::cout << std::endl;
            std::cout << expected_res << std::endl;
            a = !a;
        }

        ASSERT_NEAR((res - expected_res).norm(), 0, 5e-15);
    }
}
#endif