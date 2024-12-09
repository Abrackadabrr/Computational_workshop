//
// Created by evgen on 08.12.2024.
//

#include "fem/Parametrisation.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

#include <gtest/gtest.h>

using barycentric_v = FEM::FiniteElements::Barycentric::barycentric_v;

Types::scalar weight(const FEM::FiniteElements::Barycentric::barycentric_v & /**/) { return 1; }

Types::scalar f(const Types::point_t &x) { return 1; }

Types::point_t first_bar_coord(const Types::point_t &x, const Types::cell_t &cell) {
    const auto &nodes = cell.getNodes();
    const Types::Matrix3d &matrix = FEM::Parametrisation::getLocalParametrizationMatrix(cell).inverse();
    const Types::point_t lambda = matrix * (x - Mesh::Utils::getPoint(nodes[3]));
    return lambda;
}

class PARAMETRISATION : public ::testing::Test {
  protected:
    Types::mesh_t mesh;

    void SetUp() override {
        mesh.LoadMSH(
            "/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
            "examples/Lab2/meshes/2.msh");
        mesh.AssignGlobalID(INMOST::CELL);
        mesh.AssignGlobalID(INMOST::NODE);
    }
};

TEST_F(PARAMETRISATION, MATRIX_CONSTRUCTION) {
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto &matrix = FEM::Parametrisation::getLocalParametrizationMatrix(ielem->self());
        const auto &nodes = ielem->self()->getNodes();
        const Types::point_t n1 = Mesh::Utils::getPoint(nodes[0]);
        const Types::point_t n2 = Mesh::Utils::getPoint(nodes[1]);
        const Types::point_t n3 = Mesh::Utils::getPoint(nodes[2]);
        const Types::point_t n4 = Mesh::Utils::getPoint(nodes[3]);

        const Types::point_t n5 = 1. / 4 * (n1 + n2 + n3 + n4);

        const Types::point_t l1 = first_bar_coord(n1, ielem->self());
        const Types::point_t l2 = first_bar_coord(n2, ielem->self());
        const Types::point_t l3 = first_bar_coord(n3, ielem->self());
        const Types::point_t l4 = first_bar_coord(n4, ielem->self());
        const Types::point_t l5 = first_bar_coord(n5, ielem->self());

        ASSERT_NEAR(l1.x(), 1, 1e-15);
        ASSERT_NEAR(l1.norm(), 1, 1e-15);
        ASSERT_NEAR(l2.y(), 1, 1e-15);
        ASSERT_NEAR(l2.norm(), 1, 1e-15);
        ASSERT_NEAR(l3.z(), 1, 1e-15);
        ASSERT_NEAR(l3.norm(), 1, 1e-15);
        ASSERT_NEAR(l4.norm(), 0, 1e-15);
        std::cout << l5 << std::endl;
        ASSERT_NEAR(l5.x(), 0.25, 1e-14);
        ASSERT_NEAR(l5.y(), 0.25, 1e-14);
        ASSERT_NEAR(l5.z(), 0.25, 1e-14);
    }
}

TEST_F(PARAMETRISATION, FACE_NODES_INDEXES) {
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != mesh.EndCell(); ++ielem) {
        const auto &faces = ielem->getFaces();
        for (int i = 0; i < faces.size(); i++) {
            const auto &indexes = FEM::Parametrisation::getLocalIndexesForDOFonFace<FEM::FiniteElements::P1>(
                faces[i]->self(), ielem->self());

            for (int k = 0; k < indexes.size(); k++) {
                ASSERT_NEAR((Mesh::Utils::getPoint(faces[i].getNodes()[k]) -
                             Mesh::Utils::getPoint(ielem->getNodes()[indexes[k]]))
                                .norm(),
                            0, 0);
            }
        }
    }
}