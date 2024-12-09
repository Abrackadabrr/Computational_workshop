//
// Created by evgen on 09.12.2024.
//

#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

#include <gtest/gtest.h>

class INMOST_UNDERSTANDING : public ::testing::Test {
protected:
    Types::mesh_t mesh;

    void SetUp() override {
        mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/Computing_workshop/Workshop/computing_lib/"
                     "examples/Lab2/meshes/2.msh");
        mesh.AssignGlobalID(INMOST::CELL);
        mesh.AssignGlobalID(INMOST::NODE);
    }

    Types::scalar c(const Types::point_t& x) {
        return 0;
    }
};
#if 0
TEST_F(INMOST_UNDERSTANDING, UNDERSTAND_LOCAL_ID) {
    for (auto ielem = mesh.BeginCell(), end = mesh.EndCell(); ielem != end; ++ielem) {
        const auto& nodes = ielem->getNodes();
        ASSERT_EQ(nodes[0].LocalID(), 0);
        ASSERT_EQ(nodes[1].LocalID(), 1);
        ASSERT_EQ(nodes[2].LocalID(), 2);
        ASSERT_EQ(nodes[3].LocalID(), 3);
    }
}
#endif