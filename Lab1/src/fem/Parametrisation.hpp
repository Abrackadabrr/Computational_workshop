//
// Created by evgen on 08.12.2024.
//

#ifndef PARAMETRISATION_HPP
#define PARAMETRISATION_HPP

#include "types/MeshTypes.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

namespace FEM {
/*
 * Надо убедиться обязательно, что глобальная нумерация узлов была проинициализирована
 */
Types::array<Types::integer, 4> inline getGlobalNumber(const Types::cell_t &cell) {
    const auto &nodes = cell.getNodes();
    // Тут появляется порядок узлов, а значит и порядок соотвествующих барицентрических координат:
    // номер координаты соотвествует номеру вершины.
    // Мы считаем, что конечные элементы заданы от первых трех координат, а четвертая координата
    // тривиально выражается через первые три: l4 = 1 - l1 - l2 - l3
    return {nodes[0].GlobalID(), nodes[1].GlobalID(), nodes[2].GlobalID(), nodes[3].GlobalID()};
}

Types::Matrix3d inline getLocalParametrizationMatrix(const Types::cell_t &cell) {
    Types::Matrix3d localMatrix = Types::Matrix3d::Zero();
    const INMOST::ElementArray<INMOST::Node> &nodes = cell.getNodes();
    const auto &last_node = Mesh::Utils::getPoint(nodes[3]);
    localMatrix.col(0) = (Mesh::Utils::getPoint(nodes[0]) - last_node);
    localMatrix.col(1) = (Mesh::Utils::getPoint(nodes[1]) - last_node);
    localMatrix.col(2) = (Mesh::Utils::getPoint(nodes[2]) - last_node);
    return localMatrix;
}

template<Types::index order>
decltype(auto) getDegreesOfFreedom(const Types::cell_t & cell) {
    return cell.getNodes();
}

}

#endif //PARAMETRISATION_HPP
