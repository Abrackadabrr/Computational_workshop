//
// Created by evgen on 08.12.2024.
//

#ifndef PARAMETRISATION_HPP
#define PARAMETRISATION_HPP

#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

namespace FEM::Parametrisation {
/**
 * Матрица Якоби для аффинного преобразования lambda -> x (референсный элемент -> cell)
 * @param cell ячейка, для которой ищем параметризацию
 */
Types::Matrix3d inline getLocalParametrizationMatrix(const Types::cell_t &cell) {
    Types::Matrix3d localMatrix = Types::Matrix3d::Zero();
    const INMOST::ElementArray<INMOST::Node> &nodes = cell.getNodes();
    const auto &last_node = Mesh::Utils::getPoint(nodes[3]);
    localMatrix.col(0) = (Mesh::Utils::getPoint(nodes[0]) - last_node);
    localMatrix.col(1) = (Mesh::Utils::getPoint(nodes[1]) - last_node);
    localMatrix.col(2) = (Mesh::Utils::getPoint(nodes[2]) - last_node);
    return localMatrix;
}

/*
 * Следующие три функции работают ТОЛЬКО для элементов вида P1, поскольку
 * для них степени свободы совпадают с вершинами элементов
 */

/**
 * Возвращает пронумерованные степени свободы на ячейке cell
 */
template <typename BFT> decltype(auto) getLocalDOF(const Types::cell_t &cell) { return cell.getNodes(); }


template <typename BFT>
decltype(auto) getLocalIndexesForDOFonFace(const Types::face_t &face, const Types::cell_t &cell) {
    const auto &local_dof = Mesh::Utils::getArray<Types::node_t>(cell.getNodes());
    const auto &face_nodes = face.getNodes();

    // полныйп перебор
    Types::array<Types::index, 3> result{};
    for (int i = 0; i < face_nodes.size(); i++) {
        const Types::index index =
            std::distance(local_dof.begin(), std::find(local_dof.begin(), local_dof.end(), face_nodes[i]));
        result[i] = index;
    }
    return result;
}

/**
 * Возвращает глобальные номера степеней свободы на ячейке cell
 *
 * Сейчас это работает только для элементов P1
 */
template <typename BFT> Types::array<Types::index, BFT::n> getGlobalIndexesOfDOF(const Types::cell_t cell) {
    const auto &dof = getLocalDOF<BFT>(cell);
    assert(dof.size() == BFT::n);
    Types::array<Types::index, BFT::n> globalNumbers{};
    for (int i = 0; i < BFT::n; i++) {
        globalNumbers[i] = dof[i].GlobalID();
    }
    return globalNumbers;
}

}

#endif //PARAMETRISATION_HPP
