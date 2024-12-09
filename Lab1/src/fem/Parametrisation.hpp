//
// Created by evgen on 08.12.2024.
//

#ifndef PARAMETRISATION_HPP
#define PARAMETRISATION_HPP

#include "types/MeshTypes.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

namespace FEM::Parametrisation{
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

/**
* Возвращает пронумерованные степени свободы на ячейке cell
*/
template<typename BFT>
decltype(auto) getLocalDOF(const Types::cell_t & cell) {
    return cell.getNodes();
}

/**
* Возвращает глобальные номера степеней свободы на ячейке cell
*
* Сейчас это работает только для элементов P1
*/
template<typename BFT>
Types::array<Types::index, BFT::n> getGlobalIndexesOfDOF(const Types::cell_t cell) {
    const auto& dof = getLocalDOF<BFT>(cell);
    assert(dof.size() == BFT::n);
    Types::array<Types::index, BFT::n> globalNumbers{};
    for (int i = 0; i < BFT::n; i++) {
        globalNumbers[i] = dof[i].GlobalID();
    }
    return globalNumbers;
}

}

#endif //PARAMETRISATION_HPP
