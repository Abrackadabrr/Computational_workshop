//
// Created by evgen on 12.12.2024.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

namespace FEM::Assembly::Utils {
/*
 * Следующие три функции работают ТОЛЬКО для элементов вида P1, поскольку
 * для них степени свободы совпадают с вершинами элементов
 */

/**
 * Возвращает пронумерованные степени свободы на ячейке cell
 */
template <typename BFT> decltype(auto) getLocalDOF(const Types::cell_t &cell) {
    return Mesh::Utils::getArray<Types::node_t>(cell.getNodes());
}

template <typename BFT>
decltype(auto) getLocalIndexesForDOFonFace(const Types::face_t &face, const Types::cell_t &cell) {
    const auto &local_dof = Mesh::Utils::getArray<Types::node_t>(getLocalDOF<BFT>(cell));
    // по-хорошему тут должен быть вызов функции, собирающей все степени свободы на грани
    const auto &face_nodes = Mesh::Utils::getArray<Types::node_t>(face.getNodes());

    // полный перебор
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
template <typename BFT> Types::array<Types::index, BFT::n> getGlobalIndexesOfDOF(const Types::cell_t &cell) {
    const auto &dof = getLocalDOF<BFT>(cell);
    assert(dof.size() == BFT::n);
    Types::array<Types::index, BFT::n> globalNumbers{};
    for (int i = 0; i < BFT::n; i++) {
        globalNumbers[i] = dof[i].GlobalID();
    }
    return globalNumbers;
}

/**
 * Определить, является ли тестовая функция неправильной
 */
bool inline isIncorrect(const Types::node_t& node, const INMOST::Tag& boundary_tag) {
    const auto faces = node.getFaces();
    for (Types::index i = 0; i < faces.size(); i++) {
        if (faces[i].Integer(boundary_tag) == Mesh::BoundaryType::DIRICHLET)
            return true;
    }
    return false;
}
}
#endif //UTILS_HPP
