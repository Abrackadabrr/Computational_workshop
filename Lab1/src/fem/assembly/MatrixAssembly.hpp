//
// Created by evgen on 08.12.2024.
//

#ifndef MATRIXASSEMBLY_HPP
#define MATRIXASSEMBLY_HPP

#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"
#include "fem/Parametrisation.hpp"

namespace FEM {

namespace detail {

template<unsigned int N>
using submatrix = Eigen::Matrix<Types::scalar, N, N>;

/*
 * @tparam BFT for basis function type
 */
template<typename BFT>
submatrix<BFT::n> getSubmatrixOverCell(const Types::cell_t& cell) {
    submatrix<BFT::n> submatrix{};

    // возвращение степеней свободы, которые соотвествуют данному порядку базисных функций
    const auto& nodes = getDegreesOfFreedom<BFT::order>(cell);
    // Теперь у нас есть локальная нумерация узлов, а значит, и локальная нумерация базисных функций
    // Далее считаем вспомогательную матрицу, которая соответствует этой нумерации узлов

    return submatrix;
}

}

}

#endif //MATRIXASSEMBLY_HPP

