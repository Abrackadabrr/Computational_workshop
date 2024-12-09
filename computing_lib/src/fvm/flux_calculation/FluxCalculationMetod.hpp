//
// Created by evgen on 08.11.2024.
//

#ifndef FLUXCALCULATIONMETOD_HPP
#define FLUXCALCULATIONMETOD_HPP

#include "types/BasicTypes.hpp"

namespace FVM::Flux {
enum class FluxCalculationMetod {
    TPFA = 0,
    OMPNT = 1,
};

/**
 * Хранение значения потока в специальном виде, пригодном для заполнения матрицы
 * @tparam N - размерность вектора коэффициентов и индексов
 */
template <int N>
struct flux {
    Types::array_s<N> coef;
    Types::array<int, N> indexes;
    Types::scalar b;
};

}

#endif // FLUXCALCULATIONMETOD_HPP
