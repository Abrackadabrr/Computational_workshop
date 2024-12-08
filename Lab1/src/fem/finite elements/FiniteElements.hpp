//
// Created by evgen on 08.12.2024.
//

#ifndef P1_HPP
#define P1_HPP

#include "types/BasicTypes.hpp"

#include <types/BasicTypes.hpp>

namespace FEM::FiniteElements {
/*
 * Базисы пространства Pr для референсного тетраэдрального конечного элемента (С, \Sigma, Pr)
 */

struct Barycentric {

    // Вектор барицентрических координат для трёхмерного пространства (без последней координаты)
    using barycentric_v = Types::array_s<3>;

    Types::scalar l1, l2, l3, l4;
    Barycentric(const barycentric_v &barycentric)
        : l1(barycentric[0]), l2(barycentric[1]), l3(barycentric[2]), l4(1 - l1 - l2 - l3) {}

#if 0
    Barycentric(const std::initializer_list<Types::scalar> &barycentric)
        : l1(barycentric[0]), l2(barycentric[1]), l3(barycentric[2]){};
#endif
};

struct P1 {
private:
    template<typename F, typename G>
    struct BasisFunction {
        F function;
        G gradient;
    };
public:
    static constexpr Types::index order = 1;
    static constexpr Types::index n = 4;

    static constexpr Types::scalar first(const Barycentric &point) { return point.l1; }

    static constexpr Types::scalar second(const Barycentric &point) { return point.l2; }

    static constexpr Types::scalar third(const Barycentric &point) { return point.l3; }

    static constexpr Types::scalar fourth(const Barycentric &point) { return point.l4; }
};

struct P2 {
    static constexpr unsigned n = 10;
};

}

#endif

// P1_HPP
