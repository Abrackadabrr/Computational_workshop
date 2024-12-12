//
// Created by evgen on 08.12.2024.
//

#ifndef P1_HPP
#define P1_HPP

#include "types/BasicTypes.hpp"

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
};

struct P1 {
    static constexpr Types::index order = 1;
    static constexpr Types::index n = 4;
// функции
    static constexpr Types::scalar first(const Barycentric &point) { return point.l1; }

    static constexpr Types::scalar second(const Barycentric &point) { return point.l2; }

    static constexpr Types::scalar third(const Barycentric &point) { return point.l3; }

    static constexpr Types::scalar fourth(const Barycentric &point) { return point.l4; }
// градиенты
    static Types::Vector3d g_first(const Barycentric &point) { return {1.,0, 0}; }

    static Types::Vector3d g_second(const Barycentric &point) { return {0, 1, 0}; }

    static Types::Vector3d g_third(const Barycentric &point) { return {0, 0, 1}; }

    static Types::Vector3d g_fourth(const Barycentric &point) { return {-1, -1, -1}; }

    // thin place
    static std::function<Types::scalar(const Barycentric& point)> getFunction(Types::index i) {
        switch (i) {
            case 0:
                return first;
            case 1:
                return second;
            case 2:
                return third;
            case 3:
                return fourth;
            default:
                throw std::invalid_argument("Invalid index for getFunction in P1 finite_elements");
        }
    }

    static std::function<Types::Vector3d(const Barycentric& point)> getGradient(Types::index i) {
        switch (i) {
        case 0:
            return g_first;
        case 1:
            return g_second;
        case 2:
            return g_third;
        case 3:
            return g_fourth;
        default:
            throw std::invalid_argument("Invalid index for getGradient in P1 finite_elements");
        }
    }
};

struct P2 {
    static constexpr unsigned n = 10;
};

}

#endif

// P1_HPP
