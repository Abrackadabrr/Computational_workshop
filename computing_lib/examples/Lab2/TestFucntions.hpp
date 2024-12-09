//
// Created by evgen on 22.11.2024.
//

#ifndef TESTFUCNTIONS_HPP
#define TESTFUCNTIONS_HPP

#include "types/BasicTypes.hpp"

namespace Test {
namespace Case2 {
inline Types::Matrix3d diff(const Types::point_t &point) {
    return Types::Vector3d{1, 1., 1.}.asDiagonal();
}

inline Types::scalar analytical_solution(Types::point_t point) {
    return std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) + (point.z() + 10) * (point.z() + 10) / 20;
}

inline Types::scalar rhs(const Types::point_t &point) {
    const Types::scalar exponent = std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) / 1600;
    return 2 * exponent + 1./10;
}

inline Types::Vector3d DgradU(const Types::point_t &point) {
    const Types::scalar exponent = std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) / 40;
    return -Types::Vector3d{exponent, exponent, (point.z() + 10) / 10};
}

}

namespace Case1 {
inline Types::Matrix3d diff(const Types::point_t &point) {
    return (1 + point.x() * point.x() + point.y() * point.y()) * Types::Vector3d{1, 1. / 2, 1. / 4}.asDiagonal();
}

inline Types::scalar analytical_solution(Types::point_t point) {
    return std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) + (point.z() + 10) * (point.z() + 10) / 20;
}

inline Types::scalar rhs(const Types::point_t &point) {
    const Types::scalar exponent = std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) / 40;
    const Types::scalar d_x = (1 + point.x() * point.x() + point.y() * point.y());
    return -((1. / 40) * ((3. / 2) * exponent + 1) * d_x + exponent * (2 * point.x() + point.y()));
}

inline Types::Vector3d DgradU(const Types::point_t &point) {
    const Types::scalar exponent = std::exp((point.x() + point.y()) / 40) * std::exp(1. / 2) / 40;
    const Types::scalar d_x = (1 + point.x() * point.x() + point.y() * point.y());
    return d_x * Types::Vector3d{exponent, 1. / 2 * exponent, (point.z() + 10) / 40};
}
}
} // namespace Test

#endif // TESTFUCNTIONS_HPP
