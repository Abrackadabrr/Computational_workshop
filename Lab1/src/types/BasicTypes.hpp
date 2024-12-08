//
// Created by evgen on 15.10.2024.
//

#ifndef BASICTYPES_HPP
#define BASICTYPES_HPP

#include <array>
#include <complex>

#include <Eigen/Sparse>

namespace Types {

using scalar = double;
using complex = std::complex<double>;
template <int N> using array_s = std::array<scalar, N>;
template <typename t, int N> using array = std::array<t, N>;

using SparseMatrixXd = Eigen::SparseMatrix<scalar>;
using VectorXd = Eigen::VectorXd;
using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;

struct SLAE {
    SparseMatrixXd A;
    VectorXd b;
};

template <typename T> using vector = std::vector<T>;
using point_t = Vector3d;

} // namespace Types

// some operators //
template <long unsigned int N>
inline Types::array_s<N> operator*(const Types::array_s<N> &lhs, const Types::array_s<N> &rhs) {
    Types::array_s<N> result;
    for (long unsigned int i = 0; i < N; i++) {
        result[i] = lhs[i] * rhs[i];
    }
    return result;
}

template <long unsigned int N> inline Types::array_s<N> operator*(Types::scalar lhs, const Types::array_s<N> &rhs) {
    Types::array_s<N> result;
    for (long unsigned int i = 0; i < N; i++) {
        result[i] = lhs * rhs[i];
    }
    return result;
}

template <long unsigned int N> inline Types::array_s<N> operator*(const Types::array_s<N> &lhs, Types::scalar rhs) {
    return lhs * rhs;
}

template <long unsigned int N>
inline Types::array_s<N> operator+(const Types::array_s<N> &lhs, const Types::array_s<N> &rhs) {
    Types::array_s<N> result;
    for (long unsigned int i = 0; i < N; i++) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

template <long unsigned N> void print(std::ostream &out, const Types::array_s<N> &array) {
    for (auto i : array)
        out << i << ',';
}

template <long unsigned N1> void print(std::ostream &out, const Types::array<Types::array_s<3>, N1> &array) {

    for (auto i : array) {
        out << '[';
        print(out, i);
        out << "]," << std::endl;
    }
}

#endif // BASICTYPES_HPP
