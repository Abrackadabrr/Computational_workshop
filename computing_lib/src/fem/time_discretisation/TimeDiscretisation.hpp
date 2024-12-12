//
// Created by evgen on 12.12.2024.
//

#ifndef TIMEDISCRETISATION_HPP
#define TIMEDISCRETISATION_HPP

#include "fem/assembly/MatrixAssembly.hpp"

namespace FEM::Time {
namespace detail {
template <typename Quadrature, typename BFT>
decltype(auto) getSubmatrixBeforeTimeDerivativeIsSemiDiscreteProblem(const Types::cell_t &cell) {
    return Assembly::detail::getSubmatrixAlgebraic<Quadrature, BFT>(cell,
                                                                    [](const Types::point_t & /**/) { return 1; });
}
}

template <typename Quadrature, typename BFT>
decltype(auto) getMatrix(const Types::mesh_t &mesh) {

}


} // namespace FEM::Time

#endif // TIMEDISCRETISATION_HPP
