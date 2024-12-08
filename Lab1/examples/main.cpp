#include "function_traits/Extractions.hpp"
#include "integration/Gaussian.hpp"
#include "integration/Quadrature.hpp"
#include "mesh/MeshUtils.hpp"
#include "types/BasicTypes.hpp"

#include <inmost.h>
#include <iostream>

Types::scalar area(const Types::array_s<3> &x) { return 1; }

Types::scalar f1(const Types::array_s<3> &x) { return std::pow(x[0] + 1, 10); }

Types::scalar f2_1(const Types::array_s<3> &x) { return std::pow(x[0] + x[1] + 21, 1. / 3); }
Types::array_s<3> f2_2(const Types::array_s<3> &x) {
    return (1. / std::cosh(x[0] + x[1] + x[2]) + 1) * Types::array_s<3>{1, 1, 1};
}

Types::scalar f3(const Types::array_s<3> &x) { return std::log(x[0] + 11); }

template <typename Callable, int... indexes>
std::array<typename Extractions::function_info<Callable>::return_type, sizeof...(indexes)>
getSequence(const Callable &f, INMOST::Mesh &mesh, Mesh::CellType cellType,
            const std::integer_sequence<int, indexes...> & /**/) {
    return {Math::Integration::integrate<Math::Integration::QuadratureTypes::GaussianPoints<2, indexes>>(f, mesh,
                                                                                                         cellType)...};
}

int main() {
    INMOST::Mesh mesh{};
    mesh.LoadMSH("/home/evgen/Education/MasterDegree/5_level/Chapter_1/INM_RAS/"
                 "Computing_workshop/Workshop/Lab1/ventrical.msh");

    const auto res = getSequence(f2_2, mesh, Mesh::CellType::endocard, std::integer_sequence<int, 1, 2, 3, 4, 5>{});

    std::cout.precision(16);
    print(std::cout, res);
}
