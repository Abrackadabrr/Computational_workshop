//
// Created by evgen on 14.12.2024.
//

#ifndef SCHEMEPARAMETERS_HPP
#define SCHEMEPARAMETERS_HPP
#include <types/BasicTypes.hpp>

namespace Aux {

struct TimeIntegrationParameters {
    Types::scalar timeStep;
    Types::index N;
};

}

#endif //SCHEMEPARAMETERS_HPP
