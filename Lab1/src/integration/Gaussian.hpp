//
// Created by evgen on 15.10.2024.
//

#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

#include "integration/Node.hpp"
#include "types/BasicTypes.hpp"

namespace Math::Integration::QuadratureTypes {
template <unsigned int dimention_, unsigned int order_> struct BaseGaussianPoints {
    static constexpr unsigned dim = dimention_;
    static constexpr unsigned order = order_;
};

template <unsigned int dimention_, unsigned int order_> struct GaussianPoints {};

// --- 1D --- //

template <> struct GaussianPoints<1, 1> : BaseGaussianPoints<1, 1> {
    using Base = BaseGaussianPoints<1, 1>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order> nodes = {Points::Node<Base::dim>{{1. / 2}, 1}};
};

template <> struct GaussianPoints<1, 2> : BaseGaussianPoints<1, 2> {
    using Base = BaseGaussianPoints<1, 2>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order> nodes = {Points::Node<Base::dim>{{0}, 1. / 2},
                                                                                 Points::Node<Base::dim>{{1}, 1. / 2}};
};

template <> struct GaussianPoints<1, 3> : BaseGaussianPoints<1, 3> {
    using Base = BaseGaussianPoints<1, 3>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order> nodes = {
        Points::Node<Base::dim>{{0.}, 1. / 6}, Points::Node<Base::dim>{{1. / 2}, 2. / 3},
        Points::Node<Base::dim>{{1.}, 1. / 6}};
};

template <> struct GaussianPoints<1, 5> : BaseGaussianPoints<1, 5> {
    using Base = BaseGaussianPoints<1, 5>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order - 2> nodes = {
        Points::Node<Base::dim>{{1. / 2 - (1. / 2) * std::sqrt(3. / 5)}, 5. / 18},
        Points::Node<Base::dim>{{1. / 2 + (1. / 2) * std::sqrt(3. / 5)}, 5. / 18},
        Points::Node<Base::dim>{{1. / 2}, 4. / 9},
    };
};

template <> struct GaussianPoints<1, 7> : BaseGaussianPoints<1, 7> {
    using Base = BaseGaussianPoints<1, 7>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order - 3> nodes = {
        Points::Node<Base::dim>{{1. / 2 - (1. / 2) * std::sqrt(3. / 7 - (2. / 7) * std::sqrt(6. / 5))},
                                1. / 4 + std::sqrt(30) / 72},
        Points::Node<Base::dim>{{1. / 2 + (1. / 2) * std::sqrt(3. / 7 - (2. / 7) * std::sqrt(6. / 5))},
                                1. / 4 + std::sqrt(30) / 72},
        Points::Node<Base::dim>{{1. / 2 - (1. / 2) * std::sqrt(3. / 7 + (2. / 7) * std::sqrt(6. / 5))},
                                1. / 4 - std::sqrt(30) / 72},
        Points::Node<Base::dim>{{1. / 2 + (1. / 2) * std::sqrt(3. / 7 + (2. / 7) * std::sqrt(6. / 5))},
                                1. / 4 - std::sqrt(30) / 72},
    };
};

template <> struct GaussianPoints<1, 9> : BaseGaussianPoints<1, 9> {
    using Base = BaseGaussianPoints<1, 9>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order - 4> nodes = {
        Points::Node<Base::dim>{{1. / 2}, 64. / 225},
        Points::Node<Base::dim>{{1. / 2 - (1. / 6) * std::sqrt(5 - 2. * std::sqrt(10. / 7))},
                                (322 + 13 * std::sqrt(70)) / 1800},
        Points::Node<Base::dim>{{1. / 2 + (1. / 6) * std::sqrt(5 - 2. * std::sqrt(10. / 7))},
                                (322 + 13 * std::sqrt(70)) / 1800},
        Points::Node<Base::dim>{{1. / 2 - (1. / 6) * std::sqrt(5 + 2. * std::sqrt(10. / 7))},
                                (322 - 13 * std::sqrt(70)) / 1800},
        Points::Node<Base::dim>{{1. / 2 + (1. / 6) * std::sqrt(5 + 2. * std::sqrt(10. / 7))},
                                (322 - 13 * std::sqrt(70)) / 1800},
    };
};

template <> struct GaussianPoints<1, 17> : BaseGaussianPoints<1, 17> {
    using Base = BaseGaussianPoints<1, 17>;
    static constexpr Types::scalar k2 = 0.968160239507626;
    static constexpr Types::scalar k3 = 0.836031107326636;
    static constexpr Types::scalar k4 = 0.613371432700590;
    static constexpr Types::scalar k5 = 0.324253423403809;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order - 8> nodes = {
        Points::Node<Base::dim>{{1. / 2}, 0.165119677500630},

        Points::Node<Base::dim>{{1. / 2 - k2 / 2}, 0.040637194180787},
        Points::Node<Base::dim>{{1. / 2 + k2 / 2}, 0.040637194180787},

        Points::Node<Base::dim>{{1. / 2 - k3 / 2}, 0.090324080347429},
        Points::Node<Base::dim>{{1. / 2 + k3 / 2}, 0.090324080347429},

        Points::Node<Base::dim>{{1. / 2 - k4 / 2}, 0.130305348201468},
        Points::Node<Base::dim>{{1. / 2 + k4 / 2}, 0.130305348201468},

        Points::Node<Base::dim>{{1. / 2 - k5 / 2}, 0.156173538520002},
        Points::Node<Base::dim>{{1. / 2 + k5 / 2}, 0.156173538520002},
    };
};

// --- 2D --- //

template <> struct GaussianPoints<2, 1> : BaseGaussianPoints<2, 1> {
    using Base = BaseGaussianPoints<2, 1>;
    static constexpr Types::array<Points::Node<Base::dim>, Base::order> nodes = {
        Points::Node<Base::dim>{{1. / 3, 1. / 3}, 1.},
    };
};

template <> struct GaussianPoints<2, 2> : BaseGaussianPoints<2, 2> {
    using Base = BaseGaussianPoints<2, 2>;
    static constexpr Types::array<Points::Node<Base::dim>, 3> nodes = {
        Points::Node<Base::dim>{{1. / 6, 1. / 6}, 1. / 3},
        Points::Node<Base::dim>{{1. / 6, 2. / 3}, 1. / 3},
        Points::Node<Base::dim>{{2. / 3, 1. / 6}, 1. / 3},
    };
};

template <> struct GaussianPoints<2, 3> : BaseGaussianPoints<2, 3> {
    using Base = BaseGaussianPoints<2, 3>;
    static constexpr Types::scalar l1 = 0.2319333685530305;
    static constexpr Types::scalar l2 = 0.6590276223740922;
    static constexpr Types::scalar l3 = 1 - l2 - l1;
    static constexpr Types::array<Points::Node<Base::dim>, 6> nodes = {
        Points::Node<Base::dim>{{l1, l2}, 1. / 6}, Points::Node<Base::dim>{{l2, l1}, 1. / 6},
        Points::Node<Base::dim>{{l1, l3}, 1. / 6}, Points::Node<Base::dim>{{l3, l1}, 1. / 6},
        Points::Node<Base::dim>{{l3, l2}, 1. / 6}, Points::Node<Base::dim>{{l2, l3}, 1. / 6},
    };
};

template <> struct GaussianPoints<2, 4> : BaseGaussianPoints<2, 4> {
    using Base = BaseGaussianPoints<2, 4>;
    static constexpr Types::scalar l1 = 0.44594849091596489;
    static constexpr Types::scalar w1 = 0.22338158967801147;
    static constexpr Types::scalar l2 = 0.091576213509770854;
    static constexpr Types::scalar w2 = 0.10995174365532188;
    static constexpr Types::array<Points::Node<Base::dim>, 6> nodes = {
        Points::Node<Base::dim>{{l1, l1}, w1},         Points::Node<Base::dim>{{l1, 1 - 2 * l1}, w1},
        Points::Node<Base::dim>{{1 - 2 * l1, l1}, w1},

        Points::Node<Base::dim>{{l2, l2}, w2},         Points::Node<Base::dim>{{l2, 1 - 2 * l2}, w2},
        Points::Node<Base::dim>{{1 - 2 * l2, l2}, w2},
    };
};

template <> struct GaussianPoints<2, 5> : BaseGaussianPoints<2, 5> {
    using Base = BaseGaussianPoints<2, 5>;
    static constexpr Types::scalar l1 = (6 + std::sqrt(15)) / 21;
    static constexpr Types::scalar w1 = (155 + std::sqrt(15)) / 1200;
    static constexpr Types::scalar w2 = (155 - std::sqrt(15)) / 1200;
    static constexpr Types::scalar l2 = (6 - std::sqrt(15)) / 21;
    static constexpr Types::array<Points::Node<Base::dim>, 7> nodes = {
        Points::Node<Base::dim>{{1. / 3, 1. / 3}, 0.225},

        Points::Node<Base::dim>{{l1, l1}, w1},
        Points::Node<Base::dim>{{l1, 1 - 2 * l1}, w1},
        Points::Node<Base::dim>{{1 - 2 * l1, l1}, w1},

        Points::Node<Base::dim>{{l2, l2}, w2},
        Points::Node<Base::dim>{{l2, 1 - 2 * l2}, w2},
        Points::Node<Base::dim>{{1 - 2 * l2, l2}, w2},
    };
};

// --- 3D --- //

template <> struct GaussianPoints<3, 1> : BaseGaussianPoints<3, 1> {
    using Base = BaseGaussianPoints<3, 1>;
    static constexpr Types::array<Points::Node<Base::dim>, 1> nodes = {Points::Node<Base::dim>{{0.25, 0.25, 0.25}, 1}};
};

template <> struct GaussianPoints<3, 2> : BaseGaussianPoints<3, 2> {
    using Base = BaseGaussianPoints<3, 2>;
    static constexpr Types::array<Points::Node<Base::dim>, 4> nodes = {
        Points::Node<Base::dim>{{0.1381966011250105, 0.1381966011250105, 0.1381966011250105}, 0.25},
        Points::Node<Base::dim>{{0.5854101966249685, 0.1381966011250105, 0.1381966011250105}, 0.25},
        Points::Node<Base::dim>{{0.1381966011250105, 0.5854101966249685, 0.1381966011250105}, 0.25},
        Points::Node<Base::dim>{{0.1381966011250105, 0.1381966011250105, 0.5854101966249685}, 0.25},
    };
};

template <> struct GaussianPoints<3, 3> : BaseGaussianPoints<3, 3> {
    using Base = BaseGaussianPoints<3, 3>;
    static constexpr Types::array<Points::Node<Base::dim>, 8> nodes = {
        Points::Node<Base::dim>{{0.3286811466653490, 0.3286811466653490, 0.3286811466653490}, 0.1274913115575064},
        Points::Node<Base::dim>{{0.013956560003953067, 0.3286811466653490, 0.3286811466653490}, 0.1274913115575064},
        Points::Node<Base::dim>{{0.3286811466653490, 0.013956560003953067, 0.3286811466653490}, 0.1274913115575064},
        Points::Node<Base::dim>{{0.3286811466653490, 0.3286811466653490, 0.013956560003953067}, 0.1274913115575064},

        Points::Node<Base::dim>{{0.1119207275092915, 0.1119207275092915, 0.1119207275092915}, 0.1225086884424935},
        Points::Node<Base::dim>{{0.6642378174721255, 0.1119207275092915, 0.1119207275092915}, 0.1225086884424935},
        Points::Node<Base::dim>{{0.1119207275092915, 0.6642378174721255, 0.1119207275092915}, 0.1225086884424935},
        Points::Node<Base::dim>{{0.1119207275092915, 0.1119207275092915, 0.6642378174721255}, 0.1225086884424935}};
};

template <> struct GaussianPoints<3, 5> : BaseGaussianPoints<3, 5> {
    using Base = BaseGaussianPoints<3, 5>;
    static constexpr Types::scalar l1 = 0.3108859192633006;
    static constexpr Types::scalar l2 = 0.0927352503108912;
    static constexpr Types::scalar l3 = 0.04550370412564965;
    static constexpr Types::array<Points::Node<Base::dim>, 14> nodes = {
        Points::Node<Base::dim>{{l1, l1, l1}, 0.1274913115575064},
        Points::Node<Base::dim>{{1 - 3 * l1, l1, l1}, 0.1274913115575064},
        Points::Node<Base::dim>{{l1, 1 - 3 * l1, l1}, 0.1274913115575064},
        Points::Node<Base::dim>{{l1, l1, 1 - 3 * l1}, 0.1274913115575064},

        Points::Node<Base::dim>{{l2, l2, l2}, 0.07349304311636194},
        Points::Node<Base::dim>{{1 - 3 * l2, l2, l2}, 0.07349304311636194},
        Points::Node<Base::dim>{{l2, 1 - 3 * l2, l2}, 0.07349304311636194},
        Points::Node<Base::dim>{{l2, l2, 1 - 3 * l2}, 0.07349304311636194},

        Points::Node<Base::dim>{{l3, l3, 1. / 2 - l3}, 0.04254602077708146},
        Points::Node<Base::dim>{{l3, 1. / 2 - l3, l3}, 0.04254602077708146},
        Points::Node<Base::dim>{{1. / 2 - l3, l3, l3}, 0.04254602077708146},
        Points::Node<Base::dim>{{l3, 1. / 2 - l3, 1. / 2 - l3}, 0.04254602077708146},
        Points::Node<Base::dim>{{1. / 2 - l3, l3, 1. / 2 - l3}, 0.04254602077708146},
        Points::Node<Base::dim>{{1. / 2 - l3, 1. / 2 - l3, l3}, 0.04254602077708146},
    };
};

} // namespace Math::Integration::QuadratureTypes

#endif // GAUSSIAN_HPP
