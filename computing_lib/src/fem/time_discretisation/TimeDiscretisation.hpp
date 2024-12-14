//
// Created by evgen on 12.12.2024.
//

#ifndef TIMEDISCRETISATION_HPP
#define TIMEDISCRETISATION_HPP

#include "fem/Parametrisation.hpp"
#include "fem/assembly/Integration.hpp"
#include "fem/assembly/Utils.hpp"
#include "fem/assembly/MatrixAssembly.hpp"
#include "fem/finite_elements/FiniteElements.hpp"

#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"


namespace FEM::Time {
namespace detail {
template <typename Quadrature, typename BFT>
decltype(auto) getSubmatrixBeforeTimeDerivativeIsSemiDiscreteProblem(const Types::cell_t &cell) {
    return Assembly::detail::getSubmatrixAlgebraic<Quadrature, BFT>(cell,
                                                                    [](const Types::point_t & /**/) { return 1; });
}
} // namespace detail

/**
 *
 * @tparam Quadrature Квадратура, которой интегрируем по объему
 * @tparam BFT тип конечных элементов
 * @param mesh полная сетка
 * @return матрица при производной
 */
template <typename Quadrature, typename BFT>
Types::SparseMatrixXd getMatrix(Types::mesh_t &mesh) {
    // Мы хотим собрать матрицу, примерно такую же как и в случае с матрицей для пространственного оператора.
    // Имеется в виду, что мы хотим выкинуть уравнения, которые получены с помощью тестовых
    // функций, не принадлежащих нужному пространству в обобщенной постановке.
    // При этом, так как на границе Дирихле нет диф уравнения, то тут мы будем прямо занулять целую строку и столбец,
    // не заботясь о правой части, так как о ней мы уже позаботились в пространственном операторе.
    // Матрица получится вырожденной, но это не проблема для неявных смех расчета.

    Types::SparseMatrixXd matrix{FiniteElements::getNumberOfDOF<BFT>(mesh), FiniteElements::getNumberOfDOF<BFT>(mesh)};

    // триплеты для задания разреженной матрицы
    Types::vector<Eigen::Triplet<Types::scalar>> triplets;
    triplets.reserve(mesh.NumberOfCells() * BFT::n * BFT::n);

    // тэг для характеристики граничного условия
    const INMOST::Tag &boundary_type = mesh.GetTag("Boundary_type");

    // основной цикл по ячейкам
    for (auto ielem = mesh.BeginCell(), end_cell = mesh.EndCell(); ielem != end_cell; ++ielem) {
        const Types::cell_t cell = ielem->self();
        const auto &dof = Assembly::Utils::getLocalDOF<BFT>(cell);
        const auto &global_indexes = Assembly::Utils::getGlobalIndexesOfDOF<BFT>(cell);
        auto A_sub = detail::getSubmatrixBeforeTimeDerivativeIsSemiDiscreteProblem<Quadrature, BFT>(cell);

        // обработка граничного условия Дирихле
        // Нам нужно проитерироваться по всем "неправильным тестовым функциям"
        // Давайте проитерируемся по всем тестовым функциям (задана степенями свободы) и поймем какие из них
        // "неправильные"
        for (Types::index dof_index = 0; dof_index < dof.size(); ++dof_index) {
            // Далее, мы хотим как-то оценить правильность этой степени свободы
            if (Assembly::Utils::isIncorrect(dof[dof_index], boundary_type)) {
                // делаем модификацию
                const Types::index i0 = dof_index;
                // 1) Все строки с неправильными тестовыми функциями должны быть заменены в
                // нулевые с нулем на диагонали
                for (int k = 0; k < BFT::n; k++) {
                    A_sub(i0, k) = 0;
                }
                for (int k = 0; k < BFT::n; k++) {
                        A_sub(k, i0) = 0;
                }
                // 2) Тут ещё придется вычесть из правой части кусок при нестационарном условии дирихле,
                // но это не сегодня
            }
        }
        for (int i = 0; i < BFT::n; i++) {
            for (int k = 0; k < BFT::n; k++) {
                triplets.push_back(Eigen::Triplet<Types::scalar>(global_indexes[i], global_indexes[k], A_sub(i, k)));
            }
        }
    }
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return matrix;
}


} // namespace FEM::Time

#endif // TIMEDISCRETISATION_HPP
