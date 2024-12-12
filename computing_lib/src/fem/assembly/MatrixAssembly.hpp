//
// Created by evgen on 08.12.2024.
//

#ifndef MATRIXASSEMBLY_HPP
#define MATRIXASSEMBLY_HPP

#include "fem/Parametrisation.hpp"
#include "fem/assembly/Integration.hpp"
#include "fem/finite_elements/FiniteElements.hpp"
#include "fem/assembly/Utils.hpp"

#include "types/BasicTypes.hpp"
#include "types/MeshTypes.hpp"

#include "integration/Gaussian.hpp"
#include "integration/Quadrature.hpp"

namespace FEM::Assembly {
namespace detail {

template <unsigned int N> using submatrix = Eigen::Matrix<Types::scalar, N, N>;
template <unsigned int N> using subvector = Eigen::Vector<Types::scalar, N>;

template <typename BFT> Types::index getNumberOfDOF(const Types::mesh_t &mesh);

template <> inline Types::index getNumberOfDOF<FiniteElements::P1>(const Types::mesh_t &mesh) {
    return mesh.NumberOfNodes();
}

/**
 * Рассчитывает локальное выражение базисной функции на ячейке в декартовых координатах
 * @tparam BFT
 */
template <typename BFT>
decltype(auto) getBasisFunctionCartesianParametrisation(Types::index i, const Types::cell_t &cell) {
    const Types::Matrix3d &matrix = Parametrisation::getLocalParametrizationMatrix(cell).inverse();
    const Types::point_t &last_node = Mesh::Utils::getPoint(cell.getNodes()[3]);
    return [index = i, matrix, last_node](const Types::point_t &x) -> Types::scalar {
        const Types::Vector3d lambdas = matrix * (x - last_node);
        return BFT::getFunction(index)(
            FiniteElements::Barycentric::barycentric_v{lambdas.x(), lambdas.y(), lambdas.z()});
    };
}

/*
 * @tparam BFT for basis function type
 */
template <typename Quadrature, typename BFT, typename Callable_D, typename Callable_algebraic>
submatrix<BFT::n> getSubmatrix(const Types::cell_t &cell, // ячейка, по которой рассчитываем вспомогательную матрицу
                               const Callable_D &D, // тензор эллиптического оператора
                               const Callable_algebraic &c // коэффициент при алгебраическом члене
) {
    const Types::index n_Dof = BFT::n;
    submatrix<n_Dof> submatrix{};
    const Types::Matrix3d &A_inverse = Parametrisation::getLocalParametrizationMatrix(cell).inverse();

    // расчет интегралов от функций

    const auto tensor = [&A_inverse, &D](const Types::point_t &x) -> Types::Matrix3d {
        return A_inverse * D(x) * A_inverse.transpose();
    };

    for (Types::index i = 0; i < n_Dof; i++) {
        for (Types::index j = 0; j < n_Dof; j++) {

            const auto scalar_weight = [&](const FiniteElements::Barycentric &b) {
                return BFT::getFunction(i)(b) * BFT::getFunction(j)(b);
            };

            submatrix(i, j) = Integration::integrate_bilinear_form_over_cell<Quadrature>(
                                  cell, tensor, BFT::getGradient(i), BFT::getGradient(j)) +
                              Integration::integrate_over_cell_with_weight<Quadrature>(cell, c, scalar_weight);
        }
    }
    return submatrix;
}

/*
 * @tparam BFT тип базисных функций
 */
template <typename Quadrature, typename BFT, typename Callable, typename Callable_neumann>
subvector<BFT::n> getSubrhs(const Types::cell_t &cell, // ячейка для расчета вспомог. вектора
                            const INMOST::Tag &boundary_type, // тэг, заданный на face'ах
                            const Callable &rhs,              // функция правой части
                            const Callable_neumann &neumann // скалярные граничные условия неймана, зависят от
                                                            // точки x и от грани типа Types::face_t
) {
    subvector<BFT::n> subrhs{};
    for (int i = 0; i < BFT::n; i++) {
        subrhs[i] = Integration::integrate_over_cell_with_weight<Quadrature>(cell, rhs, BFT::getFunction(i));
    }

    // обработка граничного условия Неймана (и Робина)
    if (cell.Boundary()) {
        const auto &faces = cell.getFaces();
        for (auto i_face = faces.begin(), end = faces.end(); i_face != end; ++i_face) {
            const Types::face_t &face = i_face->getAsFace();
            if (face.Integer(boundary_type) == Mesh::BoundaryType::NEUMANN) {
                for (int i = 0; i < BFT::n; i++) {

                    const auto actual_neumann_with_weight = [&i, &cell, &face, &neumann](const Types::point_t &x) {
                        const auto basis_function = getBasisFunctionCartesianParametrisation<BFT>(i, cell);
                        return neumann(x, face) * basis_function(x);
                    };

                    subrhs[i] += Math::Integration::integrate_function_over_cell<
                        Math::Integration::QuadratureTypes::GaussianPoints<2, 3>>(actual_neumann_with_weight, face);
                }
            }
        }
    }

    return subrhs;
}
} // namespace detail

template <typename Quadrature, typename BFT, typename rsh_t, typename d_tensor_t, typename algebraic_t,
          typename neumann_t, typename dirichlet_t>
Types::SLAE getSLAE(Types::mesh_t &mesh, const rsh_t &rhs, const d_tensor_t &D, const algebraic_t &c,
                    const neumann_t &neumann, const dirichlet_t &dirichlet) {
    // результаты, которые будем заполнять
    Types::SparseMatrixXd matrix{detail::getNumberOfDOF<BFT>(mesh), detail::getNumberOfDOF<BFT>(mesh)};
    Types::VectorXd rhs_res = Types::VectorXd::Zero(detail::getNumberOfDOF<BFT>(mesh));

    // триплеты для задания разреженной матрицы
    Types::vector<Eigen::Triplet<Types::scalar>> triplets;
    triplets.reserve(mesh.NumberOfCells() * BFT::n * BFT::n);

    // тэг для характеристики граничного условия
    const INMOST::Tag &boundary_type = mesh.GetTag("Boundary_type");

    // основной цикл по ячейкам
    for (auto ielem = mesh.BeginCell(), end_cell = mesh.EndCell(); ielem != end_cell; ++ielem) {
        const Types::cell_t cell = ielem->self();
        const auto &dof = Utils::getLocalDOF<BFT>(cell);
        const auto &global_indexes = Utils::getGlobalIndexesOfDOF<BFT>(cell);
        auto A_sub = detail::getSubmatrix<Quadrature, BFT>(cell, D, c);

        const auto &A_sub_copy = A_sub;
        auto rhs_sub = detail::getSubrhs<Quadrature, BFT>(cell, boundary_type, rhs, neumann);

        // обработка граничного условия Дирихле
        // Нам нужно проитерироваться по всем "неправильным тестовым функциям"
        // Давайте проитерируемся по всем тестовым функциям (задана степенями свободы) и поймем какие из них "неправильные"
        for (Types::index dof_index = 0; dof_index < dof.size(); ++dof_index) {
            // Далее, мы хотим как-то оценить правильность этой степени свободы
            if (Utils::isIncorrect(dof[dof_index], boundary_type)) {
                // делаем модификацию
                const Types::index i0 = dof_index;
                // 1) Все строки с неправильными тестовыми функциями должны быть заменены в
                // нулевые с единицей на диагонали
                for (int k = 0; k < BFT::n; k++) {
                    A_sub(i0, k) = 0;
                }
                A_sub(i0, i0) = 1;
                // 2) и соотвествующие правые части должны быть равны значению на границе дирихле
                rhs_sub[i0] = dirichlet(Mesh::Utils::getPoint(dof[i0]));
                // Замена уравнений произведена, единтсвенное -- осталось несимметричной матрица
#if 0
                // Из оставшихся правильных уравнений нужно выкинуть все слагаемые с известными величинами
                for (Types::index k = 0; k < BFT::n; k++) {
                    // если я НЕ нашел индекс в массиве local_indexes, то тестовая функция была правильной
                    if (std::find(local_indexes.begin(), local_indexes.end(), k) == local_indexes.end()) {
                        for (Types::index p = 0; p < local_indexes.size(); p++) {
                            // переношу все в правую часть и зануляю соотвествующие слагаемые в матрице
                            rhs_sub[k] -= A_sub_copy(k, local_indexes[p]) * rhs_sub[local_indexes[p]];
                            A_sub(k, local_indexes[p]) = 0;
                        }
                    }
                }
#endif
                // Все
                // Уравнения неправильные выкинуты вообще
                // Уравнения правильные изменены с учетом известных данных
#if 0
                    std::cout << rhs_sub << std::endl << std::endl;
                    std::cout << A_sub << std::endl;
                    std::cout << "// ---------------- //" << std::endl;
#endif
            }
        }

        for (int i = 0; i < BFT::n; i++) {
            rhs_res[global_indexes[i]] += rhs_sub[i];
            for (int k = 0; k < BFT::n; k++) {
                triplets.push_back(Eigen::Triplet<Types::scalar>(global_indexes[i], global_indexes[k], A_sub(i, k)));
            }
        }
    }
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    return {matrix, rhs_res};
}

} // namespace FEM::Assembly

// умный кусок кода с модификацией матрицы и без копирования
#if 0
// std::cout << "DB" << std::endl;
                    // локальные индексы тестовых функций, которые не лежат в необходимом пространстве
                    const auto& local_indexes = Parametrisation::getLocalIndexesForDOFonFace<BFT>(face, cell);
                    for (int i = 0; i < local_indexes.size(); i++) {
                        // текущий индекс тестовой функции, не лежащей в нужном пространстве
                        // все вычисления с ней нужно выкинуть
                        Types::index i0 = local_indexes[i];
                        // для этого локального индекса у нас в принципе нет глобального уравнения
                        // (функция вне нужного пространства)
                        // поэтому с помощью это штуки мы хотим сыммитировать уравнение u_p = g_D(x_p)
                        // поэтому в правую часть мы поставим значение на границе дирихле,
                        // а в матрице просто сделаем нули всюду, кроме элемента A(i0,i0)
                        const Types::scalar dirichlet_value = dirichlet(Mesh::Utils::getPoint(nodes[i0]));
                        rhs_sub[i0] = dirichlet_value;
                        // 1) Модифицируем строку в матрице (чтобы уравнение в итоговой большой матрице
                        // имело вид N * u_i0_glob = N * g_D(x_i0_glob))
                        for (int k = 0; k < BFT::n; k++) {
                            if (k != i0)
                                A_sub(i0, k) = 0;
                        }
                        A_sub(i0, i0) = 1;

                        // Супер, мы расправились с неправильными уравнениями теми, которые были с
                        // тестовой функцией i0
                        // Но теперь у нас матрица получилась несимметрическая, что плохо
                        // Нам нужно выкинуть в правую часть значения, которые нам заранее известны из
                        // уравнений вида u_p = g_D(x_p)
                        // Слагаемые, которые предлагается выкинуть под интегралом содержат функции
                        // psi_p (p такой же как и в x_p)
                        // В нашем методе это буквально те же самые функции, что и не лежащие в необходимом
                        // пространстве, а поэтому их индексы нам известны -- это local_indexes.
                        // Давайте руками выкинем слагаемые, которые нам нужно убрать в правую часть
                        // В каждом уравнении для "правильных тестовых функций" с номером k
                        // я выбрасываю значение интеграла от неё с "неправильной базисной функцией",
                        // то есть с функцией, которая является одной из разложения для v, помноженным на значение
                        // степени свободы, для этой "неправильной базисной функции". При этом оно есть
                        // dirichlet_value, поскольку "неправильные базисные" и "неправильные тестовые" совпадают
                        // Если бы они различались, то нужно было бы сделать два цикла
                        // При этом уравнение с индексом i0 я пропускаю, потому что оно с "неправильной тестовой функцией"
                        // При этом очевидно я модифицирую и уравнения с ещё какими-то "неправильными тестовыми функциями",
                        // однако они потом точно так же изменятся при преобразовании
                        for (int k = 0; k < BFT::n; k++) {
                            // вообще по-хорошему тут k не должно равняться всем индексам из массива local_indexes
                            // но тут потом неправильные уравнения точно так же заменятся на правильные
                            // и так как в матрице нули будут в непр
                            if (k != i0)
                                rhs_sub[k] -= A_sub(k, i0) * dirichlet_value;
                        }
                        // Далее я зануляю соотвествующее слагаемое во вспомогательной матрице,
                        // тем самым убирая это слагаемое из расчета коэффициента в итоговой матрице
                        for (int k = 0; k < BFT::n; k++) {
                            if (k != i0)
                                A_sub(k, i0) = 0;
                        }
                    }
#endif

#endif //MATRIXASSEMBLY_HPP
