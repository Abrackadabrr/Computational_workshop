#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <exception>

#include "inmost.h"

using namespace INMOST;

template<typename FuncTP>
double general_mean(Element& e, FuncTP func, double t){
    using FuncType = decltype(func);
    struct MyMeanFunc: public MeanFunc{
        FuncType f;
        MyMeanFunc(FuncType f): f{f} {}
        Storage::real func(Storage::real* x, Storage::real t) const override { return f(x, t); }
    };
    MyMeanFunc mmf(std::move(func));
    return e.Mean(mmf, t);
}

int main(int argc, char* argv[]){
    std::string mesh_file = "../../mesh/ventrical_new.msh";
    
    Mesh::Initialize(&argc, &argv);
    Solver::Initialize(&argc, &argv, nullptr);
    {
        Mesh m("ventricle");
        m.Load(mesh_file);
        Tag gmsh_lbl = m.GetTag("GMSH_TAGS");

        auto get_label_check = [lbl = gmsh_lbl](int val){
            return [lbl, val](Face f){
                return f.HaveData(lbl) && f.Integer(lbl) == val;
            };
        };
        auto is_base = get_label_check(1);
        auto is_endo = get_label_check(2);
        auto is_epi = get_label_check(3);

        double dx = 1, dy = 0.1, dz = 0.01;
        auto compute_D = [dx, dy, dz](const double* X/*[3]*/, rMatrix& D) -> void{
            double d = 1 + X[0]*X[0] + X[1]*X[1];
            D.Resize(3, 3);
            D.Zero();
            D(0, 0) = d*dx;
            D(1, 1) = d*dy;
            D(2, 2) = d*dz;
        };
        auto g_D = [](const double* X){
            double x = X[0], y = X[1], z = X[2];
            return exp((x + y + 20)/40) + (z + 10)*(z + 10)/20;
        };
        auto g_N = [](const double* X, Face f){
            double N[3];
            f.OrientedUnitNormal(f.BackCell(), N);
            double gN_val = 0;
            /*посчитайте -D * \nabla U * N*/
            return gN_val;
        };
        auto r = [](const double* X){
            double r_val = 0;
            /*посчитайте -\div D * \nabla U */
            return r_val;
        };
        auto is_neumann = [is_base](Face f){ return is_base(f); };
        auto is_dirichlet = [is_endo, is_epi](Face f){ return is_endo(f) || is_epi(f); };

        long NC = m.NumberOfCells();
        Sparse::Matrix A("A", 0, NC);
        Sparse::Vector b("b", 0, NC);
        
        rMatrix tmp;
        for (auto fi = m.BeginFace(); fi != m.EndFace(); ++fi){
            Face f = fi->getAsFace();
            auto lambda_cf = [compute_D, &Dc = tmp](Cell c, Face f){
                double xf[3], xc[3], N[3], Ns[3];
                f.Barycenter(xf);
                c.Barycenter(xc);
                f.OrientedUnitNormal(c, N);
                compute_D(xc, Dc);
                raMatrix Nm(shell<double>(N, 3), 1, 3),
                         Nsm(shell<double>(Ns, 3), 1, 3),
                         xfm(shell<double>(xf, 3), 1, 3),
                         xcm(shell<double>(xc, 3), 1, 3);
                Nsm = Dc * Nm;
                return Nsm.DotProduct(Nsm) / Nsm.DotProduct(xfm - xcm);
            };
            if (!f.Boundary()){
                Cell c1 = f.BackCell(), c2 = f.FrontCell();
                double lmb1 = lambda_cf(c1, f), lmb2 = lambda_cf(c2, f);
                double area = f.Area();
                double tau = area * lmb1 * lmb2 / (lmb1 + lmb2);
                auto i1 = c1.LocalID(), i2 = c2.LocalID();
                A[i1][i2] -= tau, A[i2][i1] -= tau;
                A[i1][i1] += tau, A[i2][i2] += tau;
            } else {
                Cell c = f.BackCell();
                auto i = c.LocalID();
                if (is_dirichlet(f)){
                    double xf[3];
                    f.Barycenter(xf);
                    double lmb = lambda_cf(c, f);
                    double area = f.Area();
                    A[i][i] += area * lmb;
                    b[i] += area * lmb * g_D(xf);
                } else {
                    b[i] -= general_mean(f, [g_N, f](double* X, double t){ return g_N(X, f); }, 0);
                }
            }
        }
        for (auto ci = m.BeginCell(); ci != m.EndCell(); ++ci){
            Cell c = ci->getAsCell();
            b[c.LocalID()] += general_mean(c, [r](double* X, double t){ return r(X); }, 0);
        }

        Solver solv(Solver::INNER_ILU2);
        solv.SetMatrix(A);
        Sparse::Vector U("U", 0, NC);
        bool is_solved = solv.Solve(b, U);

        Tag u_tag = m.CreateTag("U", DATA_REAL, CELL, NONE, 1);
        for (auto ci = m.BeginCell(); ci != m.EndCell(); ++ci)
            ci->Real(u_tag) = U[ci->LocalID()];

        m.Save("result.vtk");
    }
    Solver::Finalize();
    Mesh::Finalize();

    return 0;
}