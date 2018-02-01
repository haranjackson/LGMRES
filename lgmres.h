#ifndef LGMRES_H
#define LGMRES_H

#include "include/eigen3/Eigen"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Mat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef Eigen::Ref<Mat> Matr;
typedef Eigen::Ref<Vec> Vecr;
typedef Eigen::HouseholderQR<Mat> DecQR;
typedef std::function<Vec(Vecr)> VecFunc;


class System {
  Mat A;
  Mat M;

public:
  Vec call(Vecr x) { return A * x; }
  Vec prec(Vecr x) { return M * x; }
  void set_A(Mat A_) { A = A_; }
  void set_M(Mat M_) { M = M_; }
};

Vec lgmres(VecFunc matvec, VecFunc psolve, Vecr b, Vec x,
           std::vector<Vec> &outer_v, const double tol, const int maxiter,
           const int inner_m, const unsigned int outer_k);

Vec lgmres_wrapper(Matr A, Vecr b, Vecr x0, Matr M, const double tol,
                   const int maxiter, const int inner_m, const int outer_k,
                   std::vector<Vec> &outer_v);

#endif // LGMRES_H
