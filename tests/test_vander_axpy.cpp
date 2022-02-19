#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include "wls.hpp"

#include <gtest/gtest.h>

class MatrixMap {
 public:
  MatrixMap(double *a, const int m, const int n, const int stride)
      : _a(a), _m(m), _n(n), _stride(stride) {}

  inline double operator()(const int i, const int j) const {
    return _a[i + j * _stride];
  }

  inline double &operator()(const int i, const int j) {
    return _a[i + j * _stride];
  }

 private:
  double *_a;
  int     _m, _n, _stride;
};

static std::vector<double> gen_rand_array(const int n) {
  std::vector<double>                    v(n);
  std::mt19937_64                        eng(std::time(0));
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto &vv : v) vv = dist(eng);
  return v;
}

TEST(VANDER, NO_STRIDE) {
  const auto X     = gen_rand_array(200);
  auto       Y     = gen_rand_array(200);
  auto       Y_bak = Y;
  double     alpha = 2.0;
  wls::vander_axpy(alpha, X.data(), 20, 10, Y.data());
  for (int i = 0; i < 200; ++i) Y_bak[i] += alpha * X[i];
  for (int i = 0; i < 200; ++i) EXPECT_NEAR(Y[i], Y_bak[i], 1e-12);
}

TEST(VANDER, STRIDE) {
  const auto X      = gen_rand_array(200);
  auto       Y      = gen_rand_array(200);
  auto       Y_bak  = Y;
  double     alpha  = 2.0;
  const int  stride = 20;
  wls::vander_axpy(alpha, X.data(), 10, 10, Y.data(), stride);
  MatrixMap       Y_mat(Y_bak.data(), 10, 10, stride);
  const MatrixMap X_mat((double *)X.data(), 10, 10, stride),
      Y_axpy(Y.data(), 10, 10, stride);
  for (int j = 0; j < 10; ++j)
    for (int i = 0; i < 10; ++i) Y_mat(i, j) += alpha * X_mat(i, j);
  for (int j = 0; j < 10; ++j)
    for (int i = 0; i < 10; ++i) EXPECT_NEAR(Y_mat(i, j), Y_axpy(i, j), 1e-12);
}
