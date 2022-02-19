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

TEST(FORWARDSOLVE_T, NO_STRIDE) {
  const auto      R_raw = gen_rand_array(100);
  const MatrixMap R((double *)R_raw.data(), 10, 10, 10);
  auto            bs  = gen_rand_array(20);
  auto            bs2 = bs;
  wls::rrqr_rtsolve(R_raw.data(), 10, 10, 10, 2, bs.data(), 10);
  for (int k = 0; k < 2; ++k) {
    auto *b = bs2.data() + k * 10;
    for (int i = 0; i < 10; ++i) {
      double s(0);
      for (int j = 0; j < i; ++j) s += R(j, i) * b[j];
      (b[i] -= s) /= R(i, i);
    }
  }
  for (int i = 0; i < 20; ++i)
    EXPECT_NEAR(bs[i], bs2[i], 1e-10 * std::abs(bs[i]));
}

TEST(FORWARDSOLVE_T, STRIDE) {
  // matrix's stride is 20
  const auto      R_raw = gen_rand_array(200);
  const MatrixMap R((double *)R_raw.data(), 10, 10, 20);
  // rhs stride is 15
  auto bs  = gen_rand_array(30);
  auto bs2 = bs;
  wls::rrqr_rtsolve(R_raw.data(), 10, 10, 20, 2, bs.data(), 15);
  for (int k = 0; k < 2; ++k) {
    auto *b = bs2.data() + k * 15;
    for (int i = 0; i < 10; ++i) {
      double s(0);
      for (int j = 0; j < i; ++j) s += R(j, i) * b[j];
      (b[i] -= s) /= R(i, i);
    }
    const auto *bb = bs.data() + k * 15;
    for (int i = 0; i < 10; ++i)
      EXPECT_NEAR(b[i], bb[i], 1e-10 * std::abs(b[i]));
  }
}
