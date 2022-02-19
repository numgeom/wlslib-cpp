#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include "wls.hpp"

#include <gtest/gtest.h>

class MatrixMap {
 public:
  MatrixMap(double *a, const int m, const int n) : _a(a), _m(m), _n(n) {}

  inline double operator()(const int i, const int j) const {
    return _a[i + j * _m];
  }

 private:
  double *_a;
  int     _m, _n;
};

static std::vector<double> gen_rand_array(const int n) {
  std::vector<double>                    v(n);
  std::mt19937_64                        eng(std::time(0));
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (auto &vv : v) vv = dist(eng);
  return v;
}

static std::vector<std::vector<double>> gen_vander_ref_1d(
    const std::vector<double> &us, const int degree) {
  std::vector<std::vector<double>> V(us.size());
  for (int i = 0; i < (int)us.size(); ++i) V[i].resize(degree + 1);
  for (int i = 0; i < (int)us.size(); ++i)
    for (int j = 0; j < degree + 1; ++j) V[i][j] = std::pow(us[i], j);
  return V;
}

static std::vector<std::vector<double>> gen_vander_ref_2d(
    const std::vector<double> &us, const int degree) {
  const int                        n = us.size() / 2;
  std::vector<std::vector<double>> V(n);
  for (int i = 0; i < n; ++i) V[i].resize(wls::compute_ncols<2>(degree));
  for (int i = 0; i < n; ++i) V[i][0] = 1.0;
  int cur = 1;
  for (int d = 1; d < degree + 1; ++d)
    for (int j = 0; j < d + 1; ++j) {
      for (int i = 0; i < n; ++i)
        V[i][cur] = std::pow(us[i << 1], d - j) * std::pow(us[(i << 1) + 1], j);
      ++cur;
    }
  return V;
}

static std::vector<std::vector<double>> gen_vander_ref_3d(
    const std::vector<double> &us, const int degree) {
  const int                        n = us.size() / 3;
  std::vector<std::vector<double>> V(n);
  for (int i = 0; i < n; ++i) V[i].resize(wls::compute_ncols<3>(degree));
  for (int i = 0; i < n; ++i) V[i][0] = 1.0;
  int cur = 1;
  for (int u = 1; u < degree + 1; ++u)
    for (int v = 0; v < u + 1; ++v)
      for (int w = 0; w < v + 1; ++w) {
        for (int i = 0; i < n; ++i)
          V[i][cur] = std::pow(us[3 * i], u - v) *
                      std::pow(us[3 * i + 1], v - w) *
                      std::pow(us[3 * i + 2], w);
        ++cur;
      }
  return V;
}

TEST(VANDER, ONE_D) {
  for (int degree = 1; degree <= 6; ++degree) {
    const int           n  = 10;
    const auto          us = gen_rand_array(n);
    std::vector<double> V(n * wls::compute_ncols<1>(degree));
    wls::gen_vander_univar0(n, us.data(), degree, V.data());
    const auto V_ref = gen_vander_ref_1d(us, degree);
    const auto VMap  = MatrixMap(V.data(), n, wls::compute_ncols<1>(degree));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < wls::compute_ncols<1>(degree); ++j)
        ASSERT_NEAR(VMap(i, j), V_ref[i][j], 1e-12);
  }
}

TEST(VANDER, TWO_D) {
  for (int degree = 1; degree <= 6; ++degree) {
    const int           n  = 10;
    const auto          us = gen_rand_array(2 * n);
    std::vector<double> V(n * wls::compute_ncols<2>(degree));
    wls::gen_vander_bivar0(n, us.data(), degree, V.data());
    const auto V_ref = gen_vander_ref_2d(us, degree);
    const auto VMap  = MatrixMap(V.data(), n, wls::compute_ncols<2>(degree));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < wls::compute_ncols<2>(degree); ++j)
        ASSERT_NEAR(VMap(i, j), V_ref[i][j], 1e-12);
  }
}

TEST(VANDER, THREE_D) {
  for (int degree = 1; degree <= 6; ++degree) {
    const int           n  = 10;
    const auto          us = gen_rand_array(3 * n);
    std::vector<double> V(n * wls::compute_ncols<3>(degree));
    wls::gen_vander_trivar0(n, us.data(), degree, V.data());
    const auto V_ref = gen_vander_ref_3d(us, degree);
    const auto VMap  = MatrixMap(V.data(), n, wls::compute_ncols<3>(degree));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < wls::compute_ncols<3>(degree); ++j)
        ASSERT_NEAR(VMap(i, j), V_ref[i][j], 1e-12);
  }
}
