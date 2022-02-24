#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <random>
#include <vector>

#include "wls.hpp"

#include <gtest/gtest.h>

static std::mt19937_64 eng(std::time(0));

static void fill_rand_array(double *v, const int n) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < n; ++i) v[i] = dist(eng);
}

static void fill_rand_array(std::complex<double> *v, const int n) {
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  for (int i = 0; i < n; ++i) {
    v[i].real(dist(eng));
    v[i].imag(dist(eng));
  }
}

TEST(BLAS1, MULTIPLY_ASSIGN) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n  = d(eng);
  double *                           v  = wls::alloc_mat<double>(n, 1);
  double *                           v0 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0, n);
  const double alpha = 1.253;
  wls::vec_scalar_multiply(v0, n, alpha, v);
  for (int i = 0; i < n; ++i) EXPECT_NEAR(v[i], alpha * v0[i], 1e-15);
  wls::free_mat(v);
  wls::free_mat(v0);
}

TEST(BLAS1, PROD) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n    = d(eng);
  double *                           v    = wls::alloc_mat<double>(n, 1);
  double *                           v0_1 = wls::alloc_mat<double>(n, 1),
         *v0_2                            = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0_1, n);
  fill_rand_array(v0_2, n);
  wls::vec_prod(v0_1, v0_2, n, v);
  for (int i = 0; i < n; ++i) EXPECT_NEAR(v[i], v0_1[i] * v0_2[i], 1e-15);
  wls::free_mat(v);
  wls::free_mat(v0_1);
  wls::free_mat(v0_2);
}

TEST(BLAS1, AXPY) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1);
  double *v0 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0, n);
  fill_rand_array(v1, n);
  std::copy_n(v1, n, v2);
  const double alpha = 1.253;
  wls::vec_axpy(v0, n, alpha, v1);
  for (int i = 0; i < n; ++i) EXPECT_NEAR(v1[i], alpha * v0[i] + v2[i], 1e-15);
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v0);
}

TEST(BLAS1, AXBYPZ) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1);
  double *v0_1 = wls::alloc_mat<double>(n, 1),
         *v0_2 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0_1, n);
  fill_rand_array(v0_2, n);
  fill_rand_array(v1, n);
  std::copy_n(v1, n, v2);
  const double alpha = 1.253, beta = -0.0312;
  wls::vec_axbypz(v0_1, v0_2, n, alpha, beta, v1);
  for (int i = 0; i < n; ++i)
    EXPECT_NEAR(v1[i], alpha * v0_1[i] + beta * v0_2[i] + v2[i], 1e-15);
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v0_1);
  wls::free_mat(v0_2);
}

TEST(BLAS1, SCALE) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n  = d(eng);
  double *                           v1 = wls::alloc_mat<double>(n, 1);
  double *                           v2 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v1, n);
  std::copy_n(v1, n, v2);
  const double alpha = 1.253;
  wls::vec_scale(n, alpha, v1);
  for (int i = 0; i < n; ++i) EXPECT_NEAR(v1[i], alpha * v2[i], 1e-15);
  wls::free_mat(v1);
  wls::free_mat(v2);
}

extern "C" double WLS_FC(ddot, DDOT)(lapack_int *, double *, lapack_int *,
                                     double *, lapack_int *);
TEST(BLAS1, DOT_REAL) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n  = d(eng);
  double *                           v1 = wls::alloc_mat<double>(n, 1);
  double *                           v2 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v1, n);
  fill_rand_array(v2, n);
  const auto v_dot  = wls::vec_dot(n, v1, v2);
  lapack_int nn     = n;
  lapack_int inc    = 1;
  const auto v_dot2 = WLS_FC(ddot, DDOT)(&nn, v1, &inc, v2, &inc);
  EXPECT_NEAR(v_dot, v_dot2, 1.e-10);
  wls::free_mat(v1);
  wls::free_mat(v2);
}

typedef struct {
  double real;
  double imag;
} blas_complex_t;
extern "C" blas_complex_t WLS_FC(zdotc, ZDOTC)(lapack_int *, blas_complex_t *,
                                               lapack_int *, blas_complex_t *,
                                               lapack_int *);
TEST(BLAS1, DOT_COMPLEX) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  auto *v1 = wls::alloc_mat<std::complex<double>>(n, 1);
  auto *v2 = wls::alloc_mat<std::complex<double>>(n, 1);
  fill_rand_array(v1, n);
  fill_rand_array(v2, n);
  const auto v_dot  = wls::vec_dot(n, v1, v2);
  lapack_int nn     = n;
  lapack_int inc    = 1;
  const auto v_dot2 = WLS_FC(zdotc, ZDOTC)(&nn, (blas_complex_t *)v1, &inc,
                                           (blas_complex_t *)v2, &inc);
  EXPECT_NEAR(v_dot.real(), v_dot2.real, 1e-13);
  EXPECT_NEAR(v_dot.imag(), v_dot2.imag, 1e-13);
  wls::free_mat(v1);
  wls::free_mat(v2);
}