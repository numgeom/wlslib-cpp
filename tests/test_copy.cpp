#include <algorithm>
#include <cmath>
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

TEST(COPY, ONE2ONE) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n  = d(eng);
  double *                           v  = wls::alloc_mat<double>(n, 1);
  double *                           v0 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0, n);
  wls::vec_copy11(v0, n, v);
  for (int i = 0; i < n; ++i) EXPECT_EQ(v0[i], v[i]) << i;
  wls::free_mat(v);
  wls::free_mat(v0);
}

TEST(COPY, ONE2TWO) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1);
  double *v0 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0, n);
  wls::vec_copy12(v0, n, v1, v2);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0[i], v1[i]) << i;
    EXPECT_EQ(v0[i], v2[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v0);
}

TEST(COPY, ONE2THREE) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1),
         *v3 = wls::alloc_mat<double>(n, 1);
  double *v0 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0, n);
  wls::vec_copy13(v0, n, v1, v2, v3);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0[i], v1[i]) << i;
    EXPECT_EQ(v0[i], v2[i]) << i;
    EXPECT_EQ(v0[i], v3[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v3);
  wls::free_mat(v0);
}

TEST(COPY, TWO2TWO) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1);
  double *v0_1 = wls::alloc_mat<double>(n, 1),
         *v0_2 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0_1, n);
  fill_rand_array(v0_2, n);
  wls::vec_copy22(v0_1, v0_2, n, v1, v2);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0_1[i], v1[i]) << i;
    EXPECT_EQ(v0_2[i], v2[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v0_1);
  wls::free_mat(v0_2);
}

TEST(COPY, THREE2THREE) {
  std::uniform_int_distribution<int> d(1, 100);
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1),
         *v3   = wls::alloc_mat<double>(n, 1);
  double *v0_1 = wls::alloc_mat<double>(n, 1),
         *v0_2 = wls::alloc_mat<double>(n, 1),
         *v0_3 = wls::alloc_mat<double>(n, 1);
  fill_rand_array(v0_1, n);
  fill_rand_array(v0_2, n);
  fill_rand_array(v0_3, n);
  wls::vec_copy33(v0_1, v0_2, v0_3, n, v1, v2, v3);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0_1[i], v1[i]) << i;
    EXPECT_EQ(v0_2[i], v2[i]) << i;
    EXPECT_EQ(v0_3[i], v3[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v3);
  wls::free_mat(v0_1);
  wls::free_mat(v0_2);
  wls::free_mat(v0_3);
}
