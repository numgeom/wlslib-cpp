#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#include "wls.hpp"

#include <gtest/gtest.h>

TEST(ASSIGN, ONE) {
  std::uniform_int_distribution<int> d(WLS_ALIGN / sizeof(double), 100);
  std::mt19937_64                    eng(std::time(0));
  const int                          n = d(eng);
  double *                           v = wls::alloc_mat<double>(n, 1);
  wls::vec_assign1(100.0, n, v);
  std::vector<double> v0(n);
  std::fill_n(v0.begin(), n, 100.0);
  for (int i = 0; i < n; ++i) EXPECT_EQ(v[i], v0[i]) << i;
  wls::free_mat(v);
}

TEST(ASSIGN, TWO) {
  std::uniform_int_distribution<int> d(1, 100);
  std::mt19937_64                    eng(std::time(0));
  const int                          n  = d(eng);
  double *                           v1 = wls::alloc_mat<double>(n, 1);
  double *                           v2 = wls::alloc_mat<double>(n, 1);
  wls::vec_assign2(100.0, n, v1, v2);
  std::vector<double> v0(n);
  std::fill_n(v0.begin(), n, 100.0);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0[i], v1[i]) << i;
    EXPECT_EQ(v0[i], v2[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
}

TEST(ASSIGN, THREE) {
  std::uniform_int_distribution<int> d(1, 100);
  std::mt19937_64                    eng(std::time(0));
  const int                          n = d(eng);
  double *v1 = wls::alloc_mat<double>(n, 1), *v2 = wls::alloc_mat<double>(n, 1),
         *v3 = wls::alloc_mat<double>(n, 1);
  wls::vec_assign3(100.0, n, v1, v2, v3);
  std::vector<double> v0(n);
  std::fill_n(v0.begin(), n, 100.0);
  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(v0[i], v1[i]) << i;
    EXPECT_EQ(v0[i], v2[i]) << i;
    EXPECT_EQ(v0[i], v3[i]) << i;
  }
  wls::free_mat(v1);
  wls::free_mat(v2);
  wls::free_mat(v3);
}