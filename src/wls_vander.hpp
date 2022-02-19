/*
    This file is part of wlslib project

    Copyright (C) 2020 NumGeom Group at Stony Brook University

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*!
 * @file wls_vander.hpp
 *
 * @brief Construction of generalized and confluent Vandermonde matrices.
 * @todo Add support for confluent Vandermonde system with order > 0
 */

#ifndef WLS_VANDER_HPP_
#define WLS_VANDER_HPP_

#include <stdlib.h>

#include <algorithm>

#include "wls_config.hpp"

/*!
 * @def GEN_VANDER_UNIVAR0(__nv, __us, __degree, __V)
 *
 * @brief macro kernel for computing Vandermonde systems with single var
 *
 * @param[in] __nv number of nodes, integer type
 * @param[in] __us single independent variables, real pointer
 * @param[in] __degree degree of basis polynomial, integer type
 * @param[out] __V Vandermonde system output, real pointer
 */
#define GEN_VANDER_UNIVAR0(__nv, __us, __degree, __V)              \
  do {                                                             \
    std::fill_n(__V, __nv, T(1));                                  \
    const T *v_prev = __V;                                         \
    for (int i = 0; i < __degree; ++i) {                           \
      __V += __nv;                                                 \
      for (int j = 0; j < __nv; ++j) __V[j] = __us[j] * v_prev[j]; \
      v_prev = __V;                                                \
    }                                                              \
  } while (false)

/*!
 * @def GEN_VANDER_BIVAR0(__nv, __us, __degree, __V)
 *
 * @brief macro kernel for computing Vandermonde systems with two vars
 *
 * @param[in] __nv number of nodes, integer type
 * @param[in] __us single independent variables, real pointer, stride is 2
 * @param[in] __degree degree of basis polynomial, integer type
 * @param[out] __V Vandermonde system output, real pointer
 */
#define GEN_VANDER_BIVAR0(__nv, __us, __degree, __V)                          \
  do {                                                                        \
    std::fill_n(__V, __nv, T(1));                                             \
    const T *v_prev(__V);                                                     \
    for (int deg(0); deg < __degree; ++deg) {                                 \
      const int stop = deg + 1;                                               \
      for (int u(0); u < stop; ++u) {                                         \
        __V += __nv;                                                          \
        for (int j(0); j < __nv; ++j) __V[j] = __us[j << 1] * v_prev[j];      \
        v_prev += __nv;                                                       \
      }                                                                       \
      __V += __nv;                                                            \
      const T *vv_prev = v_prev - __nv;                                       \
      for (int j(0); j < __nv; ++j) __V[j] = __us[(j << 1) + 1] * vv_prev[j]; \
    }                                                                         \
  } while (false)

/*!
 * @def GEN_VANDER_TRIVAR0(__nv, __us, __degree, __V)
 *
 * @brief macro kernel for computing Vandermonde systems with three vars
 *
 * @param[in] __nv number of nodes, integer type
 * @param[in] __us single independent variables, real pointer, stride is 3
 * @param[in] __degree degree of basis polynomial, integer type
 * @param[out] __V Vandermonde system output, real pointer
 */
#define GEN_VANDER_TRIVAR0(__nv, __us, __degree, __V)                          \
  do {                                                                         \
    std::fill_n(__V, __nv, T(1));                                              \
    const T *v_prev(__V);                                                      \
    for (int deg(0); deg < __degree; ++deg) {                                  \
      const int stop1 = compute_ncols<2>(deg);                                 \
      for (int u(0); u < stop1; ++u) {                                         \
        __V += __nv;                                                           \
        for (int j(0); j < __nv; ++j) __V[j] = __us[(j << 2) - j] * v_prev[j]; \
        v_prev += __nv;                                                        \
      }                                                                        \
      const int stop2   = deg + 1;                                             \
      const T * vv_prev = v_prev - __nv * stop2;                               \
      for (int v(0); v < stop2; ++v) {                                         \
        __V += __nv;                                                           \
        for (int j(0); j < __nv; ++j)                                          \
          __V[j] = __us[(j << 2) - j + 1] * vv_prev[j];                        \
        vv_prev += __nv;                                                       \
      }                                                                        \
      __V += __nv;                                                             \
      vv_prev -= __nv;                                                         \
      for (int j(0); j < __nv; ++j)                                            \
        __V[j] = __us[(j << 2) - j + 2] * vv_prev[j];                          \
    }                                                                          \
  } while (false)

namespace wls {
/*!
 * @brief helper function to determine number of columns
 *
 * @tparam Dim coordinate dimension, must be in (1,2,3)
 * @param[in] degree degree of monomial degree
 * @sa compute_nrows
 */
template <int Dim>
inline constexpr int compute_ncols(const int degree) {
  static_assert(Dim > 0, "must be positive dimension");
  return compute_ncols<Dim - 1>(degree) * (degree + Dim) / Dim;
}

// specialize for dim == 1
template <>
inline constexpr int compute_ncols<1>(const int degree) {
  return degree + 1;
}

/*!
 * @brief helper routine to scale a matrix in row-wise
 *
 * @tparam T value type, either \a double or \a float
 * @param[in] m number of rows
 * @param[in] n number of columns
 * @param[in] w row weights, length of at least \a n
 * @param[in,out] A matrix to be scaled
 */
template <class T>
inline void scale_mat(const int m, const int n, const T *w, T *A) {
  // NOTE assume A is column-major and stride (aka leading dimension) is n
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < m; ++i) *A++ *= w[i];
}

/*!
 * @addtogroup interface
 * User interface
 * @{
 */

/*!
 * @brief Helper function to get the alignment size
 */
inline constexpr std::size_t get_alignment() { return WLS_ALIGN; }

/*!
 * @brief Helper function to determine aligned size
 *
 * @param[in] n input arbitrary size
 */
inline constexpr int determine_aligned_size(const int n) {
  return WLS_ALIGN + (n & (~(WLS_ALIGN - 1)));
}

/*!
 * @brief Allocate aligned 2D matrix data with *column-major*
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] m number of rows
 * @param[in] n number of columns
 *
 * This function will create data for a 2D column-major matrix with row and
 * column sizes of \a m and \a n, respectively. The memory guarantees to be
 * aligned for *each column* with aligned size \a WLS_ALIGN whose default value
 * is 32. Notice that the stride of the matrix will be given by
 * \ref determine_aligned_size
 *
 * @sa determine_aligned_size
 */
template <class T>
inline T *alloc_mat(const std::size_t m, const std::size_t n) {
  // require C11
  return reinterpret_cast<T *>(
      aligned_alloc(WLS_ALIGN, n * determine_aligned_size(m * sizeof(T))));
}

/*!
 * @brief Safely free memory
 *
 * @param[in] a pointer to be freed.
 */
inline void free_mat(void *a) {
  if (a) free(a);
}

/*!
 * @brief compute (weighted) Vandermonde system with single var
 *
 * @tparam Degree degree of basis polynomial
 * @tparam T value type, either \a double or \a float
 * @param[in] nv number of points
 * @param[in] us variable values
 * @param[out] V output Vandermonde system, stored in Fortran order
 * @param[in] w (optional) row weights, default is identity
 */
template <int Degree, class T>
inline void gen_vander_univar0(const int nv, const T *us, T *V,
                               const T *w = nullptr) {
  static_assert(Degree > -1, "degree must be nonnegative");
  GEN_VANDER_UNIVAR0(nv, us, Degree, V);
  if (w) scale_mat(nv, compute_ncols<1>(Degree), w, V);
}

/*!
* @brief Compute (weighted) Vandermonde systems with runtime degrees
*
* @param[in] nv number of points
* @param[in] us variable values
* @param[in] degree degree of basis polynomial
* @param[out] V output Vandermonde system, stored in Fortran order
* @param[in] w (optional) row weights, default is identity
* @note This functions calls \ref gen_vander_univar0<Degree,T> to optimize
///       commonly used degree choices.
*/
template <class T>
inline void gen_vander_univar0(const int nv, const T *us, const int degree,
                               T *V, const T *w = nullptr) {
  switch (degree) {
    case 2:
      gen_vander_univar0<2>(nv, us, V, w);
      break;
    case 3:
      gen_vander_univar0<3>(nv, us, V, w);
      break;
    case 4:
      gen_vander_univar0<4>(nv, us, V, w);
      break;
    default:
      GEN_VANDER_UNIVAR0(nv, us, degree, V);
      if (w) scale_mat(nv, compute_ncols<1>(degree), w, V);
      break;
  }
}

/*!
 * @brief compute (weighted) Vandermonde system with two independent vars
 *
 * @tparam Degree degree of basis polynomial
 * @tparam T value type, either \a double or \a float
 * @param[in] nv number of points
 * @param[in] us variable values, length is nv*2, stride is 2
 * @param[out] V output Vandermonde system, stored in Fortran order
 * @param[in] w (optional) row weights, default is identity
 */
template <int Degree, class T>
inline void gen_vander_bivar0(const int nv, const T *us, T *V,
                              const T *w = nullptr) {
  static_assert(Degree > -1, "degree must be nonnegative");
  GEN_VANDER_BIVAR0(nv, us, Degree, V);
  if (w) scale_mat(nv, compute_ncols<2>(Degree), w, V);
}

/*!
* @brief Compute (weighted) Vandermonde systems with runtime degrees
*
* @param[in] nv number of points
* @param[in] us variable values, length is nv*2, stride is 2
* @param[in] degree degree of basis polynomial
* @param[out] V output Vandermonde system, stored in Fortran order
* @param[in] w (optional) row weights, default is identity
* @note This functions calls \ref gen_vander_bivar0<Degree,T> to optimize
///       commonly used degree choices.
*/
template <class T>
inline void gen_vander_bivar0(const int nv, const T *us, const int degree, T *V,
                              const T *w = nullptr) {
  switch (degree) {
    case 2:
      gen_vander_bivar0<2>(nv, us, V, w);
      break;
    case 3:
      gen_vander_bivar0<3>(nv, us, V, w);
      break;
    case 4:
      gen_vander_bivar0<4>(nv, us, V, w);
      break;
    default:
      GEN_VANDER_BIVAR0(nv, us, degree, V);
      if (w) scale_mat(nv, compute_ncols<2>(degree), w, V);
      break;
  }
}

/*!
 * @brief compute (weighted) Vandermonde system with three independent vars
 *
 * @tparam Degree degree of basis polynomial
 * @tparam T value type, either \a double or \a float
 * @param[in] nv number of points
 * @param[in] us variable values, length is nv*3, stride is 3
 * @param[out] V output Vandermonde system, stored in Fortran order
 * @param[in] w (optional) row weights, default is identity
 */
template <int Degree, class T>
inline void gen_vander_trivar0(const int nv, const T *us, T *V,
                               const T *w = nullptr) {
  static_assert(Degree > -1, "degree must be nonnegative");
  GEN_VANDER_TRIVAR0(nv, us, Degree, V);
  if (w) scale_mat(nv, compute_ncols<3>(Degree), w, V);
}

/*!
* @brief Compute (weighted) Vandermonde systems with runtime degrees
*
* @param[in] nv number of points
* @param[in] us variable values, length is nv*3, stride is 3
* @param[in] degree degree of basis polynomial
* @param[out] V output Vandermonde system, stored in Fortran order
* @param[in] w (optional) row weights, default is identity
* @note This functions calls \ref gen_vander_trivar0<Degree,T> to optimize
///       commonly used degree choices.
*/
template <class T>
inline void gen_vander_trivar0(const int nv, const T *us, const int degree,
                               T *V, const T *w = nullptr) {
  switch (degree) {
    case 2:
      gen_vander_trivar0<2>(nv, us, V);
      break;
    case 3:
      gen_vander_trivar0<3>(nv, us, V);
      break;
    case 4:
      gen_vander_trivar0<4>(nv, us, V);
      break;
    default:
      GEN_VANDER_TRIVAR0(nv, us, degree, V);
      if (w) scale_mat(nv, compute_ncols<3>(degree), w, V);
      break;
  }
}

/*! @} */

}  // namespace wls

#endif  // WLS_VANDER_HPP_
