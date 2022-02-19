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
 * @file wls_lapack.hpp
 *
 * @brief LAPACK wrapper functions
 *
 */

#ifndef WLS_LAPACK_HPP_
#define WLS_LAPACK_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <type_traits>
#include <vector>

#include "wls_config.hpp"
#include "rrqr_trunc.h"

// wrapping LAPACK Fortran 77 interfaces
extern "C" {
// double precision

// qrcp from lapack
void WLS_FC(dgeqp3, DGEQP3)(lapack_int *, lapack_int *, double *, lapack_int *,
                            lapack_int *, double *, double *, lapack_int *,
                            lapack_int *);
// triangular matrix condition number
void WLS_FC(dtrcon, DTRCON)(char *, char *, char *, lapack_int *, double *,
                            lapack_int *, double *, double *, lapack_int *,
                            lapack_int *);
// implicitly apply op(Q)*rhs or rhs*op(Q)
void WLS_FC(dormqr, DORMQR)(char *, char *, lapack_int *, lapack_int *,
                            lapack_int *, double *, lapack_int *, double *,
                            double *, lapack_int *, double *, lapack_int *,
                            lapack_int *);
// triangular solve
void WLS_FC(dtrtrs, DTRTRS)(char *, char *, char *, lapack_int *, lapack_int *,
                            double *, lapack_int *, double *, lapack_int *,
                            lapack_int *);
// increment 2-norm condition number estimator
void WLS_FC(dlaic1, DLAIC1)(lapack_int *, lapack_int *, double *, double *,
                            double *, double *, double *, double *, double *);
// axpy for double precision
void WLS_FC(daxpy, DAXPY)(lapack_int *, double *, double *, lapack_int *,
                          double *, lapack_int *);

// single precision
// qrcp from lapack
void WLS_FC(sgeqp3, SGEQP3)(lapack_int *, lapack_int *, float *, lapack_int *,
                            lapack_int *, float *, float *, lapack_int *,
                            lapack_int *);
// triangular matrix condition number
void WLS_FC(strcon, STRCON)(char *, char *, char *, lapack_int *, float *,
                            lapack_int *, float *, float *, lapack_int *,
                            lapack_int *);
// implicitly apply op(Q)*rhs or rhs*op(Q)
void WLS_FC(sormqr, SORMQR)(char *, char *, lapack_int *, lapack_int *,
                            lapack_int *, float *, lapack_int *, float *,
                            float *, lapack_int *, float *, lapack_int *,
                            lapack_int *);
// triangular solve
void WLS_FC(strtrs, STRTRS)(char *, char *, char *, lapack_int *, lapack_int *,
                            float *, lapack_int *, float *, lapack_int *,
                            lapack_int *);
// increment 2-norm condition number estimator
void WLS_FC(slaic1, SLAIC1)(lapack_int *, lapack_int *, float *, float *,
                            float *, float *, float *, float *, float *);
// axpy for single precision
void WLS_FC(saxpy, SAXPY)(lapack_int *, float *, float *, lapack_int *, float *,
                          lapack_int *);
}

namespace wls {

// anonymous namespace for overloading LAPACK interface
namespace {
// qrcp
inline void geqp3(lapack_int *m, lapack_int *n, double *a, lapack_int *lda,
                  lapack_int *jpvt, double *tau, double *work,
                  lapack_int *lwork, lapack_int *info) {
  WLS_FC(dgeqp3, DGEQP3)(m, n, a, lda, jpvt, tau, work, lwork, info);
}
inline void geqp3(lapack_int *m, lapack_int *n, float *a, lapack_int *lda,
                  lapack_int *jpvt, float *tau, float *work, lapack_int *lwork,
                  lapack_int *info) {
  WLS_FC(sgeqp3, SGEQP3)(m, n, a, lda, jpvt, tau, work, lwork, info);
}

// tri cond#
inline void trcon(char *norm, char *uplo, char *diag, lapack_int *n, double *a,
                  lapack_int *lda, double *rcond, double *work,
                  lapack_int *lwork, lapack_int *info) {
  WLS_FC(dtrcon, DTRCON)
  (norm, uplo, diag, n, a, lda, rcond, work, lwork, info);
}
inline void trcon(char *norm, char *uplo, char *diag, lapack_int *n, float *a,
                  lapack_int *lda, float *rcond, float *work, lapack_int *lwork,
                  lapack_int *info) {
  WLS_FC(strcon, STRCON)
  (norm, uplo, diag, n, a, lda, rcond, work, lwork, info);
}

// apply Q
inline void ormqr(char *side, char *trans, lapack_int *m, lapack_int *n,
                  lapack_int *k, double *a, lapack_int *lda, double *tau,
                  double *c, lapack_int *ldc, double *work, lapack_int *lwork,
                  lapack_int *info) {
  WLS_FC(dormqr, DORMQR)
  (side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

inline void ormqr(char *side, char *trans, lapack_int *m, lapack_int *n,
                  lapack_int *k, float *a, lapack_int *lda, float *tau,
                  float *c, lapack_int *ldc, float *work, lapack_int *lwork,
                  lapack_int *info) {
  WLS_FC(sormqr, SORMQR)
  (side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info);
}

// tri solve
inline void trtrs(char *uplo, char *trans, char *diag, lapack_int *n,
                  lapack_int *nrhs, double *a, lapack_int *lda, double *b,
                  lapack_int *ldb, lapack_int *info) {
  WLS_FC(dtrtrs, DTRTRS)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
}

inline void trtrs(char *uplo, char *trans, char *diag, lapack_int *n,
                  lapack_int *nrhs, float *a, lapack_int *lda, float *b,
                  lapack_int *ldb, lapack_int *info) {
  WLS_FC(strtrs, STRTRS)(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
}

// 2-norm incremental cond. # estimator
inline void laic1(lapack_int *job, lapack_int *j, double *x, double *sest,
                  double *w, double *gamma, double *sestpr, double *s,
                  double *c) {
  WLS_FC(dlaic1, DLAIC1)
  (job, j, x, sest, w, gamma, sestpr, s, c);
}

inline void laic1(lapack_int *job, lapack_int *j, float *x, float *sest,
                  float *w, float *gamma, float *sestpr, float *s, float *c) {
  WLS_FC(slaic1, SLAIC1)
  (job, j, x, sest, w, gamma, sestpr, s, c);
}

// axpy
inline void axpy(lapack_int *n, double *alpha, double *x, lapack_int *incx,
                 double *y, lapack_int *incy) {
  WLS_FC(daxpy, DAXPY)(n, alpha, x, incx, y, incy);
}

inline void axpy(lapack_int *n, float *alpha, float *x, lapack_int *incx,
                 float *y, lapack_int *incy) {
  WLS_FC(saxpy, SAXPY)(n, alpha, x, incx, y, incy);
}
}  // namespace

/*!
 * @addtogroup interface
 * @{
 */

/*!
 * @brief query for the work space that is needed for best performace in QRCP
 *
 * @param[in] m number of rows
 * @param[in] n number of columns
 * @return optimal buffer size
 */
inline int query_work_size(const int m, const int n) {
  const int NB        = 32;
  int       work_size = std::max(m, n + 1) * NB + std::max(4 * n, 64 * 65);

#ifdef DEBUG
  lapack_int query = -1, one = 1;
  char       left('L'), notran('N'), tran('T');

  lapack_int mm(m), nn(n), info;
  double     size[3];
  // query size for qrcp
  geqp3(&mm, &nn, nullptr, &mm, nullptr, nullptr, size, &query, &info);
  // We allocate extra space for permutation vector with proper alignment
  size[0] += static_cast<int>(size[0]) % 1 + 2 * n;

  // query size for Q*x
  ormqr(&left, &notran, &mm, &one, &nn, nullptr, &mm, nullptr, nullptr, &mm,
        size + 1, &query, &info);
  size[1] += static_cast<int>(size[1]) % 1;
  // query size for Q'*x
  ormqr(&left, &tran, &mm, &one, &nn, nullptr, &mm, nullptr, nullptr, &mm,
        size + 2, &query, &info);
  size[2] += static_cast<int>(size[2]) % 1;

  int msz1 = static_cast<int>(std::max(std::max(size[0], size[1]), size[2]));
  // NOTE: 3*n is used in triangular condition num. estimator
  assert(work_size >= msz1);
#endif

  return work_size;
}

/*!
 * @brief determine numerical rank of an upper triangular matrix
 *
 * @tparam T Value type, either \a double or \a float
 * @param[in] n Dimension of the triangular matrix
 * @param[in] R Triangular matrix
 * @param[in] stride Column stride of \a R
 * @param[in] cond_thres Threshold of 2-norm condition number
 * @param[in] work Work space of size at least 6*n
 * @return Numerical rank if nonnegative, negative value for LAPACK error.
 * @note The parameter \a iwork must be type of \a lapack_int
 *
 */
template <class T>
inline int determine_rank(const int n, const T *R, const int stride,
                          const T cond_thres, T *work) {
  char       norm('1'), uplo('U'), diag('N');
  lapack_int job_min(2), job_max(1);

  // NOTE: work has size of at least 6*n in type T. We use first 3*n for
  // floating-point work space and start int work space at 4*n for alignment
  lapack_int *iwork = reinterpret_cast<lapack_int *>(&work[4 * n]);

  // first, determine whether the matrix is numerically rank deficient
  lapack_int n_(n), ldr(stride), info(0);
  T          rcond;
  trcon(&norm, &uplo, &diag, &n_, (T *)R, &ldr, &rcond, work, iwork, &info);

  if (info) return info;

  lapack_int rank(n);
  if (rcond * cond_thres < T(10.)) {
    // if 1-norm based condition number if larger than 0.1 times the threshold,
    // then we incrementally estimate the 2-norm condition numbers

    T *x = work, *y = work + n;
    x[0] = y[0] = T(1);
    T smax(std::abs(*R)), smin(smax);  // initial singular values
    T s1, c1, s2, c2, sminpr, smaxpr;
    for (rank = 0; rank < n; ++rank) {
      // current column, and w+rank is the diagonal entry
      T *w = const_cast<T *>(R + ldr * rank);
      // estimate smaller
      laic1(&job_min, &rank, x, &smin, w, w + rank, &sminpr, &s1, &c1);
      // estimate larger
      laic1(&job_max, &rank, y, &smax, w, w + rank, &smaxpr, &s2, &c2);
      if (smaxpr <= sminpr * cond_thres) {
        // still well-conditioned
        for (int i = 0; i < rank; ++i) {
          x[i] *= s1;
          y[i] *= s2;
        }
        x[rank] = c1;
        y[rank] = c2;
        smin    = sminpr;
        smax    = smaxpr;
        continue;
      }
      break;
    }
  }
  return rank;
}

/*!
 * @brief core function for rank-revealing QR with column pivoting
 *
 * @tparam T Value type, either \a double or \a float
 * @param[in] A   Weighted generalized Vandermode matrix
 * @param[in] cond_thres Threshold of 2-norm condition number
 * @param[in] m  Number of rows
 * @param[in] n  Number of columns
 * @param[out] B Upon output, this is the internal structure for QRCP with tau
 * @param[in,out] jpvt Column pivoting array (1-based index), size of n;
 *                    if jpvt[0]>0 at input, it will be used to permute A.
 * @param[out] work Internal work space
 * @param[in]  work_size Size of internal floating-point work space,
 *                       should equal to query_work_size
 * @param[in] stride (optional) Stride for A and B
 * @return rank If nonnegative, numerical rank of the matrix.
 *         If negative, an LAPACK error.
 *
 * @sa query_work_size, rrqr_factor
 *
 */
template <class T>
inline int rrqr_factor_nodag(const T *A, const T cond_thres, const int m,
                             const int n, T *B, int *jpvt, T *work,
                             const int work_size, const int stride = 0) {
  lapack_int lda(stride <= 0 ? m : stride);
  T *        tau = &B[lda * n];

  // first step, copy A to B
  // NOTE: the code assumes same stride for A and B
  if (jpvt[0] == 0) {
    if (m != lda)
      for (int j = 0; j < n; ++j) std::copy_n(&A[lda * j], m, &B[lda * j]);
    else
      std::copy_n(A, lda * n, B);
  } else {
    // we have a prior permutation
    for (int j = 0; j < n; ++j)
      std::copy_n(&A[lda * (jpvt[j] - 1)], m, &B[lda * j]);
  }

  // second step, init column pivoting
  lapack_int *jpvt_(nullptr);

  if (sizeof(int) == sizeof(lapack_int)) {
    jpvt_ = reinterpret_cast<lapack_int *>(jpvt);
  } else {
    // NOTE: this block of code should be optimized out if integer type is
    // consistent during compilation.
    // use last 2n entries as buffer space for jpvt_
    jpvt_ = reinterpret_cast<lapack_int *>(&work[work_size - 2 * n]);
  }

  jpvt_[0] = 1;  // fix leading column to enforce partition of unity
  std::fill_n(jpvt_ + 1, n - 1, 0);  // mark other columns to be free

  // third step, call QRCP
  lapack_int m_(m), n_(n), work_size_(work_size - 2 * n), info(0);
  geqp3(&m_, &n_, B, &lda, jpvt_, tau, work, &work_size_, &info);
  if (info) return info;

  // copy column pivoting to the user array if necessary
  if (sizeof(int) != sizeof(lapack_int)) std::copy_n(jpvt_, n, jpvt);

  // fourth step, determine whether the matrix is rank deficient
  return determine_rank(std::min(m, n), B, lda, cond_thres, work);
}

/*!
 * @brief core function for rank-revealing QR with column pivoting (with trunc)
 *
 * @tparam T Value type, either \a double or \a float
 * @param[in] A   Weighted generalized Vandermode matrix
 * @param[in] cond_thres Threshold of 2-norm condition number
 * @param[in] m  Number of rows
 * @param[in] n  Number of columns
 * @param[out] B Upon output, this is the internal structure for QRCP with tau
 * @param[in,out] jpvt Column pivoting array (1-based index), size of n;
 *                    if jpvt[0]>0 at input, it will be used to permute A.
 * @param[out] work Internal work space
 * @param[in]  work_size Size of internal floating-point work space,
 *                       should equal to query_work_size
 * @param[in] stride (optional) Stride for A and B
 * @param[in] dag (optional) DAG for columns in Vandermonde system
 * @param[in] dim (optional) Dimension of DAG (1:3), if \a dag is given, then
 *                \a dim must be provided.
 * @param[in] nlvls_dag (optional) Number of nodes/levels in the DAG, should be
 *                      no smaller than \a n.
 * @return rank If nonnegative, numerical rank of the matrix.
 *         If negative, an LAPACK error.
 *
 * @sa query_work_size, rrqr_factor_nodag
 *
 */
template <class T>
inline int rrqr_factor(const T *A, const T cond_thres, const int m, const int n,
                       T *B, int *jpvt, T *work, const int work_size,
                       const int stride = 0, const unsigned char *dag = nullptr,
                       const int dim = 0, const int nlvls_dag = 0) {
  if (n) *jpvt = 0;
  int rank =
      rrqr_factor_nodag(A, cond_thres, m, n, B, jpvt, work, work_size, stride);
  if (rank < 0) return rank;
  if (rank == n || !dag) return rank;
  const int nlvls_dag_(nlvls_dag < n ? n : nlvls_dag);
  if (dim < 1 || dim > 3) return -111;  // error inputs
  // create coder dag interface
  coder::array<unsigned char, 2U> dag_array;
  dag_array.set(const_cast<unsigned char *>(dag), dim, nlvls_dag_);
  coder::array<int, 1U> jpvt_array;
  jpvt_array.set(jpvt, n);
  coder::array<int, 2U> work_array;
  work_array.set(reinterpret_cast<int *>(work), 4, n);
#ifdef DEBUG
  if (dag_array.is_owner()) return -222;  // should not be owner
#endif
  int n1 = n;
  for (;;) {
    if (rrqr_trunc(dag_array, &n1, rank, jpvt_array, work_array)) {
      // if we have permutation, continue factorizing QRCP with updated columns
      rank = rrqr_factor_nodag(A, cond_thres, m, n1, B, jpvt, work, work_size,
                               stride);
      continue;
    }
    break;
  }
  return rank;
}

/*!
 * @brief AXPY for Vandermonde systems
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] alpha Scalar coefficient
 * @param[in] X Rhs Vandermonde system
 * @param[in] m Number of rows
 * @param[in] n Number of columns
 * @param[in] Y Lhs Vandermonde system
 * @param[in] stride (optional) stride in columns, default is \a m (row size)
 */
template <class T>
inline int vander_axpy(const T alpha, const T *X, const int m, const int n,
                       T *Y, const int stride = 0) {
  static lapack_int inc(1);

  if (alpha != T(0)) {
    lapack_int ld(stride <= 0 ? m : stride);
    if (ld == m) {
      lapack_int N(m * n);
      axpy(&N, const_cast<T *>(&alpha), const_cast<T *>(X), &inc, Y, &inc);
    } else {
      lapack_int m_(m);
      for (int j = 0; j < n; ++j)
        axpy(&m_, const_cast<T *>(&alpha), const_cast<T *>(X + j * ld), &inc,
             Y + ld * j, &inc);
    }
  }
  return 0;
}

/*!
 * @brief General interface for triangular solve used in accessing R in QRCP
 */
template <class T>
inline int rrqr_trisolve(const T *R, const int n, const int rank,
                         const int stride_r, const int nrhs, T *bs,
                         const int stride_bs, char trans) {
  if (nrhs <= 0) return 0;  // quick return if neccessary

  char       uplo = 'U', diag = 'N';  // configuration flags
  lapack_int info(0), rank_(rank), nrhs_(nrhs), ldr(stride_r), ldbs(stride_bs);
  if (rank_ < 0 || rank_ > n) rank_ = n;
  trtrs(&uplo, &trans, &diag, &rank_, &nrhs_, const_cast<T *>(R), &ldr, bs,
        &ldbs, &info);
  if (rank_ != n)
    for (int k = 0; k < nrhs; ++k)
      std::fill(bs + k * ldbs + rank_, bs + k * ldbs + n, T(0));
  return info;
}

/*!
 * @brief Perform backward substitution of R\bs
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] R Upper triangular matrix from rrqr_factor
 * @param[in] n Dimension of the triangular matrix
 * @param[in] rank Numerical rank of \a R, no larger than \a n
 * @param[in] stride_r Stride in column for \a R
 * @param[in] nrhs Number of right-hand sides
 * @param[in,out] bs Upon input, this stores the RHS values, on successful
 *                   return, this array will be overwritten by the solution.
 * @param[in] stride_bs Stride in column for \a bs
 * @return Zero for successful return, while negative values indicate error
 *         occurs in LAPACK's ?trtrs
 *
 * @sa rrqr_factor, rrqr_rtsolve
 */

template <class T>
inline int rrqr_rsolve(const T *R, const int n, const int rank,
                       const int stride_r, const int nrhs, T *bs,
                       const int stride_bs) {
  char no_tran = 'N';
  return rrqr_trisolve(R, n, rank, stride_r, nrhs, bs, stride_bs, no_tran);
}

/*!
 * @brief Perform forward subsitution of R'\bs
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] R Upper triangular matrix from rrqr_factor
 * @param[in] n Dimension of the triangular matrix
 * @param[in] rank Numerical rank of \a R
 * @param[in] stride_r Stride in column for \a R
 * @param[in] nrhs Number of right-hand sides
 * @param[in,out] bs Upon input, this stores the RHS values, on successful
 *                   return, this array will be overwritten by the solution.
 * @param[in] stride_bs Stride in column for \a bs
 * @return Zero for successful return, while negative values indicate error
 *         occurs in LAPACK's ?trtrs
 *
 * @sa rrqr_factor, rrqr_rsolve
 */
template <class T>
inline int rrqr_rtsolve(const T *R, const int n, const int rank,
                        const int stride_r, const int nrhs, T *bs,
                        const int stride_bs) {
  char tran = 'T';
  return rrqr_trisolve(R, n, rank, stride_r, nrhs, bs, stride_bs, tran);
}

/*!
 * @brief Common interface for accessing Q in QRCP
 */
template <class T>
inline int rrqr_apply_q(const T *QR, const int m, const int n, const int rank,
                        const int stride_qr, const int nrhs, T *bs,
                        const int stride_bs, T *work, const int work_size,
                        char trans) {
  if (nrhs <= 0) return 0;  // quick return if neccessary

  char       left = 'L';
  lapack_int m_(m), nrhs_(nrhs), rank_(rank), ldqr(stride_qr), ldbs(stride_bs),
      work_size_(work_size), info(0);
  if (rank_ < 0 || rank_ > n) rank_ = n;
  T *tau = const_cast<T *>(&QR[ldqr * n]);  // decode first entries in HS
  ormqr(&left, &trans, &m_, &nrhs_, &rank_, const_cast<T *>(QR), &ldqr, tau, bs,
        &ldbs, work, &work_size_, &info);
  return info;
}

/*!
 * @brief Implicitly form Q*b for Q arising from rrqr_factor
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] QR RRQR from \ref rrqr_factor
 * @param[in] m Number of rows in \a QR
 * @param[in] n Number of columns in \a QR
 * @param[in] rank Numerical rank of \a QR
 * @param[in] stride_qr Stride in columns for QR
 * @param[in] nrhs Number of right-hand sides
 * @param[in,out] bs Upon input, this stores the RHS values for the first
 *                   \a rank entries for each of the RHS vector. On successful
 *                   return, The first \a m entries will be overwritten by the
 *                   solutions. Therefore, the dimension of \a bs is m*nrhs.
 * @param[in] stride_bs Stride in columns for \a bs
 * @param[out] work Work space
 * @param[in] work_size Length of work space, returned by \ref query_work_size
 * @return Zero for successful return, while negative values indicate error
 *         occurs in LAPACK's ?ormqr
 *
 * @sa rrqr_qtmulti
 */
template <class T>
inline int rrqr_qmulti(const T *QR, const int m, const int n, const int rank,
                       const int stride_qr, const int nrhs, T *bs,
                       const int stride_bs, T *work, const int work_size) {
  char no_tran = 'N';
  return rrqr_apply_q(QR, m, n, rank, stride_qr, nrhs, bs, stride_bs, work,
                      work_size, no_tran);
}

/*!
 * @brief Implicitly form Q'*b for Q arising from rrqr_factor
 *
 * @tparam T Value type, e.g., \a double or \a float
 * @param[in] QR RRQR from \ref rrqr_factor
 * @param[in] m Number of rows in \a QR
 * @param[in] n Number of columns in \a QR
 * @param[in] rank Numerical rank of \a QR
 * @param[in] stride_qr Stride in columns for QR
 * @param[in] nrhs Number of right-hand sides
 * @param[in,out] bs Upon input, this stores the RHS values for the first
 *                   \a m entries for each of the RHS vector. On successful
 *                   return, The first \a rank entries will be overwritten by
 *                   the solutions. Therefore, the dimension of \a bs is m*nrhs.
 * @param[in] stride_bs Stride in columns for \a bs
 * @param[out] work Work space
 * @param[in] work_size Length of work space, returned by \ref query_work_size
 * @return Zero for successful return, while negative values indicate error
 *         occurs in LAPACK's ?ormqr
 *
 * @sa rrqr_qmulti
 */
template <class T>
inline int rrqr_qtmulti(const T *QR, const int m, const int n, const int rank,
                        const int stride_qr, const int nrhs, T *bs,
                        const int stride_bs, T *work, const int work_size) {
  char tran = 'T';
  return rrqr_apply_q(QR, m, n, rank, stride_qr, nrhs, bs, stride_bs, work,
                      work_size, tran);
}

/*!
 * @brief Vectorized assignment of one vector
 *
 * @note \a v must be aligned wrt WLS_ALING (default is 32) and satisfied
 *        restrict requirement.
 */
template <class T>
inline void vec_assign1(const T a, const int n,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, v);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) v[ld + j] = a;
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) v[i] = a;
}

/*!
 * @brief Vectorized assignment of two vectors
 *
 * @note \a v1 and \a v2 must be aligned wrt WLS_ALING (default is 32) and
 *        satisfied restrict requirement.
 */
template <class T>
inline void vec_assign2(const T a, const int n,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v1,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v2) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, v1);
  WLS_ASSERT_ALIGNED(T, v2);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) v1[ld + j] = a;
    for (int j = 0; j < ALIGN; ++j) v2[ld + j] = a;
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    v1[i] = a;
    v2[i] = a;
  }
}

/*!
 * @brief Vectorized assignment of three vectors
 *
 * @note \a v1, \a v2 and \a v3 must be aligned wrt WLS_ALING (default is 32)
 *        and satisfied restrict requirement.
 */
template <class T>
inline void vec_assign3(const T a, const int n,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v1,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v2,
                        T *WLS_RESTRICT WLS_ALIGN_ATTR v3) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, v1);
  WLS_ASSERT_ALIGNED(T, v2);
  WLS_ASSERT_ALIGNED(T, v3);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) v1[ld + j] = a;
    for (int j = 0; j < ALIGN; ++j) v2[ld + j] = a;
    for (int j = 0; j < ALIGN; ++j) v3[ld + j] = a;
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    v1[i] = a;
    v2[i] = a;
    v3[i] = a;
  }
}

/*!
 * @brief Perform vectorized y=x
 *
 */
template <class T>
inline void vec_copy11(const T *WLS_RESTRICT WLS_ALIGN_ATTR x, const int n,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, y);
  WLS_ASSERT_ALIGNED(const T, x);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y[ld + j] = x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) y[i] = x[i];
}

/*!
 * @brief Perform vectorized y1=x and y2=x
 *
 */
template <class T>
inline void vec_copy12(const T *WLS_RESTRICT WLS_ALIGN_ATTR x, const int n,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y1,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y2) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, y1);
  WLS_ASSERT_ALIGNED(T, y2);
  WLS_ASSERT_ALIGNED(const T, x);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y1[ld + j] = x[ld + j];
    for (int j = 0; j < ALIGN; ++j) y2[ld + j] = x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    y1[i] = x[i];
    y2[i] = x[i];
  }
}

/*!
 * @brief Perform vectorized y1=x, y2=x and y3=x
 *
 */
template <class T>
inline void vec_copy13(const T *WLS_RESTRICT WLS_ALIGN_ATTR x, const int n,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y1,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y2,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y3) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(T, y1);
  WLS_ASSERT_ALIGNED(T, y2);
  WLS_ASSERT_ALIGNED(T, y3);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y1[ld + j] = x[ld + j];
    for (int j = 0; j < ALIGN; ++j) y2[ld + j] = x[ld + j];
    for (int j = 0; j < ALIGN; ++j) y3[ld + j] = x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    y1[i] = x[i];
    y2[i] = x[i];
    y3[i] = x[i];
  }
}

/*!
 * @brief Perform vectorized y1=x1 and y2=x2
 *
 */
template <class T>
inline void vec_copy22(const T *WLS_RESTRICT WLS_ALIGN_ATTR x1,
                       const T *WLS_RESTRICT WLS_ALIGN_ATTR x2, const int n,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y1,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y2) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, y1);
  WLS_ASSERT_ALIGNED(T, y2);
  WLS_ASSERT_ALIGNED(const T, x1);
  WLS_ASSERT_ALIGNED(const T, x2);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y1[ld + j] = x1[ld + j];
    for (int j = 0; j < ALIGN; ++j) y2[ld + j] = x2[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    y1[i] = x1[i];
    y2[i] = x2[i];
  }
}

/*!
 * @brief Perform vectorized y1=x1, y2=x2 and y3=x3
 *
 */
template <class T>
inline void vec_copy33(const T *WLS_RESTRICT WLS_ALIGN_ATTR x1,
                       const T *WLS_RESTRICT WLS_ALIGN_ATTR x2,
                       const T *WLS_RESTRICT WLS_ALIGN_ATTR x3, const int n,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y1,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y2,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR y3) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x1);
  WLS_ASSERT_ALIGNED(const T, x2);
  WLS_ASSERT_ALIGNED(const T, x3);
  WLS_ASSERT_ALIGNED(T, y1);
  WLS_ASSERT_ALIGNED(T, y2);
  WLS_ASSERT_ALIGNED(T, y3);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y1[ld + j] = x1[ld + j];
    for (int j = 0; j < ALIGN; ++j) y2[ld + j] = x2[ld + j];
    for (int j = 0; j < ALIGN; ++j) y3[ld + j] = x3[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    y1[i] = x1[i];
    y2[i] = x2[i];
    y3[i] = x3[i];
  }
}

/*!
 * @brief Asigned the scaled vector to another one, i.e., y=alpha*x
 *
 */
template <class T>
inline void vec_scalar_multiply(const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
                                const int n, const T alpha,
                                T *WLS_RESTRICT WLS_ALIGN_ATTR y) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(T, y);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y[ld + j] = alpha * x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) y[i] = alpha * x[i];
}

/*!
 * @brief Scale in-place, i.e., x=alpha*x
 *
 */
template <class T>
inline void vec_scale(const int n, const T alpha,
                      T *WLS_RESTRICT WLS_ALIGN_ATTR x) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(T, x);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) x[ld + j] *= alpha;
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) x[i] *= alpha;
}

/*!
 * @brief component-wise scaling, i.e., y=x.*y
 *
 */
template <class T>
inline void vec_scale(const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
                      T *WLS_RESTRICT WLS_ALIGN_ATTR y, const int n) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(T, y);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) y[ld + j] *= x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) y[i] *= x[i];
}

/*!
 * @brief vector-wise product, i.e., z=x.*y
 *
 */
template <class T>
inline void vec_prod(const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
                     const T *WLS_RESTRICT WLS_ALIGN_ATTR y, const int n,
                     T *WLS_RESTRICT WLS_ALIGN_ATTR z) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(const T, y);
  WLS_ASSERT_ALIGNED(T, z);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) z[ld + j] = y[ld + j] * x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) z[i] = y[i] * x[i];
}

/*!
 * @brief Perfom dot/inner product of two real vectors, i.e., v=x^T*y
 *
 */
template <class T>
inline typename std::enable_if<std::is_floating_point<T>::value, T>::type
vec_dot(const int n, const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
        const T *WLS_RESTRICT WLS_ALIGN_ATTR y) {
  T v(0);
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(const T, y);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) v += y[ld + j] * x[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) v += y[i] * x[i];
  return v;
}

/*!
 * @brief Perform dot/inner product of two complex vectors, with conjugate of
 *        the first vector, i.e., v=x^H*y
 *
 * Notice that this routine works with any struct/class with public fields of
 * "re" (real part) and "im" (imaginary part).
 */
template <class T>
inline typename std::enable_if<!std::is_floating_point<T>::value, T>::type
vec_dot(const int n, const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
        const T *WLS_RESTRICT WLS_ALIGN_ATTR y) {
  T v;
  v.re = v.im = decltype(v.re)(0);
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(const T, y);
  // the inner loop should be unrolled and vectorized automatically
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j) {
      (v.re += x[ld + j].re * y[ld + j].re) += x[ld + j].im * y[ld + j].im;
      (v.im += x[ld + j].re * y[ld + j].im) -= x[ld + j].im * y[ld + j].re;
    }
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) {
    (v.re += x[i].re * y[i].re) += x[i].im * y[i].im;
    (v.im += x[i].re * y[i].im) -= x[i].im * y[i].re;
  }
  return v;
}

/*!
 * @brief Perform dot/inner product of two complex vectors with \a std::complex
 *
 * Notice that the built-in arithmatic operations with \a std::complex are
 * relatively slow, thus we internally use faster routine with struct of fields
 * "re" and "im".
 */
template <class T>
inline std::complex<T> vec_dot(
    const int n, const std::complex<T> *WLS_RESTRICT WLS_ALIGN_ATTR x,
    const std::complex<T> *WLS_RESTRICT WLS_ALIGN_ATTR y) {
  typedef struct {
    T re;
    T im;
  } __cast_type;
  auto v = vec_dot(n, reinterpret_cast<const WLS_ALIGN_ATTR __cast_type *>(x),
                   reinterpret_cast<const WLS_ALIGN_ATTR __cast_type *>(y));
  return *reinterpret_cast<std::complex<T> *>(&v);
}

/*!
 * @brief Perform vectorized AXPY, i.e., y=y+alpha*x
 *
 */
template <class T>
inline void vec_axpy(const T *WLS_RESTRICT WLS_ALIGN_ATTR x, const int n,
                     const T alpha, T *WLS_RESTRICT WLS_ALIGN_ATTR y,
                     const bool zero_y = false) {
  if (zero_y) return vec_scalar_multiply(x, n, alpha, y);
  if (alpha != T(0)) {
#if defined(__GNUC__) && __GNUC__ > 5
    static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
    const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
    WLS_ASSERT_ALIGNED(const T, x);
    WLS_ASSERT_ALIGNED(T, y);
    // the inner loop should be unrolled and vectorized automatically
    for (int i = 0; i < n0; ++i) {
      const int ld = i * ALIGN;
      for (int j = 0; j < ALIGN; ++j) y[ld + j] += alpha * x[ld + j];
    }
#else
    static constexpr int offset_start = 0;
#endif
    for (int i = offset_start; i < n; ++i) y[i] += alpha * x[i];
  }
}

/*!
 * @brief Perform vectorized AXBYPZ, i.e., z=z+alpha*x+beta*y
 *
 */
template <class T>
inline void vec_axbypz(const T *WLS_RESTRICT WLS_ALIGN_ATTR x,
                       const T *WLS_RESTRICT WLS_ALIGN_ATTR y, const int n,
                       const T alpha, const T beta,
                       T *WLS_RESTRICT WLS_ALIGN_ATTR z) {
#if defined(__GNUC__) && __GNUC__ > 5
  static constexpr int ALIGN = WLS_ALIGN / sizeof(T);
  const int            n0(n / ALIGN), offset_start(n0 * ALIGN);
  // the inner loop should be unrolled and vectorized automatically
  WLS_ASSERT_ALIGNED(const T, x);
  WLS_ASSERT_ALIGNED(const T, y);
  WLS_ASSERT_ALIGNED(T, z);
  for (int i = 0; i < n0; ++i) {
    const int ld = ALIGN * i;
    for (int j = 0; j < ALIGN; ++j)
      (z[ld + j] += alpha * x[ld + j]) += beta * y[ld + j];
  }
#else
  static constexpr int offset_start = 0;
#endif
  for (int i = offset_start; i < n; ++i) (z[i] += alpha * x[i]) += beta * y[i];
}

/*! @} */

}  // namespace wls

#endif  // WLS_LAPACK_HPP_
