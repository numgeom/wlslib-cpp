//
// wls_internal.cpp
//
// Code generation for function 'wls_internal'
//

// Include files
#include "wls_internal.h"
#include "wls_internal_types.h"
#include "coder_array.h"
#include "wls_lapack.hpp"
#include <cmath>
#include <cstring>
#include <stdio.h>

// Variable Definitions
namespace wls
{
  static const double dv[7]{ 333.33333333333331, 1000.0, 3333.3333333333335,
    10000.0, 100000.0, 1.0E+6, 1.0E+7 };

  static const char cv[128]{ '\x00', '\x01', '\x02', '\x03', '\x04', '\x05',
    '\x06', '\x07', '\x08', '	', '\x0a', '\x0b', '\x0c', '\x0d', '\x0e', '\x0f',
    '\x10', '\x11', '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18',
    '\x19', '\x1a', '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', ' ', '!', '\"', '#',
    '$', '%', '&', '\'', '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2',
    '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A',
    'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
    'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_',
    '`', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
    'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '{', '|', '}',
    '~', '\x7f' };

  static const signed char iv[9]{ 8, 9, 9, 10, 11, 11, 12, 12, 13 };

  static const signed char iv1[9]{ 14, 15, 15, 16, 17, 17, 18, 18, 19 };

  static const signed char iv2[9]{ 5, 6, 7, 8, 9, 10, 0, 0, 0 };

  static const signed char iv3[18]{ 7, 10, 6, 0, 8, 0, 6, 0, 5, 10, 9, 0, 8, 0,
    9, 0, 5, 7 };

  static const signed char iv4[9]{ 2, 3, 4, 0, 0, 0, 0, 0, 0 };

  static const signed char iv5[9]{ 3, 4, 2, 4, 2, 3, 0, 0, 0 };

  static const signed char iv6[9]{ 5, 6, 8, 6, 7, 9, 8, 9, 10 };
}

// Function Declarations
namespace wls
{
  static inline
  void compute_distances(const ::coder::array<double, 2U> &us, int
    npoints, boolean_T infnorm, ::coder::array<double, 1U> &ws);
  static inline
  double find_kth_shortest_dist(::coder::array<double, 1U> &arr, int k,
    int l, int r);
  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, int order, const ::coder::array<double, 1U> &weights, ::coder::array<
    double, 2U> &V);
  static inline
  void gen_vander(const double us_data[], const int us_size[2], int
    degree, ::coder::array<double, 2U> &V);
  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, int order, const double hs_inv_data[], const int hs_inv_size[2], ::
    coder::array<double, 2U> &V);
  static inline
  void gen_vander_1d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag);
  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const ::coder::array<double, 1U> &weights, ::coder::
    array<double, 2U> &V);
  static inline
  void gen_vander_2d(const double us_data[], int degree, ::coder::array<
    double, 2U> &V);
  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const double hs_inv_data[], const int hs_inv_size[2],
    ::coder::array<double, 2U> &V);
  static inline
  void gen_vander_2d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag);
  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const ::coder::array<double, 1U> &weights, ::coder::
    array<double, 2U> &V);
  static inline
  void gen_vander_3d(const double us_data[], int degree, ::coder::array<
    double, 2U> &V);
  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const double hs_inv_data[], const int hs_inv_size[2],
    ::coder::array<double, 2U> &V);
  static inline
  void gen_vander_3d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag);
  static inline
  int rrqr_factor(const ::coder::array<double, 2U> &A, double thres, int
    m, int n, ::coder::array<double, 2U> &QR, ::coder::array<int, 1U> &p, ::
    coder::array<double, 1U> &work, const ::coder::array<unsigned char, 2U> &dag);
  static inline
  void rrqr_qmulti(const ::coder::array<double, 2U> &QR, int m, int n,
    int rank, ::coder::array<double, 2U> &bs, int nrhs, ::coder::array<double,
    1U> &work);
  static inline
  void rrqr_rtsolve(const ::coder::array<double, 2U> &QR, int n, int rank,
    ::coder::array<double, 2U> &bs, int nrhs);
  static inline
  void wls_buhmann_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const ::coder::array<double, 1U> &params_sh, const ::
    coder::array<double, 2U> &params_pw, ::coder::array<double, 1U> &ws);
  static inline
  void wls_eno_weights(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 2U> &us_unscaled, const ::coder::
    array<double, 1U> &params_sh, const ::coder::array<double, 2U> &params_pw, ::
    coder::array<double, 1U> &ws);
  static inline
  void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const ::coder::array<double, 1U> &params_sh, const ::
    coder::array<double, 2U> &params_pw, ::coder::array<double, 1U> &ws);
  static inline
  void wls_resize(WlsObject *b_wls, int dim, int npoints, int degree, int
    order, boolean_T use_dag);
}

// Function Definitions
namespace wls
{
  static inline
  void compute_distances(const ::coder::array<double, 2U> &us, int
    npoints, boolean_T infnorm, ::coder::array<double, 1U> &ws)
  {
    //  Compute norms for each point using 2-norm or inf-norm
    if (!infnorm) {
      //  Compute 2-norm
      for (int i{0}; i < npoints; i++) {
        double d;
        double r2;
        int b_i;
        d = us[us.size(1) * i];
        r2 = d * d;
        b_i = us.size(1);
        for (int j{2}; j <= b_i; j++) {
          d = us[(j + us.size(1) * i) - 1];
          r2 += d * d;
        }

        ws[i] = std::sqrt(r2);
      }
    } else {
      //  Compute inf-norm for tensor-product
      for (int i{0}; i < npoints; i++) {
        double r;
        int b_i;
        r = std::abs(us[us.size(1) * i]);
        b_i = us.size(1);
        for (int j{2}; j <= b_i; j++) {
          double r1;
          r1 = std::abs(us[(j + us.size(1) * i) - 1]);
          if (r1 > r) {
            r = r1;
          }
        }

        ws[i] = r;
      }
    }
  }

  static inline
  double find_kth_shortest_dist(::coder::array<double, 1U> &arr, int k,
    int l, int r)
  {
    double dist;
    double val;
    int i;
    int j;

    //  Find the kth smallest number in arr(l:r).
    //
    //       [dist, arr] = find_kth_shortest_dist(arr, k);
    //       [dist, arr] = find_kth_shortest_dist(arr, k, start, end)
    if (k < l) {
      k = l;
    }

    if (k > r) {
      k = r;
    }

    val = arr[l - 1];
    i = l;
    j = r;
    while (i <= j) {
      double d;
      double d1;
      int exitg1;
      do {
        exitg1 = 0;
        d = arr[i - 1];
        if (d < val) {
          i++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      do {
        exitg1 = 0;
        d1 = arr[j - 1];
        if (d1 > val) {
          j--;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (i <= j) {
        arr[i - 1] = d1;
        arr[j - 1] = d;
        i++;
        j--;
      }
    }

    if (k <= j) {
      dist = find_kth_shortest_dist(arr, k, l, j);
    } else if (k >= i) {
      dist = find_kth_shortest_dist(arr, k, i, r);
    } else {
      dist = val;
    }

    return dist;
  }

  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, int order, const double hs_inv_data[], const int hs_inv_size[2], ::
    coder::array<double, 2U> &V)
  {
    //  Wrapper function for computing confluent Vandermonde matrix in 1D, 2D, or 3D.
    //
    //     V = gen_vander(us)
    //     V = gen_vander(us, npoints)
    //     V = gen_vander(us, npoints, degree)
    //     V = gen_vander(us, npoints, degree, order)
    //     V = gen_vander(us, npoints, degree, order, weights)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv, V)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-d, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix
    //              Use -1, -2, and -4 for grad, Laplacian, and bi-Laplacian.
    //     weights: Weights for all points (n-by-1, where n>=1; use zeros(0,1)
    //              or omit it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //  The order argument must passed in with coder.ignoreConst(<expr>) in
    //  the caller order to avoid generation of local buffers.
    //
    //  See also
    //     gen_vander_1d, gen_vander_2d, gen_vander_3d
    //  Compute Vandermonde system
    switch (us.size(1)) {
     case 1:
      {
        double h_inv_;
        int b_npoints;
        int i;
        int i1;
        int i2;
        int iPnt;
        int nrblks;
        int nrows;
        int r;
        int stride;
        boolean_T b;
        boolean_T b1;
        boolean_T flag;
        b_npoints = npoints - 1;

        //  Generate (confluent) Vandermonde matrix in 1D.
        //
        //     V = gen_vander_1d(us)
        //     V = gen_vander_1d(us, npoints)
        //     V = gen_vander_1d(us, npoints, degree, order)
        //     V = gen_vander_1d(us, npoints, degree, order, weights)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V, stride)
        //
        //  Parameters
        //  ----------
        //     us:      Local coordinates of points (n-by-1, where n>=npoints)
        //     npoints: Number of points. Use 0 for default (size(us, 1))
        //     degree:  Maximum degree of monomials (default is 2)
        //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
        //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively
        //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
        //              it to use unit weights)
        //     h_inv:   Inverse of radius for scaling rows in CVM (size 1-by-0 or 1-by-1)
        //     V:       Vandermonde matrix (must be preallocated if present at input)
        //     stride:  number of rows in each row block in V (0 for default)
        //
        //  Returns
        //  -------
        //     V:      confluent Vandermonde matrix
        //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
        //
        //  Notes
        //  -----
        //     It is more efficient to provide weights here to incorporate them when
        //     constructing the Vandermonde matrix (linear-time overhead in npoints)
        //     than scaling the matrix afterwards (linear in npoints*nmonomials).
        //
        //     Entries in each row of V are in ascending degrees. For  the confluent
        //     Vandermonde matrix, row blocks are in increasing order of derivatives.
        //
        //     For example, if order==1, V has the following content:
        //         weights(1) * [1, u1, u1^2, u1^3, u1^4, ...]
        //         weights(2) * [1, u2, u2^2, u2^3, u2^4, ...]
        //         ...
        //         h_inv * weights(1) * [0, 1, 2u1, 3u1^2, 4u1^3, ...]
        //         h_inv * weights(2) * [0, 1, 2u2, 3u2^2, 4u2^3, ...]
        //         ...
        //
        //  See also gen_vander_2d, gen_vander_3d
        //  Handle input arguments
        if (npoints == 0) {
          b_npoints = us.size(0) - 1;
        } else {
          flag = (npoints <= us.size(0));

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "Input us is too small.");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Input us is too small.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        }

        flag = (degree >= 0);

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(flag, "Degree must be nonnegative");

#else //MATLAB_MEX_FILE

        if (!flag) {
          fprintf(stderr, "Degree must be nonnegative\n");
          fflush(stderr);
        }

        assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
          flag = true;
        } else {
          flag = false;
        }

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

        if (!flag) {
          fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
          fflush(stderr);
        }

        assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        if (hs_inv_size[1] == 0) {
          h_inv_ = 1.0;
        } else {
          h_inv_ = hs_inv_data[0];
        }

        stride = us.size(0);
        nrblks = order + 1;

        //  Number of row blocks
        switch (order) {
         case -1:
          nrblks = 2;
          break;

         case -2:
          nrblks = 3;
          break;

         case -4:
          nrblks = 4;
          break;
        }

        nrows = us.size(0) * nrblks;
        if ((V.size(1) != nrows) || (V.size(0) != degree + 1)) {
          V.set_size(degree + 1, nrows);
        }

        //  Compute rows corresponding to function values
        if (degree != 0) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          i = us.size(1) * us.size(0);
          i1 = 0;
          for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
            if (b1 || (iPnt >= i)) {
              i1 = 0;
              b = true;
            } else if (b) {
              b = false;
              i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              i2 = us.size(1) * us.size(0) - 1;
              if (i1 > MAX_int32_T - us.size(1)) {
                i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i1 += us.size(1);
                if (i1 > i2) {
                  i1 -= i2;
                }
              }
            }

            V[iPnt] = 1.0;
            V[iPnt + V.size(1)] = us[i1];
          }
        } else {
          for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
            V[iPnt] = 1.0;
          }
        }

        if (0 > order) {
          i = order;
        } else {
          i = 0;
        }

        i = (degree + i) + 1;
        for (int ii{2}; ii <= i; ii++) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          i1 = us.size(1) * us.size(0);
          i2 = 0;
          for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
            if (b1 || (iPnt >= i1)) {
              i2 = 0;
              b = true;
            } else if (b) {
              b = false;
              i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              int i3;
              i3 = us.size(1) * us.size(0) - 1;
              if (i2 > MAX_int32_T - us.size(1)) {
                i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i2 += us.size(1);
                if (i2 > i3) {
                  i2 -= i3;
                }
              }
            }

            V[iPnt + V.size(1) * (ii - 1)] = V[iPnt + V.size(1) * (ii - 2)] *
              us[i2];
          }
        }

        //  Add row blocks corresponding to kth derivatives
        r = us.size(0);
        if (order >= 0) {
          for (int k{0}; k < order; k++) {
            int j;
            for (j = 0; j <= k; j++) {
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = k + 1; j <= degree; j++) {
              double s;
              s = h_inv_ * static_cast<double>(j);
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                i = r + iPnt;
                V[i + V.size(1) * j] = V[(i - stride) + V.size(1) * (j - 1)] * s;
              }
            }

            r += stride;
          }
        } else {
          double s;
          int j;

          //     %% computing negative orders
          if (-order > 2) {
            i = 2;
          } else {
            i = -order;
          }

          for (int k{0}; k < i; k++) {
            for (j = 0; j <= k; j++) {
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = k + 1; j <= degree; j++) {
              s = h_inv_ * static_cast<double>(j);
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                i1 = r + iPnt;
                V[i1 + V.size(1) * j] = V[(i1 - stride) + V.size(1) * (j - 1)] *
                  s;
              }
            }

            r += stride;
          }

          //     %% Calculate Biharmonic if order = -4
          if (order == -4) {
            for (j = 0; j < 4; j++) {
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = 2; j <= degree; j++) {
              s = h_inv_ * static_cast<double>(j);
              for (iPnt = 0; iPnt <= b_npoints; iPnt++) {
                i = r + iPnt;
                V[i + V.size(1) * j] = V[(i - stride) + V.size(1) * (j - 2)] * s
                  * (s - 1.0);
              }
            }
          }
        }
      }
      break;

     case 2:
      gen_vander_2d(us, npoints, degree, order, hs_inv_data, hs_inv_size, V);
      break;

     default:
      gen_vander_3d(us, npoints, degree, order, hs_inv_data, hs_inv_size, V);
      break;
    }
  }

  static inline
  void gen_vander(const double us_data[], const int us_size[2], int
    degree, ::coder::array<double, 2U> &V)
  {
    //  Wrapper function for computing confluent Vandermonde matrix in 1D, 2D, or 3D.
    //
    //     V = gen_vander(us)
    //     V = gen_vander(us, npoints)
    //     V = gen_vander(us, npoints, degree)
    //     V = gen_vander(us, npoints, degree, order)
    //     V = gen_vander(us, npoints, degree, order, weights)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv, V)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-d, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix
    //              Use -1, -2, and -4 for grad, Laplacian, and bi-Laplacian.
    //     weights: Weights for all points (n-by-1, where n>=1; use zeros(0,1)
    //              or omit it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //  The order argument must passed in with coder.ignoreConst(<expr>) in
    //  the caller order to avoid generation of local buffers.
    //
    //  See also
    //     gen_vander_1d, gen_vander_2d, gen_vander_3d
    //  Compute Vandermonde system
    switch (us_size[1]) {
     case 1:
      {
        int i;

        //  Generate (confluent) Vandermonde matrix in 1D.
        //
        //     V = gen_vander_1d(us)
        //     V = gen_vander_1d(us, npoints)
        //     V = gen_vander_1d(us, npoints, degree, order)
        //     V = gen_vander_1d(us, npoints, degree, order, weights)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V, stride)
        //
        //  Parameters
        //  ----------
        //     us:      Local coordinates of points (n-by-1, where n>=npoints)
        //     npoints: Number of points. Use 0 for default (size(us, 1))
        //     degree:  Maximum degree of monomials (default is 2)
        //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
        //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively
        //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
        //              it to use unit weights)
        //     h_inv:   Inverse of radius for scaling rows in CVM (size 1-by-0 or 1-by-1)
        //     V:       Vandermonde matrix (must be preallocated if present at input)
        //     stride:  number of rows in each row block in V (0 for default)
        //
        //  Returns
        //  -------
        //     V:      confluent Vandermonde matrix
        //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
        //
        //  Notes
        //  -----
        //     It is more efficient to provide weights here to incorporate them when
        //     constructing the Vandermonde matrix (linear-time overhead in npoints)
        //     than scaling the matrix afterwards (linear in npoints*nmonomials).
        //
        //     Entries in each row of V are in ascending degrees. For  the confluent
        //     Vandermonde matrix, row blocks are in increasing order of derivatives.
        //
        //     For example, if order==1, V has the following content:
        //         weights(1) * [1, u1, u1^2, u1^3, u1^4, ...]
        //         weights(2) * [1, u2, u2^2, u2^3, u2^4, ...]
        //         ...
        //         h_inv * weights(1) * [0, 1, 2u1, 3u1^2, 4u1^3, ...]
        //         h_inv * weights(2) * [0, 1, 2u2, 3u2^2, 4u2^3, ...]
        //         ...
        //
        //  See also gen_vander_2d, gen_vander_3d
        //  Handle input arguments
        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(true, "Input us is too small.");

#else //MATLAB_MEX_FILE

        assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(true, "Degree must be nonnegative");

#else //MATLAB_MEX_FILE

        assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(true, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

        assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        //  Number of row blocks
        if ((V.size(1) != 1) || (V.size(0) != degree + 1)) {
          V.set_size(degree + 1, 1);
        }

        //  Compute rows corresponding to function values
        V[0] = 1.0;
        V[V.size(1)] = us_data[0];
        i = degree + 1;
        for (int ii{2}; ii <= i; ii++) {
          V[V.size(1) * (ii - 1)] = V[V.size(1) * (ii - 2)] * us_data[0];
        }

        //  Add row blocks corresponding to kth derivatives
      }
      break;

     case 2:
      gen_vander_2d(us_data, degree, V);
      break;

     default:
      gen_vander_3d(us_data, degree, V);
      break;
    }
  }

  static inline
  void gen_vander(const ::coder::array<double, 2U> &us, int npoints, int
    degree, int order, const ::coder::array<double, 1U> &weights, ::coder::array<
    double, 2U> &V)
  {
    //  Wrapper function for computing confluent Vandermonde matrix in 1D, 2D, or 3D.
    //
    //     V = gen_vander(us)
    //     V = gen_vander(us, npoints)
    //     V = gen_vander(us, npoints, degree)
    //     V = gen_vander(us, npoints, degree, order)
    //     V = gen_vander(us, npoints, degree, order, weights)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander(us, npoints, degree, order, weights, hs_inv, V)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-d, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix
    //              Use -1, -2, and -4 for grad, Laplacian, and bi-Laplacian.
    //     weights: Weights for all points (n-by-1, where n>=1; use zeros(0,1)
    //              or omit it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //  The order argument must passed in with coder.ignoreConst(<expr>) in
    //  the caller order to avoid generation of local buffers.
    //
    //  See also
    //     gen_vander_1d, gen_vander_2d, gen_vander_3d
    //  Compute Vandermonde system
    switch (us.size(1)) {
     case 1:
      {
        int i;
        int i1;
        int i2;
        int iPnt;
        int nrblks;
        int nrows;
        int r;
        int stride;
        boolean_T b;
        boolean_T b1;
        boolean_T flag;

        //  Generate (confluent) Vandermonde matrix in 1D.
        //
        //     V = gen_vander_1d(us)
        //     V = gen_vander_1d(us, npoints)
        //     V = gen_vander_1d(us, npoints, degree, order)
        //     V = gen_vander_1d(us, npoints, degree, order, weights)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V)
        //     V = gen_vander_1d(us, npoints, degree, order, weights, h_inv, V, stride)
        //
        //  Parameters
        //  ----------
        //     us:      Local coordinates of points (n-by-1, where n>=npoints)
        //     npoints: Number of points. Use 0 for default (size(us, 1))
        //     degree:  Maximum degree of monomials (default is 2)
        //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
        //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively
        //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
        //              it to use unit weights)
        //     h_inv:   Inverse of radius for scaling rows in CVM (size 1-by-0 or 1-by-1)
        //     V:       Vandermonde matrix (must be preallocated if present at input)
        //     stride:  number of rows in each row block in V (0 for default)
        //
        //  Returns
        //  -------
        //     V:      confluent Vandermonde matrix
        //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
        //
        //  Notes
        //  -----
        //     It is more efficient to provide weights here to incorporate them when
        //     constructing the Vandermonde matrix (linear-time overhead in npoints)
        //     than scaling the matrix afterwards (linear in npoints*nmonomials).
        //
        //     Entries in each row of V are in ascending degrees. For  the confluent
        //     Vandermonde matrix, row blocks are in increasing order of derivatives.
        //
        //     For example, if order==1, V has the following content:
        //         weights(1) * [1, u1, u1^2, u1^3, u1^4, ...]
        //         weights(2) * [1, u2, u2^2, u2^3, u2^4, ...]
        //         ...
        //         h_inv * weights(1) * [0, 1, 2u1, 3u1^2, 4u1^3, ...]
        //         h_inv * weights(2) * [0, 1, 2u2, 3u2^2, 4u2^3, ...]
        //         ...
        //
        //  See also gen_vander_2d, gen_vander_3d
        //  Handle input arguments
        flag = (npoints <= us.size(0));

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(flag, "Input us is too small.");

#else //MATLAB_MEX_FILE

        if (!flag) {
          fprintf(stderr, "Input us is too small.\n");
          fflush(stderr);
        }

        assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        flag = (degree >= 0);

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(flag, "Degree must be nonnegative");

#else //MATLAB_MEX_FILE

        if (!flag) {
          fprintf(stderr, "Degree must be nonnegative\n");
          fflush(stderr);
        }

        assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
          flag = true;
        } else {
          flag = false;
        }

        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

        if (!flag) {
          fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
          fflush(stderr);
        }

        assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        stride = us.size(0);
        nrblks = order + 1;

        //  Number of row blocks
        switch (order) {
         case -1:
          nrblks = 2;
          break;

         case -2:
          nrblks = 3;
          break;

         case -4:
          nrblks = 4;
          break;
        }

        nrows = us.size(0) * nrblks;
        if ((V.size(1) != nrows) || (V.size(0) != degree + 1)) {
          V.set_size(degree + 1, nrows);
        }

        //  Compute rows corresponding to function values
        if (weights.size(0) == 0) {
          if (degree != 0) {
            b = true;
            b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
            i = us.size(1) * us.size(0);
            i1 = 0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              if (b1 || (iPnt >= i)) {
                i1 = 0;
                b = true;
              } else if (b) {
                b = false;
                i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i2 = us.size(1) * us.size(0) - 1;
                if (i1 > MAX_int32_T - us.size(1)) {
                  i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
                } else {
                  i1 += us.size(1);
                  if (i1 > i2) {
                    i1 -= i2;
                  }
                }
              }

              V[iPnt] = 1.0;
              V[iPnt + V.size(1)] = us[i1];
            }
          } else {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[iPnt] = 1.0;
            }
          }
        } else if (degree != 0) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          i = us.size(1) * us.size(0);
          i1 = 0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            if (b1 || (iPnt >= i)) {
              i1 = 0;
              b = true;
            } else if (b) {
              b = false;
              i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              i2 = us.size(1) * us.size(0) - 1;
              if (i1 > MAX_int32_T - us.size(1)) {
                i1 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i1 += us.size(1);
                if (i1 > i2) {
                  i1 -= i2;
                }
              }
            }

            V[iPnt] = weights[iPnt];
            V[iPnt + V.size(1)] = us[i1] * weights[iPnt];
          }
        } else {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[iPnt] = weights[iPnt];
          }
        }

        if (0 > order) {
          i = order;
        } else {
          i = 0;
        }

        i = (degree + i) + 1;
        for (int ii{2}; ii <= i; ii++) {
          b = true;
          b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
          i1 = us.size(1) * us.size(0);
          i2 = 0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            if (b1 || (iPnt >= i1)) {
              i2 = 0;
              b = true;
            } else if (b) {
              b = false;
              i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
            } else {
              int i3;
              i3 = us.size(1) * us.size(0) - 1;
              if (i2 > MAX_int32_T - us.size(1)) {
                i2 = iPnt % us.size(0) * us.size(1) + iPnt / us.size(0);
              } else {
                i2 += us.size(1);
                if (i2 > i3) {
                  i2 -= i3;
                }
              }
            }

            V[iPnt + V.size(1) * (ii - 1)] = V[iPnt + V.size(1) * (ii - 2)] *
              us[i2];
          }
        }

        //  Add row blocks corresponding to kth derivatives
        r = us.size(0);
        if (order >= 0) {
          for (int k{0}; k < order; k++) {
            int j;
            for (j = 0; j <= k; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = k + 1; j <= degree; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i = r + iPnt;
                V[i + V.size(1) * j] = V[(i - stride) + V.size(1) * (j - 1)] *
                  static_cast<double>(j);
              }
            }

            r += stride;
          }
        } else {
          int j;

          //     %% computing negative orders
          if (-order > 2) {
            i = 2;
          } else {
            i = -order;
          }

          for (int k{0}; k < i; k++) {
            for (j = 0; j <= k; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = k + 1; j <= degree; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = r + iPnt;
                V[i1 + V.size(1) * j] = V[(i1 - stride) + V.size(1) * (j - 1)] *
                  static_cast<double>(j);
              }
            }

            r += stride;
          }

          //     %% Calculate Biharmonic if order = -4
          if (order == -4) {
            for (j = 0; j < 4; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(r + iPnt) + V.size(1) * j] = 0.0;
              }
            }

            for (j = 2; j <= degree; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i = r + iPnt;
                V[i + V.size(1) * j] = V[(i - stride) + V.size(1) * (j - 2)] *
                  static_cast<double>(j) * (static_cast<double>(j) - 1.0);
              }
            }
          }
        }
      }
      break;

     case 2:
      gen_vander_2d(us, npoints, degree, order, weights, V);
      break;

     default:
      gen_vander_3d(us, npoints, degree, order, weights, V);
      break;
    }
  }

  static inline
  void gen_vander_1d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag)
  {
    int u0;
    int u1;
    boolean_T flag;

    //  Build a dag for Vandermonde matrix in 1D.
    //
    //     dag = gen_vander_1d_dag(degree)
    //     dag = gen_vander_1d_dag(degree, dag)
    //
    //  Parameters
    //  ----------
    //     degree:  Maximum degree of monomials (default is 2).
    //
    //  Returns
    //  -------
    //     dag:     A direct acyclic graph stored in a 1-by-(degree+2) M-array
    //        (in column major, or a (degree+2)-by-1 M-array if row major).
    //        Each monomial points to its "child" that is one degree higher.
    //        Each of its entry stores the difference between the index
    //        of its child and its own index. The last entry is used store
    //        a signature, so that it can be recomputed when needed.
    //
    //     dag_colMajor is an optional output in column-major for testing.
    //
    //  See also gen_vander_1d, gen_vander_2d_dag, gen_vander_3d_dag, rrqr_trunc
    //  Handle input arguments
    if ((dag.size(0) == 0) || (dag.size(1) == 0)) {
      u1 = 0;
    } else {
      u0 = dag.size(0);
      u1 = dag.size(1);
      if (u0 >= u1) {
        u1 = u0;
      }
    }

    flag = (u1 >= degree + 2);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "DAG must be preallocated.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "DAG must be preallocated.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    for (int i{0}; i < degree; i++) {
      u1 = dag.size(0);
      dag[i % u1 * dag.size(1) + i / u1] = 1U;
    }

    u1 = dag.size(0);
    dag[degree % u1 * dag.size(1) + degree / u1] = 0U;

    //  a leaf has no child
    //  Use last entry as signature
    u1 = dag.size(1) * dag.size(0) - 1;
    u0 = dag.size(0);
    dag[u1 % u0 * dag.size(1) + u1 / u0] = static_cast<unsigned char>(degree +
      127);
  }

  static inline
  void gen_vander_2d(const double us_data[], int degree, ::coder::array<
    double, 2U> &V)
  {
    int c;
    int i;
    int ncols;

    //  Generate generalized/confluent Vandermonde matrix in 2D.
    //
    //     V = gen_vander_2d(us)
    //     V = gen_vander_2d(us, npoints)
    //     V = gen_vander_2d(us, npoints, degree, order)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-2, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal triangle,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal triangle.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2]
    //        ...
    //
    //     For degree=-2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2, u1^2*v1, u1*v1^2, u1^2*v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2, u2^2*v2, u2*v2^2, u2^2*v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0,    2u1*v1,  v1^2,    2u1*v1^2]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0,    2u2*v2,  v2^2,    2u2*v2^2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1,    u1^2,    2u1*v1,  2u1^2*v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2,    u2^2,    2u2*v1,  2u2^2*v2]
    //        ...
    //
    //  See also gen_vander_1d, gen_vander_3d
    //  Handle input arguments
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "Input us is too small.");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  Number of row blocks
    ncols = (degree + 1) * (degree + 2) / 2;

    //  Allocate storage for V
    if ((V.size(1) != 1) || (V.size(0) != ncols)) {
      V.set_size(ncols, 1);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    V[0] = 1.0;
    V[V.size(1)] = us_data[0];
    V[V.size(1) * 2] = us_data[1];
    c = 3;
    if (degree < 0) {
      i = -degree;
    } else {
      i = degree;
    }

    for (int deg{2}; deg <= i; deg++) {
      for (int j{0}; j < deg; j++) {
        V[V.size(1) * c] = V[V.size(1) * (c - deg)] * us_data[0];
        c++;
      }

      V[V.size(1) * c] = V[V.size(1) * ((c - deg) - 1)] * us_data[1];
      c++;
    }

    //  Compute the bi-degree terms if degree<0
    //  compute higher order confluent Vandermonde matrix blocks incrementally
  }

  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const ::coder::array<double, 1U> &weights, ::coder::
    array<double, 2U> &V)
  {
    int b_degree;
    int c;
    int deg;
    int i;
    int i1;
    int iPnt;
    int j;
    int k;
    int ncols;
    int nrblks;
    int nrows;
    int stride;
    int x;
    boolean_T flag;

    //  Generate generalized/confluent Vandermonde matrix in 2D.
    //
    //     V = gen_vander_2d(us)
    //     V = gen_vander_2d(us, npoints)
    //     V = gen_vander_2d(us, npoints, degree, order)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-2, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal triangle,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal triangle.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2]
    //        ...
    //
    //     For degree=-2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2, u1^2*v1, u1*v1^2, u1^2*v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2, u2^2*v2, u2*v2^2, u2^2*v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0,    2u1*v1,  v1^2,    2u1*v1^2]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0,    2u2*v2,  v2^2,    2u2*v2^2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1,    u1^2,    2u1*v1,  2u1^2*v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2,    u2^2,    2u2*v1,  2u2^2*v2]
    //        ...
    //
    //  See also gen_vander_1d, gen_vander_3d
    //  Handle input arguments
    flag = (npoints <= us.size(0));

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Input us is too small.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Input us is too small.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    stride = us.size(0);
    nrblks = (order + 1) * (order + 2) / 2;

    //  Number of row blocks
    switch (order) {
     case -1:
      nrblks = 3;
      break;

     case -2:
      nrblks = 5;
      break;

     case -4:
      if (degree > 0) {
        nrblks = 8;
      } else {
        nrblks = 11;
      }
      break;
    }

    if (degree >= 0) {
      ncols = (degree + 1) * (degree + 2) / 2;
    } else {
      ncols = (1 - degree) * (1 - degree);
    }

    //  Allocate storage for V
    nrows = us.size(0) * nrblks;
    if ((V.size(1) != nrows) || (V.size(0) != ncols)) {
      V.set_size(ncols, nrows);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    if (weights.size(0) == 0) {
      if (degree != 0) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
          V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
          V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
        }
      } else {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
        }
      }
    } else if (degree != 0) {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
      }
    } else {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
      }
    }

    c = 3;
    if (0 > order) {
      i = order;
    } else {
      i = 0;
    }

    if (-degree > 0) {
      b_degree = 1;
    } else if (-degree < 0) {
      b_degree = -1;
    } else {
      b_degree = 0;
    }

    x = b_degree * i;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    if (x < 0) {
      x = 0;
    }

    i1 = b_degree - x;
    for (deg = 2; deg <= i1; deg++) {
      for (j = 0; j < deg; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }

      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
    }

    //  Compute the bi-degree terms if degree<0
    i1 = -degree;
    i = 1 - i;
    for (deg = i1; deg >= i; deg--) {
      for (k = 0; k < deg; k++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }
    }

    //  compute higher order confluent Vandermonde matrix blocks incrementally
    if (order != 0) {
      double scaleu;
      double scalev;
      int i2;
      int len;
      int offset;

      //  This is an optimized version of update_vander_ordern for first-order CVM
      flag = (degree != 0);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  Compute derivative with respect to u
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        i = stride + iPnt;
        V[i] = 0.0;
        V[i + V.size(1)] = V[iPnt];
        V[i + V.size(1) * 2] = 0.0;
      }

      c = 3;
      if (0 > order + 1) {
        i = order + 1;
      } else {
        i = 0;
      }

      if (-degree > 0) {
        b_degree = 1;
      } else if (-degree < 0) {
        b_degree = -1;
      } else {
        b_degree = 0;
      }

      x = b_degree * i;
      if (degree < 0) {
        b_degree = -degree;
      } else {
        b_degree = degree;
      }

      if (x < 0) {
        x = 0;
      }

      i2 = b_degree - x;
      for (deg = 2; deg <= i2; deg++) {
        scaleu = deg;
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] *
            static_cast<double>(deg);
        }

        c++;
        for (j = 0; j <= deg - 2; j++) {
          scaleu--;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)]
              * scaleu;
          }

          c++;
        }

        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(stride + iPnt) + V.size(1) * c] = 0.0;
        }

        c++;
      }

      //  Compute the bi-degree terms if degree<0
      for (len = i1; len >= i; len--) {
        scaleu = 1 - degree;
        for (k = 0; k < len; k++) {
          scaleu--;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - len)]
              * scaleu;
          }

          c++;
        }
      }

      //  Compute derivative with respect to v
      offset = us.size(0) + us.size(0);
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        i2 = offset + iPnt;
        V[i2] = 0.0;
        V[i2 + V.size(1)] = 0.0;
        V[i2 + V.size(1) * 2] = V[iPnt];
      }

      c = 3;
      if (-degree > 0) {
        b_degree = 1;
      } else if (-degree < 0) {
        b_degree = -1;
      } else {
        b_degree = 0;
      }

      if (0 > order + 1) {
        i2 = order + 1;
      } else {
        i2 = 0;
      }

      x = b_degree * i2;
      if (degree < 0) {
        b_degree = -degree;
      } else {
        b_degree = degree;
      }

      if (x < 0) {
        x = 0;
      }

      i2 = b_degree - x;
      for (deg = 2; deg <= i2; deg++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(offset + iPnt) + V.size(1) * c] = 0.0;
        }

        c++;
        for (j = 0; j < deg; j++) {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg)
              - 1)] * (static_cast<double>(j) + 1.0);
          }

          c++;
        }
      }

      //  Compute the bi-degree terms if degree<0
      deg = -degree;
      for (len = i1; len >= i; len--) {
        deg++;
        scalev = (deg + degree) - 1;
        for (k = 0; k < len; k++) {
          scalev++;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * ((c - len)
              - 1)] * scalev;
          }

          c++;
        }
      }

      //     %% compute regular orders if order > 0
      if (order > 0) {
        for (int dd{2}; dd <= order; dd++) {
          int col;
          int offset_prev;
          int row;

          //  Compute order-N CVM row blocks from order-(N-1) CVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          x = dd * (dd + 1) / 2;
          offset = x * stride;
          offset_prev = (dd - 1) * dd / 2 * stride;

          //  Compute derivative with respect to u
          for (int b_i{0}; b_i < dd; b_i++) {
            //  Initialize block to zero
            for (col = 0; col < x; col++) {
              i = offset + 1;
              i2 = offset + npoints;
              for (row = i; row <= i2; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = x;
            if (degree < 0) {
              i = -degree;
            } else {
              i = degree;
            }

            for (deg = dd; deg <= i; deg++) {
              scaleu = deg + 1;
              i2 = deg - 1;
              for (j = 0; j <= i2; j++) {
                scaleu--;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - deg)] * scaleu;
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            //  Compute the bi-degree terms if degree<0
            for (len = i1; len >= 0; len--) {
              scaleu = 1 - degree;
              for (k = 0; k < len; k++) {
                scaleu--;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - len)] * scaleu;
                }

                c++;
              }
            }

            offset += stride;
            offset_prev += stride;
          }

          //  Compute derivative with respect to v
          //  Initialize block to zero
          for (col = 0; col < x; col++) {
            i = offset + 1;
            i2 = offset + npoints;
            for (row = i; row <= i2; row++) {
              V[(row + V.size(1) * col) - 1] = 0.0;
            }
          }

          c = x;
          if (degree < 0) {
            i = -degree;
          } else {
            i = degree;
          }

          for (deg = dd; deg <= i; deg++) {
            i2 = dd - 1;
            for (j = 0; j <= i2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (j = dd; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride) +
                  iPnt) + V.size(1) * ((c - deg) - 1)] * static_cast<double>(j);
              }

              c++;
            }
          }

          //  Compute the bi-degree terms if degree<0
          deg = -degree;
          for (len = i1; len >= 0; len--) {
            deg++;
            scalev = (deg + degree) - 1;
            for (k = 0; k < len; k++) {
              scalev++;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride) +
                  iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
              }

              c++;
            }
          }
        }
      } else if (order < -1) {
        //         %% compute efficient laplacian and bi-laplacian
        switch (order) {
         case -2:
          {
            int col;
            int offset_prev;
            int row;

            //  Compute order-N CVM row blocks from order-(N-1) CVM.
            flag = (degree != 0);

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag, "");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr, "Runtime assertion error.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            offset = 3 * us.size(0);
            offset_prev = us.size(0);

            //  Compute derivative with respect to u
            if (degree < 0) {
              i = -degree;
            } else {
              i = degree;
            }

            //  Initialize block to zero
            i2 = offset + 1;
            x = offset + npoints;
            for (col = 0; col < 3; col++) {
              for (row = i2; row <= x; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = 3;
            for (deg = 2; deg <= i; deg++) {
              scaleu = deg + 1;
              i2 = deg - 1;
              for (j = 0; j <= i2; j++) {
                scaleu--;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - deg)] * scaleu;
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            //  Compute the bi-degree terms if degree<0
            for (len = i1; len >= -2; len--) {
              scaleu = 1 - degree;
              for (k = 0; k < len; k++) {
                scaleu--;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - len)] * scaleu;
                }

                c++;
              }
            }

            offset += us.size(0);

            //  Compute derivative with respect to v
            //  Initialize block to zero
            offset_prev = (us.size(0) + us.size(0)) + us.size(0);
            i = offset + 1;
            i2 = offset + npoints;
            for (col = 0; col < 3; col++) {
              for (row = i; row <= i2; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = 3;
            if (degree < 0) {
              i = -degree;
            } else {
              i = degree;
            }

            for (deg = 2; deg <= i; deg++) {
              for (j = 0; j < 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (j = 2; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride)
                    + iPnt) + V.size(1) * ((c - deg) - 1)] * static_cast<double>
                    (j);
                }

                c++;
              }
            }

            //  Compute the bi-degree terms if degree<0
            deg = -degree;
            for (len = i1; len >= -2; len--) {
              deg++;
              scalev = (deg + degree) - 1;
              for (k = 0; k < len; k++) {
                scalev++;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride)
                    + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                }

                c++;
              }
            }
          }
          break;

         case -4:
          {
            if (degree > 0) {
              int col;
              int offset_prev;
              int row;

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(true, "");

#else //MATLAB_MEX_FILE

              assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 3 * us.size(0);

              //  Compute derivative with respect to u
              //  Initialize block to zero
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= degree; deg++) {
                scaleu = deg + 1;
                i = deg - 1;
                for (j = 0; j <= i; j++) {
                  scaleu--;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(us.size(0) + iPnt) +
                      V.size(1) * (c - deg)] * scaleu;
                  }

                  c++;
                }

                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              //  Compute the bi-degree terms if degree<0
              offset += us.size(0);

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev = (us.size(0) + us.size(0)) + us.size(0);
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= degree; deg++) {
                for (j = 0; j < 2; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 2; j <= deg; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] *
                      static_cast<double>(j);
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              //  function for the Bilaplacian of P elements
              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);
              flag = (degree >= 4);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              //  compute du^4 and du^2*dv^2
              for (int terms{0}; terms < 2; terms++) {
                i = offset + 1;
                i1 = offset + npoints;
                for (col = 0; col < 10; col++) {
                  for (row = i; row <= i1; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 10;
                for (deg = 4; deg <= degree; deg++) {
                  scaleu = deg + 1;
                  i = deg - 1;
                  for (j = 0; j <= i; j++) {
                    scaleu--;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * ((c - (deg << 1)) + 1)] * scaleu * (scaleu
                        - 1.0);
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                offset += stride;
                offset_prev += stride;
              }

              //  compute dv^4
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 10; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 10;
              for (deg = 4; deg <= degree; deg++) {
                for (j = 0; j < 4; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 4; j <= deg; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - (deg << 1)) - 1)] *
                      static_cast<double>(j) * (static_cast<double>(j) - 1.0);
                  }

                  c++;
                }
              }
            } else {
              int b_i;
              int col;
              int offset_prev;
              int row;

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 3 * us.size(0);
              offset_prev = us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              //  Initialize block to zero
              i2 = offset + 1;
              x = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i2; row <= x; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= i; deg++) {
                scaleu = deg + 1;
                i2 = deg - 1;
                for (j = 0; j <= i2; j++) {
                  scaleu--;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                      + V.size(1) * (c - deg)] * scaleu;
                  }

                  c++;
                }

                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              //  Compute the bi-degree terms if degree<0
              for (len = i1; len >= -2; len--) {
                scaleu = 1 - degree;
                for (k = 0; k < len; k++) {
                  scaleu--;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                      + V.size(1) * (c - len)] * scaleu;
                  }

                  c++;
                }
              }

              offset += us.size(0);

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev = (us.size(0) + us.size(0)) + us.size(0);
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 2; deg <= i; deg++) {
                for (j = 0; j < 2; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 2; j <= deg; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] *
                      static_cast<double>(j);
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -2; len--) {
                deg++;
                scalev = (deg + degree) - 1;
                for (k = 0; k < len; k++) {
                  scalev++;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (b_i = 0; b_i < 2; b_i++) {
                //  Initialize block to zero
                i2 = offset + 1;
                x = offset + npoints;
                for (col = 0; col < 6; col++) {
                  for (row = i2; row <= x; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 6;
                for (deg = 3; deg <= i; deg++) {
                  scaleu = deg + 1;
                  i2 = deg - 1;
                  for (j = 0; j <= i2; j++) {
                    scaleu--;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - deg)] * scaleu;
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                //  Compute the bi-degree terms if degree<0
                for (len = i1; len >= -3; len--) {
                  scaleu = 1 - degree;
                  for (k = 0; k < len; k++) {
                    scaleu--;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - len)] * scaleu;
                    }

                    c++;
                  }
                }

                offset += stride;
                offset_prev += stride;
              }

              //  Compute derivative with respect to v
              //  Initialize block to zero
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 6; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 6;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 3; deg <= i; deg++) {
                for (j = 0; j < 3; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 3; j <= deg; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] *
                      static_cast<double>(j);
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -3; len--) {
                deg++;
                scalev = (deg + degree) - 1;
                for (k = 0; k < len; k++) {
                  scalev++;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = us.size(0) << 3;
              offset_prev = 5 * us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (b_i = 0; b_i < 2; b_i++) {
                //  Initialize block to zero
                i2 = offset + 1;
                x = offset + npoints;
                for (col = 0; col < 10; col++) {
                  for (row = i2; row <= x; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 10;
                for (deg = 4; deg <= i; deg++) {
                  scaleu = deg + 1;
                  i2 = deg - 1;
                  for (j = 0; j <= i2; j++) {
                    scaleu--;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - deg)] * scaleu;
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                //  Compute the bi-degree terms if degree<0
                for (len = i1; len >= -4; len--) {
                  scaleu = 1 - degree;
                  for (k = 0; k < len; k++) {
                    scaleu--;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - len)] * scaleu;
                    }

                    c++;
                  }
                }

                offset += stride;
                offset_prev += stride;
              }

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev += us.size(0);
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 10; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 10;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 4; deg <= i; deg++) {
                for (j = 0; j < 4; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 4; j <= deg; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] *
                      static_cast<double>(j);
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -4; len--) {
                deg++;
                scalev = (deg + degree) - 1;
                for (k = 0; k < len; k++) {
                  scalev++;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }
            }
          }
          break;

         default:
          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(false, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

          fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
          fflush(stderr);
          assert(false);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          break;
        }
      }
    }
  }

  static inline
  void gen_vander_2d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const double hs_inv_data[], const int hs_inv_size[2],
    ::coder::array<double, 2U> &V)
  {
    double hs_inv__idx_0;
    double hs_inv__idx_1;
    int b_degree;
    int c;
    int deg;
    int i;
    int i1;
    int iPnt;
    int j;
    int k;
    int ncols;
    int nrblks;
    int nrows;
    int stride;
    int x;
    boolean_T flag;

    //  Generate generalized/confluent Vandermonde matrix in 2D.
    //
    //     V = gen_vander_2d(us)
    //     V = gen_vander_2d(us, npoints)
    //     V = gen_vander_2d(us, npoints, degree, order)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_2d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-2, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal triangle,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal triangle.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2]
    //        ...
    //
    //     For degree=-2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, u1^2, u1*v1, v1^2, u1^2*v1, u1*v1^2, u1^2*v1^2]
    //        weights(2) * [1, u2, v2, u2^2, u2*v2, v2^2, u2^2*v2, u2*v2^2, u2^2*v2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  2u1,  v1,    0,    2u1*v1,  v1^2,    2u1*v1^2]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  2u2,  v2,    0,    2u2*v2,  v2^2,    2u2*v2^2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    u1,  2v1,    u1^2,    2u1*v1,  2u1^2*v1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    u2,  2v2,    u2^2,    2u2*v1,  2u2^2*v2]
    //        ...
    //
    //  See also gen_vander_1d, gen_vander_3d
    //  Handle input arguments
    if (npoints == 0) {
      npoints = us.size(0);
    } else {
      flag = (npoints <= us.size(0));

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "Input us is too small.");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Input us is too small.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (hs_inv_size[1] == 0) {
      hs_inv__idx_0 = 1.0;
      hs_inv__idx_1 = 1.0;
    } else {
      hs_inv__idx_0 = hs_inv_data[0];
      hs_inv__idx_1 = hs_inv_data[1];
    }

    stride = us.size(0);
    nrblks = (order + 1) * (order + 2) / 2;

    //  Number of row blocks
    switch (order) {
     case -1:
      nrblks = 3;
      break;

     case -2:
      nrblks = 5;
      break;

     case -4:
      if (degree > 0) {
        nrblks = 8;
      } else {
        nrblks = 11;
      }
      break;
    }

    if (degree >= 0) {
      ncols = (degree + 1) * (degree + 2) / 2;
    } else {
      ncols = (1 - degree) * (1 - degree);
    }

    //  Allocate storage for V
    nrows = us.size(0) * nrblks;
    if ((V.size(1) != nrows) || (V.size(0) != ncols)) {
      V.set_size(ncols, nrows);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    if (degree != 0) {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
      }
    } else {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }

    c = 3;
    if (0 > order) {
      i = order;
    } else {
      i = 0;
    }

    if (-degree > 0) {
      b_degree = 1;
    } else if (-degree < 0) {
      b_degree = -1;
    } else {
      b_degree = 0;
    }

    x = b_degree * i;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    if (x < 0) {
      x = 0;
    }

    i1 = b_degree - x;
    for (deg = 2; deg <= i1; deg++) {
      for (j = 0; j < deg; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }

      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
    }

    //  Compute the bi-degree terms if degree<0
    i1 = -degree;
    i = 1 - i;
    for (deg = i1; deg >= i; deg--) {
      for (k = 0; k < deg; k++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] * us[us.size
            (1) * iPnt];
        }

        c++;
      }
    }

    //  compute higher order confluent Vandermonde matrix blocks incrementally
    if (order != 0) {
      double scaleu;
      double scalev;
      int i2;
      int len;
      int offset;

      //  This is an optimized version of update_vander_ordern for first-order CVM
      flag = (degree != 0);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  Compute derivative with respect to u
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        i = stride + iPnt;
        V[i] = 0.0;
        V[i + V.size(1)] = V[iPnt] * hs_inv__idx_0;
        V[i + V.size(1) * 2] = 0.0;
      }

      c = 3;
      if (0 > order + 1) {
        i = order + 1;
      } else {
        i = 0;
      }

      if (-degree > 0) {
        b_degree = 1;
      } else if (-degree < 0) {
        b_degree = -1;
      } else {
        b_degree = 0;
      }

      x = b_degree * i;
      if (degree < 0) {
        b_degree = -degree;
      } else {
        b_degree = degree;
      }

      if (x < 0) {
        x = 0;
      }

      i2 = b_degree - x;
      for (deg = 2; deg <= i2; deg++) {
        scaleu = static_cast<double>(deg) * hs_inv__idx_0;
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)] *
            scaleu;
        }

        c++;
        for (j = 0; j <= deg - 2; j++) {
          scaleu -= hs_inv__idx_0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - deg)]
              * scaleu;
          }

          c++;
        }

        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(stride + iPnt) + V.size(1) * c] = 0.0;
        }

        c++;
      }

      //  Compute the bi-degree terms if degree<0
      for (len = i1; len >= i; len--) {
        scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
        for (k = 0; k < len; k++) {
          scaleu -= hs_inv__idx_0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - len)]
              * scaleu;
          }

          c++;
        }
      }

      //  Compute derivative with respect to v
      offset = us.size(0) + us.size(0);
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        i2 = offset + iPnt;
        V[i2] = 0.0;
        V[i2 + V.size(1)] = 0.0;
        V[i2 + V.size(1) * 2] = V[iPnt] * hs_inv__idx_1;
      }

      c = 3;
      if (-degree > 0) {
        b_degree = 1;
      } else if (-degree < 0) {
        b_degree = -1;
      } else {
        b_degree = 0;
      }

      if (0 > order + 1) {
        i2 = order + 1;
      } else {
        i2 = 0;
      }

      x = b_degree * i2;
      if (degree < 0) {
        b_degree = -degree;
      } else {
        b_degree = degree;
      }

      if (x < 0) {
        x = 0;
      }

      i2 = b_degree - x;
      for (deg = 2; deg <= i2; deg++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[(offset + iPnt) + V.size(1) * c] = 0.0;
        }

        c++;
        for (j = 0; j < deg; j++) {
          scalev = (static_cast<double>(j) + 1.0) * hs_inv__idx_1;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * ((c - deg)
              - 1)] * scalev;
          }

          c++;
        }
      }

      //  Compute the bi-degree terms if degree<0
      deg = -degree;
      for (len = i1; len >= i; len--) {
        deg++;
        scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
        for (k = 0; k < len; k++) {
          scalev += hs_inv__idx_1;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * ((c - len)
              - 1)] * scalev;
          }

          c++;
        }
      }

      //     %% compute regular orders if order > 0
      if (order > 0) {
        for (b_degree = 2; b_degree <= order; b_degree++) {
          int col;
          int offset_prev;
          int row;

          //  Compute order-N CVM row blocks from order-(N-1) CVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          offset = 3 * stride;
          offset_prev = stride;

          //  Compute derivative with respect to u
          if (degree < 0) {
            i = -degree;
          } else {
            i = degree;
          }

          for (int b_i{0}; b_i < 2; b_i++) {
            //  Initialize block to zero
            i2 = offset + 1;
            x = offset + npoints;
            for (col = 0; col < 3; col++) {
              for (row = i2; row <= x; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = 3;
            for (deg = 2; deg <= i; deg++) {
              scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
              i2 = deg - 1;
              for (j = 0; j <= i2; j++) {
                scaleu -= hs_inv__idx_0;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - deg)] * scaleu;
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            //  Compute the bi-degree terms if degree<0
            for (len = i1; len >= 0; len--) {
              scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
              for (k = 0; k < len; k++) {
                scaleu -= hs_inv__idx_0;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - len)] * scaleu;
                }

                c++;
              }
            }

            offset += stride;
            offset_prev += stride;
          }

          //  Compute derivative with respect to v
          //  Initialize block to zero
          i = offset + 1;
          i2 = offset + npoints;
          for (col = 0; col < 3; col++) {
            for (row = i; row <= i2; row++) {
              V[(row + V.size(1) * col) - 1] = 0.0;
            }
          }

          c = 3;
          if (degree < 0) {
            i = -degree;
          } else {
            i = degree;
          }

          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j < 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (j = 2; j <= deg; j++) {
              scalev = static_cast<double>(j) * hs_inv__idx_1;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride) +
                  iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
              }

              c++;
            }
          }

          //  Compute the bi-degree terms if degree<0
          deg = -degree;
          for (len = i1; len >= 0; len--) {
            deg++;
            scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
            for (k = 0; k < len; k++) {
              scalev += hs_inv__idx_1;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride) +
                  iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
              }

              c++;
            }
          }
        }
      } else if (order < -1) {
        //         %% compute efficient laplacian and bi-laplacian
        switch (order) {
         case -2:
          {
            int col;
            int offset_prev;
            int row;

            //  Compute order-N CVM row blocks from order-(N-1) CVM.
            flag = (degree != 0);

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag, "");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr, "Runtime assertion error.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            offset = 3 * us.size(0);
            offset_prev = us.size(0);

            //  Compute derivative with respect to u
            if (degree < 0) {
              i = -degree;
            } else {
              i = degree;
            }

            //  Initialize block to zero
            i2 = offset + 1;
            x = offset + npoints;
            for (col = 0; col < 3; col++) {
              for (row = i2; row <= x; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = 3;
            for (deg = 2; deg <= i; deg++) {
              scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
              i2 = deg - 1;
              for (j = 0; j <= i2; j++) {
                scaleu -= hs_inv__idx_0;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - deg)] * scaleu;
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            //  Compute the bi-degree terms if degree<0
            for (len = i1; len >= -2; len--) {
              scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
              for (k = 0; k < len; k++) {
                scaleu -= hs_inv__idx_0;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt) +
                    V.size(1) * (c - len)] * scaleu;
                }

                c++;
              }
            }

            offset += us.size(0);

            //  Compute derivative with respect to v
            //  Initialize block to zero
            offset_prev = (us.size(0) + us.size(0)) + us.size(0);
            i = offset + 1;
            i2 = offset + npoints;
            for (col = 0; col < 3; col++) {
              for (row = i; row <= i2; row++) {
                V[(row + V.size(1) * col) - 1] = 0.0;
              }
            }

            c = 3;
            if (degree < 0) {
              i = -degree;
            } else {
              i = degree;
            }

            for (deg = 2; deg <= i; deg++) {
              for (j = 0; j < 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (j = 2; j <= deg; j++) {
                scalev = static_cast<double>(j) * hs_inv__idx_1;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride)
                    + iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
                }

                c++;
              }
            }

            //  Compute the bi-degree terms if degree<0
            deg = -degree;
            for (len = i1; len >= -2; len--) {
              deg++;
              scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
              for (k = 0; k < len; k++) {
                scalev += hs_inv__idx_1;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev - stride)
                    + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                }

                c++;
              }
            }
          }
          break;

         case -4:
          {
            if (degree > 0) {
              int col;
              int offset_prev;
              int row;

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(true, "");

#else //MATLAB_MEX_FILE

              assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 3 * us.size(0);

              //  Compute derivative with respect to u
              //  Initialize block to zero
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= degree; deg++) {
                scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
                i = deg - 1;
                for (j = 0; j <= i; j++) {
                  scaleu -= hs_inv__idx_0;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(us.size(0) + iPnt) +
                      V.size(1) * (c - deg)] * scaleu;
                  }

                  c++;
                }

                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              //  Compute the bi-degree terms if degree<0
              offset += us.size(0);

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev = (us.size(0) + us.size(0)) + us.size(0);
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= degree; deg++) {
                for (j = 0; j < 2; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 2; j <= deg; j++) {
                  scalev = static_cast<double>(j) * hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              //  function for the Bilaplacian of P elements
              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);
              flag = (degree >= 4);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              //  compute du^4 and du^2*dv^2
              for (int terms{0}; terms < 2; terms++) {
                i = offset + 1;
                i1 = offset + npoints;
                for (col = 0; col < 10; col++) {
                  for (row = i; row <= i1; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 10;
                for (deg = 4; deg <= degree; deg++) {
                  scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
                  i = deg - 1;
                  for (j = 0; j <= i; j++) {
                    scaleu -= hs_inv__idx_0;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * ((c - (deg << 1)) + 1)] * scaleu * (scaleu
                        - hs_inv__idx_0);
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                offset += stride;
                offset_prev += stride;
              }

              //  compute dv^4
              i = offset + 1;
              i1 = offset + npoints;
              for (col = 0; col < 10; col++) {
                for (row = i; row <= i1; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 10;
              for (deg = 4; deg <= degree; deg++) {
                for (j = 0; j < 4; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 4; j <= deg; j++) {
                  scalev = static_cast<double>(j) * (hs_inv__idx_1 *
                    hs_inv__idx_1);
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - (deg << 1)) - 1)] *
                      scalev * (static_cast<double>(j) - 1.0);
                  }

                  c++;
                }
              }
            } else {
              int b_i;
              int col;
              int offset_prev;
              int row;

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 3 * us.size(0);
              offset_prev = us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              //  Initialize block to zero
              i2 = offset + 1;
              x = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i2; row <= x; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              for (deg = 2; deg <= i; deg++) {
                scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
                i2 = deg - 1;
                for (j = 0; j <= i2; j++) {
                  scaleu -= hs_inv__idx_0;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                      + V.size(1) * (c - deg)] * scaleu;
                  }

                  c++;
                }

                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              //  Compute the bi-degree terms if degree<0
              for (len = i1; len >= -2; len--) {
                scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
                for (k = 0; k < len; k++) {
                  scaleu -= hs_inv__idx_0;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                      + V.size(1) * (c - len)] * scaleu;
                  }

                  c++;
                }
              }

              offset += us.size(0);

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev = (us.size(0) + us.size(0)) + us.size(0);
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 3; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 3;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 2; deg <= i; deg++) {
                for (j = 0; j < 2; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 2; j <= deg; j++) {
                  scalev = static_cast<double>(j) * hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -2; len--) {
                deg++;
                scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
                for (k = 0; k < len; k++) {
                  scalev += hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (b_i = 0; b_i < 2; b_i++) {
                //  Initialize block to zero
                i2 = offset + 1;
                x = offset + npoints;
                for (col = 0; col < 6; col++) {
                  for (row = i2; row <= x; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 6;
                for (deg = 3; deg <= i; deg++) {
                  scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
                  i2 = deg - 1;
                  for (j = 0; j <= i2; j++) {
                    scaleu -= hs_inv__idx_0;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - deg)] * scaleu;
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                //  Compute the bi-degree terms if degree<0
                for (len = i1; len >= -3; len--) {
                  scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
                  for (k = 0; k < len; k++) {
                    scaleu -= hs_inv__idx_0;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - len)] * scaleu;
                    }

                    c++;
                  }
                }

                offset += stride;
                offset_prev += stride;
              }

              //  Compute derivative with respect to v
              //  Initialize block to zero
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 6; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 6;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 3; deg <= i; deg++) {
                for (j = 0; j < 3; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 3; j <= deg; j++) {
                  scalev = static_cast<double>(j) * hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -3; len--) {
                deg++;
                scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
                for (k = 0; k < len; k++) {
                  scalev += hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute order-N CVM row blocks from order-(N-1) CVM.
              flag = (degree != 0);

              //  Throw error if condition false
              //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

              mxAssert(flag, "");

#else //MATLAB_MEX_FILE

              if (!flag) {
                fprintf(stderr, "Runtime assertion error.\n");
                fflush(stderr);
              }

              assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

              offset = us.size(0) << 3;
              offset_prev = 5 * us.size(0);

              //  Compute derivative with respect to u
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (b_i = 0; b_i < 2; b_i++) {
                //  Initialize block to zero
                i2 = offset + 1;
                x = offset + npoints;
                for (col = 0; col < 10; col++) {
                  for (row = i2; row <= x; row++) {
                    V[(row + V.size(1) * col) - 1] = 0.0;
                  }
                }

                c = 10;
                for (deg = 4; deg <= i; deg++) {
                  scaleu = static_cast<double>(deg + 1) * hs_inv__idx_0;
                  i2 = deg - 1;
                  for (j = 0; j <= i2; j++) {
                    scaleu -= hs_inv__idx_0;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - deg)] * scaleu;
                    }

                    c++;
                  }

                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                //  Compute the bi-degree terms if degree<0
                for (len = i1; len >= -4; len--) {
                  scaleu = static_cast<double>(1 - degree) * hs_inv__idx_0;
                  for (k = 0; k < len; k++) {
                    scaleu -= hs_inv__idx_0;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(offset + iPnt) + V.size(1) * c] = V[(offset_prev + iPnt)
                        + V.size(1) * (c - len)] * scaleu;
                    }

                    c++;
                  }
                }

                offset += stride;
                offset_prev += stride;
              }

              //  Compute derivative with respect to v
              //  Initialize block to zero
              offset_prev += us.size(0);
              i = offset + 1;
              i2 = offset + npoints;
              for (col = 0; col < 10; col++) {
                for (row = i; row <= i2; row++) {
                  V[(row + V.size(1) * col) - 1] = 0.0;
                }
              }

              c = 10;
              if (degree < 0) {
                i = -degree;
              } else {
                i = 0;
              }

              for (deg = 4; deg <= i; deg++) {
                for (j = 0; j < 4; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                }

                for (j = 4; j <= deg; j++) {
                  scalev = static_cast<double>(j) * hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - deg) - 1)] * scalev;
                  }

                  c++;
                }
              }

              //  Compute the bi-degree terms if degree<0
              deg = -degree;
              for (len = i1; len >= -4; len--) {
                deg++;
                scalev = static_cast<double>((deg + degree) - 1) * hs_inv__idx_1;
                for (k = 0; k < len; k++) {
                  scalev += hs_inv__idx_1;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((offset_prev -
                      stride) + iPnt) + V.size(1) * ((c - len) - 1)] * scalev;
                  }

                  c++;
                }
              }
            }
          }
          break;

         default:
          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(false, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

          fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
          fflush(stderr);
          assert(false);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          break;
        }
      }
    }
  }

  static inline
  void gen_vander_2d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag)
  {
    int c;
    int deg;
    int i;
    int j;

    //  Build a dag for Vandermonde matrix in 2D.
    //
    //     dag = gen_vander_2d_dag(degree)
    //     dag = gen_vander_2d_dag(degree, dag)
    //
    //  Parameters
    //  ----------
    //     degree:  Maximum degree of monomials. Use negative for tensor-product monomials.
    //
    //  Returns
    //  -------
    //     dag:     A direct acyclic graph stored in a 2-by-(#monomials+1) M-array
    //        (in column major, or a (#monomials+1)-by-2 M-array if row major).
    //        Each column represents a monomial (x^i*y^j), and it column stores
    //        the offsets to the indices of its "child" monomials
    //        (i.e., x^(i+1)*y^j and x^i*y^(j+1)). The last entry is used store
    //        a signature, so that it can be recomputed when needed.
    //
    //     dag_colMajor is an optional output in column-major for testing.
    //
    //  See also gen_vander_2d, gen_vander_1d_dag, gen_vander_3d_dag, rrqr_trunc
    //  Handle input arguments
    if (degree != 0) {
      dag[0] = 1U;

      //  x-child
      dag[1] = 2U;

      //  y-child
    } else {
      dag[0] = 0U;
      dag[1] = 0U;

      //  No children
    }

    c = 2;
    if (degree < 0) {
      i = -degree;
    } else {
      i = degree;
    }

    for (deg = 2; deg <= i; deg++) {
      dag[dag.size(1) * ((c - deg) + 1)] = static_cast<unsigned char>(deg);

      //  x-child
      c += 2;
      for (j = 2; j <= deg; j++) {
        dag[dag.size(1) * (((c + j) - deg) - 2)] = static_cast<unsigned char>
          (deg);

        //  x-child
        dag[dag.size(1) * (((c + j) - deg) - 3) + 1] = static_cast<unsigned char>
          (deg + 1);

        //  y-child
      }

      c = (c + deg) - 1;
      dag[dag.size(1) * ((c - deg) - 1) + 1] = static_cast<unsigned char>(deg +
        1);

      //  y-child
    }

    //  Set the children of last row to zero
    if (degree > 0) {
      c -= degree;
      for (j = 0; j <= degree; j++) {
        dag[dag.size(1) * c] = 0U;
        dag[dag.size(1) * c + 1] = 0U;

        //  no children
        c++;
      }
    } else if (degree < 0) {
      //  Compute the bi-degree terms if degree<0
      i = -degree;
      for (deg = i; deg >= 1; deg--) {
        j = c - deg;
        dag[dag.size(1) * j] = 0U;

        //  no x-child
        dag[dag.size(1) * j + 1] = static_cast<unsigned char>(deg + 1);

        //  y-child
        dag[dag.size(1) * (j + 1)] = static_cast<unsigned char>(deg);

        //  x-child
        c++;
        for (int k{2}; k <= deg; k++) {
          j = (c + k) - deg;
          dag[dag.size(1) * (j - 1)] = static_cast<unsigned char>(deg);

          //  x-child
          dag[dag.size(1) * (j - 2) + 1] = static_cast<unsigned char>(deg + 1);

          //  y-child
        }

        c = (c + deg) - 1;
        dag[dag.size(1) * (c - deg) + 1] = 0U;

        //  no y-child
      }

      dag[dag.size(1) * c] = 0U;
      dag[dag.size(1) * c + 1] = 0U;
    }

    //  Use last entry as signature
    i = dag.size(1) * dag.size(0) - 1;
    j = dag.size(0);
    dag[i % j * dag.size(1) + i / j] = static_cast<unsigned char>(degree + 127);
  }

  static inline
  void gen_vander_3d(const double us_data[], int degree, ::coder::array<
    double, 2U> &V)
  {
    int c;
    int d;
    int i;
    int ncols;

    //  Generate generalized/confluent Vandermonde matrix in 3D.
    //
    //     V = gen_vander_3d(us)
    //     V = gen_vander_3d(us, npoints)
    //     V = gen_vander_3d(us, npoints, degree, order)
    //     V = gen_vander_3d(us, npoints, degree, order, weights)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-3, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal tetrahedron,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal tetrahedron.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1^2, u1*v1, v1^2, u1*w1, v1*w1, w1^2]
    //        weights(2) * [1, u2, v2, w2, u2^2, u2*v2, v2^2, u2*w2, v2*w2, w2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  2u1,  v1,    0,   w1,     0,     0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  2u2,  v2,    0,   w2,     0,     0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    0,  u1,  2v1,   0,     w1,    0]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    0,  u2,  2v2,   0,     w2,    0]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,    0,   0,    0,   u1,   v1,     2w1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,    0,   0,    0,   u2,   v2,     2w2]
    //
    //     For degree=-1 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1*v1, u1*w1, v1*w1, u1*v1*w1]
    //        weights(2) * [1, u2, v2, w2, u2*v2, u2*w2, v2*w2, u2*v2*w2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  v1,   w1,   0,     v1*w1]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  v2,   w2,   0,     v2*w2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,  u1,   0,    w1,    u1*w1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,  u2,   0,    w2,    u2*w2]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,   0,   u1,   v1,    u1*v1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,   0,   u2,   v2,    u2*v2]
    //
    //  See also gen_vander_1d, gen_vander_2d
    //  Handle input arguments
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;

    //  Allocate storage for V
    if ((V.size(1) != 1) || (V.size(0) != ncols)) {
      V.set_size(ncols, 1);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    V[0] = 1.0;
    V[V.size(1)] = us_data[0];
    V[V.size(1) * 2] = us_data[1];
    V[V.size(1) * 3] = us_data[2];
    c = 4;
    d = 4;
    if (degree < 0) {
      i = -degree;
    } else {
      i = degree;
    }

    for (int deg{2}; deg <= i; deg++) {
      int j;

      //  Within each level, use convention of Pascal triangle with x^deg at peak
      for (j = 0; j < deg; j++) {
        V[V.size(1) * c] = V[V.size(1) * ((c - d) + 1)] * us_data[0];
        c++;
      }

      V[V.size(1) * c] = V[V.size(1) * (c - d)] * us_data[1];
      c++;
      for (j = 0; j <= d - 2; j++) {
        V[V.size(1) * c] = V[V.size(1) * ((c - d) - deg)] * us_data[2];
        c++;
      }

      d = (d + deg) + 1;
    }

    //  Compute the tri-degree terms if degree<0
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

  }

  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const double hs_inv_data[], const int hs_inv_size[2],
    ::coder::array<double, 2U> &V)
  {
    double hs_inv_idx_0;
    double hs_inv_idx_1;
    double hs_inv_idx_2;
    int b_degree;
    int c;
    int cornerTriangle;
    int counterBottomRow;
    int d;
    int deg;
    int excess;
    int i;
    int i1;
    int iPnt;
    int j;
    int k;
    int maxLayers;
    int nTermsInLayer;
    int nTermsInPrevLayer;
    int ncols;
    int nrblks;
    int nrows;
    int p;
    int stride;
    int x;
    boolean_T flag;

    //  Generate generalized/confluent Vandermonde matrix in 3D.
    //
    //     V = gen_vander_3d(us)
    //     V = gen_vander_3d(us, npoints)
    //     V = gen_vander_3d(us, npoints, degree, order)
    //     V = gen_vander_3d(us, npoints, degree, order, weights)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-3, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal tetrahedron,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal tetrahedron.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1^2, u1*v1, v1^2, u1*w1, v1*w1, w1^2]
    //        weights(2) * [1, u2, v2, w2, u2^2, u2*v2, v2^2, u2*w2, v2*w2, w2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  2u1,  v1,    0,   w1,     0,     0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  2u2,  v2,    0,   w2,     0,     0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    0,  u1,  2v1,   0,     w1,    0]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    0,  u2,  2v2,   0,     w2,    0]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,    0,   0,    0,   u1,   v1,     2w1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,    0,   0,    0,   u2,   v2,     2w2]
    //
    //     For degree=-1 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1*v1, u1*w1, v1*w1, u1*v1*w1]
    //        weights(2) * [1, u2, v2, w2, u2*v2, u2*w2, v2*w2, u2*v2*w2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  v1,   w1,   0,     v1*w1]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  v2,   w2,   0,     v2*w2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,  u1,   0,    w1,    u1*w1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,  u2,   0,    w2,    u2*w2]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,   0,   u1,   v1,    u1*v1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,   0,   u2,   v2,    u2*v2]
    //
    //  See also gen_vander_1d, gen_vander_2d
    //  Handle input arguments
    if (npoints == 0) {
      npoints = us.size(0);
    } else {
      flag = (npoints <= us.size(0));

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (hs_inv_size[1] == 0) {
      hs_inv_idx_0 = 1.0;
      hs_inv_idx_1 = 1.0;
      hs_inv_idx_2 = 1.0;
    } else {
      hs_inv_idx_0 = hs_inv_data[0];
      hs_inv_idx_1 = hs_inv_data[1];
      hs_inv_idx_2 = hs_inv_data[2];
    }

    stride = us.size(0);
    nrblks = (order + 1) * (order + 2) * (order + 3) / 6;
    switch (order) {
     case -1:
      nrblks = 4;
      break;

     case -2:
      nrblks = 7;
      break;

     case -4:
      if (degree > 0) {
        nrblks = 13;
      } else {
        nrblks = 19;
      }
      break;
    }

    if (degree >= 0) {
      ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
    } else {
      ncols = (1 - degree) * (1 - degree) * (1 - degree);
    }

    //  Allocate storage for V
    nrows = us.size(0) * nrblks;
    if ((V.size(1) != nrows) || (V.size(0) != ncols)) {
      V.set_size(ncols, nrows);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    if (degree != 0) {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
        V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
      }
    } else {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = 1.0;
      }
    }

    c = 4;
    d = 3;
    if (0 > order) {
      i = order;
    } else {
      i = 0;
    }

    if (-degree > 0) {
      b_degree = 1;
    } else if (-degree < 0) {
      b_degree = -1;
    } else {
      b_degree = 0;
    }

    x = b_degree * i;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    if (x < 0) {
      x = 0;
    }

    i1 = b_degree - x;
    for (deg = 2; deg <= i1; deg++) {
      //  Within each level, use convention of Pascal triangle with x^deg at peak
      for (j = 0; j < deg; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)] * us[us.size(1)
            * iPnt];
        }

        c++;
      }

      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - d) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
      for (j = 0; j < d; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
            us[us.size(1) * iPnt + 2];
        }

        c++;
      }

      d = (d + deg) + 1;
    }

    //  Compute the tri-degree terms if degree<0
    if (degree < 0) {
      deg = -degree;
      maxLayers = -degree * 3 + i;

      // max number of layers needed in the Pascal tetrahedron
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      nTermsInLayer = d;

      // initializing number of elements in layer
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      i = 1 - degree;
      for (p = i; p <= maxLayers; p++) {
        int gap;

        //  Within each level, x^deg is at the peak of Pascal triangle
        //  implicitly calculating number of elements in corner Pascal triangles
        cornerTriangle = (cornerTriangle + p) + degree;
        counterBottomRow = 1;

        // counter for the bottom row to be subtracted later
        for (k = 0; k < deg; k++) {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - nTermsInLayer)] *
              us[us.size(1) * iPnt + 1];
          }

          c++;
          counterBottomRow++;
        }

        deg--;
        x = ((degree + degree) + p) - 1;
        if (x < 0) {
          x = 0;
        }

        excess += x;
        d = (d + p) + 1;

        // number of terms in Pascal tetrahedron
        nTermsInPrevLayer = nTermsInLayer;
        nTermsInLayer = d + 3 * (excess - cornerTriangle);
        gap = (nTermsInPrevLayer + counterBottomRow) - 1;
        i1 = nTermsInLayer - counterBottomRow;
        for (j = 0; j <= i1; j++) {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - gap)] *
              us[us.size(1) * iPnt + 2];
          }

          c++;
        }
      }
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (order != 0) {
      //     %% compute higher order confluent Vandermonde matrix blocks incrementally
      switch (order) {
       case 1:
        {
          double scalew;
          int balance;
          int kdegree;
          int offset;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt] * hs_inv_idx_0;
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            double scaleu;
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt] * hs_inv_idx_1;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            double scalev;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt] * hs_inv_idx_2;
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (int b_i{0}; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              int degg;

              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                int partition;
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case 2:
        {
          double scaleu;
          double scalev;
          double scalew;
          double uu2;
          double vv2;
          double ww2;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt] * hs_inv_idx_0;
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt] * hs_inv_idx_1;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt] * hs_inv_idx_2;
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          uu2 = 2.0 * hs_inv_idx_0 * hs_inv_idx_0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = uu2 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uv;

            //     %% compute du*dv
            offset += us.size(0);
            uv = hs_inv_idx_0 * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = uv * V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = hs_inv_idx_1;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev += hs_inv_idx_1;
                c++;
              }

              scalev = static_cast<double>(deg) * hs_inv_idx_1;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scalev;
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree) * hs_inv_idx_0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          vv2 = 2.0 * hs_inv_idx_1 * hs_inv_idx_1;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = vv2 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uw;
            double vw;

            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            uw = hs_inv_idx_0 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = uw * V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            vw = hs_inv_idx_1 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = vw * V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          ww2 = 2.0 * hs_inv_idx_2 * hs_inv_idx_2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = ww2 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -1:
        {
          double scalew;
          int balance;
          int kdegree;
          int offset;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt] * hs_inv_idx_0;
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            double scaleu;
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt] * hs_inv_idx_1;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            double scalev;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt] * hs_inv_idx_2;
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (int b_i{0}; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              int degg;

              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                int partition;
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -2:
        {
          double scaleu;
          double scalev;
          double scalew;
          double uu2;
          double vv2;
          double ww2;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt] * hs_inv_idx_0;
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt] * hs_inv_idx_1;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt] * hs_inv_idx_2;
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          uu2 = 2.0 * hs_inv_idx_0 * hs_inv_idx_0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = uu2 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uv;

            //     %% compute du*dv
            offset += us.size(0);
            uv = hs_inv_idx_0 * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = uv * V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = hs_inv_idx_1;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev += hs_inv_idx_1;
                c++;
              }

              scalev = static_cast<double>(deg) * hs_inv_idx_1;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scalev;
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree) * hs_inv_idx_0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          vv2 = 2.0 * hs_inv_idx_1 * hs_inv_idx_1;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = vv2 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uw;
            double vw;

            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            uw = hs_inv_idx_0 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = uw * V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            vw = hs_inv_idx_1 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = vw * V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          ww2 = 2.0 * hs_inv_idx_2 * hs_inv_idx_2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = ww2 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -4:
        {
          double b_u2v2_tmp;
          double scaleu;
          double scalev;
          double scalew;
          double u2v2;
          double u2v2_tmp;
          double u2w2;
          double u2w2_tmp;
          double uu2;
          double uu4;
          double v2w2;
          double v4;
          double vv2;
          double ww2;
          double ww4;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;
          flag = (degree > 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag,
                   "Biharnomic is only supported for Pascal-tetrahedral monomials in 3D.");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr,
                    "Biharnomic is only supported for Pascal-tetrahedral monomials in 3D.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt] * hs_inv_idx_0;
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt] * hs_inv_idx_1;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt] * hs_inv_idx_2;
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          uu2 = 2.0 * hs_inv_idx_0 * hs_inv_idx_0;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = uu2 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * scaleu;
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uv;

            //     %% compute du*dv
            offset += us.size(0);
            uv = hs_inv_idx_0 * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = uv * V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = hs_inv_idx_1;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev += hs_inv_idx_1;
                c++;
              }

              scalev = static_cast<double>(deg) * hs_inv_idx_1;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scalev;
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree) * hs_inv_idx_0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          vv2 = 2.0 * hs_inv_idx_1 * hs_inv_idx_1;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = vv2 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * scalev;
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            double uw;
            double vw;

            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            uw = hs_inv_idx_0 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = uw * V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            vw = hs_inv_idx_1 * hs_inv_idx_2;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = vw * V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * scalew;
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * scalew;
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          ww2 = 2.0 * hs_inv_idx_2 * hs_inv_idx_2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = ww2 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * scalew;
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  scalew = static_cast<double>(k + 1) * hs_inv_idx_2;
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * scalew;
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^4
          offset = us.size(0) * 7;
          uu4 = 24.0 * std::pow(hs_inv_idx_0, 4.0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 20] = uu4 * V[iPnt];
          }

          c = 35;
          d = 15;
          if (degree < 0) {
            i = -degree;
          } else {
            i = degree;
          }

          for (deg = 5; deg <= i; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2)) +
                V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu -
                hs_inv_idx_0);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2)) +
                  V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu -
                  hs_inv_idx_0);
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * (c + 1)] = V[i1 + V.size(1) * ((c - d) - deg)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          //  compute du^2dv^2
          offset = us.size(0) << 3;
          u2v2_tmp = 4.0 * (hs_inv_idx_0 * hs_inv_idx_0);
          b_u2v2_tmp = hs_inv_idx_1 * hs_inv_idx_1;
          u2v2 = u2v2_tmp * b_u2v2_tmp;
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 22] = u2v2 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            scaleu = static_cast<double>(deg) * hs_inv_idx_0;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu -
                hs_inv_idx_0);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu -= hs_inv_idx_0;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                  V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu -
                  hs_inv_idx_0);
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * (c + 1)] = V[i1 + V.size(1) * ((c - d) - deg)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          //  compute dv^4
          offset = us.size(0) * 9;
          v4 = 24.0 * std::pow(hs_inv_idx_1, 4.0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 24] = v4 * V[iPnt];
          }

          c = 34;
          d = 15;
          for (deg = 5; deg <= degree; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            scalev = hs_inv_idx_1;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * (c + 2)] = V[(iPnt + 5 * stride)
                  + V.size(1) * ((c - (d << 1)) + deg)] * scalev * (scalev -
                  hs_inv_idx_1);
              }

              scalev += hs_inv_idx_1;
              c++;
            }

            scalev = static_cast<double>(deg) * hs_inv_idx_1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 2)] = V[(5 * stride + iPnt) +
                V.size(1) * ((c - (d << 1)) + deg)] * scalev * (scalev -
                hs_inv_idx_1);
            }

            c += 3;
            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (((c - d) - deg) - 1)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            d = (d + deg) + 1;
          }

          //  compute du^2*dw^2
          offset += us.size(0);
          u2w2_tmp = hs_inv_idx_2 * hs_inv_idx_2;
          u2w2 = u2v2_tmp * u2w2_tmp;
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 29] = u2w2 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2))
                    + V.size(1) * (((c - (d << 1)) - deg) - 1)] * scalew *
                    (scalew - hs_inv_idx_2);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          //  compute dv^2*dw^2
          offset += us.size(0);
          v2w2 = 4.0 * b_u2v2_tmp * u2w2_tmp;
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 31] = v2w2 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                    V.size(1) * (((c - (d << 1)) - deg) - 1)] * scalew * (scalew
                    - hs_inv_idx_2);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          //  compute dw^4
          offset += us.size(0);
          ww4 = 24.0 * std::pow(hs_inv_idx_2, 4.0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 34] = ww4 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                scalew = (static_cast<double>(k) + 1.0) * hs_inv_idx_2;
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 6 * stride) +
                    V.size(1) * (((c - (d << 1)) - deg) - 1)] * scalew * (scalew
                    - hs_inv_idx_2);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }
        }
        break;

       default:
        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(false, "Order must be 0, 1, 2, -1, -2, or -4.");

#else //MATLAB_MEX_FILE

        fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4.\n");
        fflush(stderr);
        assert(false);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        break;
      }
    }
  }

  static inline
  void gen_vander_3d(const ::coder::array<double, 2U> &us, int npoints,
    int degree, int order, const ::coder::array<double, 1U> &weights, ::coder::
    array<double, 2U> &V)
  {
    int b_degree;
    int c;
    int cornerTriangle;
    int counterBottomRow;
    int d;
    int deg;
    int excess;
    int i;
    int i1;
    int iPnt;
    int j;
    int k;
    int maxLayers;
    int nTermsInLayer;
    int nTermsInPrevLayer;
    int ncols;
    int nrblks;
    int nrows;
    int p;
    int stride;
    int x;
    boolean_T flag;

    //  Generate generalized/confluent Vandermonde matrix in 3D.
    //
    //     V = gen_vander_3d(us)
    //     V = gen_vander_3d(us, npoints)
    //     V = gen_vander_3d(us, npoints, degree, order)
    //     V = gen_vander_3d(us, npoints, degree, order, weights)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V)
    //     V = gen_vander_3d(us, npoints, degree, order, weights, hs_inv, V, stride)
    //     V_colMajor (optional): V stored in colum-major. Useful for debugging.
    //
    //  Parameters
    //  ----------
    //     us:      Local coordinates of points (n-by-3, where n>=npoints)
    //     npoints: Number of points. Use 0 for default (size(us, 1))
    //     degree:  Maximum degree of monomials (default is 2)
    //     order:   Order of derivative in confluent Vandermonde matrix. Use -1,
    //              -2, and -4 for grad, Laplacian, and biLaplacian, respectively.
    //     weights: Weights for all points (n-by-1, where n>=1; use [] or omit
    //              it to use unit weights)
    //     hs_inv:  Inverse length for scaling rows in CVM (size 1-by-0 or 1-by-d)
    //     V:       Vandermonde matrix (must be preallocated if present at input)
    //     stride:  number of rows in each row block in V (0 for default)
    //
    //  Returns
    //  -------
    //     V:      confluent Vandermonde matrix
    //
    //  Notes
    //  -----
    //     It is more efficient to provide weights here to incorporate them when
    //     constructing the Vandermonde matrix (linear-time overhead in npoints)
    //     than scaling the matrix afterwards (linear in npoints*nmonomials).
    //
    //     Entries in each row are ordered based on the levels in Pascal tetrahedron,
    //     including tensor-product monomials. The row blocks of CVM are based on
    //     the levels in Pascal tetrahedron.
    //
    //     For example, for degree=2 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1^2, u1*v1, v1^2, u1*w1, v1*w1, w1^2]
    //        weights(2) * [1, u2, v2, w2, u2^2, u2*v2, v2^2, u2*w2, v2*w2, w2^2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  2u1,  v1,    0,   w1,     0,     0]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  2u2,  v2,    0,   w2,     0,     0]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,    0,  u1,  2v1,   0,     w1,    0]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,    0,  u2,  2v2,   0,     w2,    0]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,    0,   0,    0,   u1,   v1,     2w1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,    0,   0,    0,   u2,   v2,     2w2]
    //
    //     For degree=-1 and order=1, V looks as follows:
    //        weights(1) * [1, u1, v1, w1, u1*v1, u1*w1, v1*w1, u1*v1*w1]
    //        weights(2) * [1, u2, v2, w2, u2*v2, u2*w2, v2*w2, u2*v2*w2]
    //        ...
    //        hs_inv(1)*weights(1) * [0, 1,  0,  0,  v1,   w1,   0,     v1*w1]
    //        hs_inv(1)*weights(2) * [0, 1,  0,  0,  v2,   w2,   0,     v2*w2]
    //        ...
    //        hs_inv(2)*weights(1) * [0, 0,  1,  0,  u1,   0,    w1,    u1*w1]
    //        hs_inv(2)*weights(2) * [0, 0,  1,  0,  u2,   0,    w2,    u2*w2]
    //        ...
    //        hs_inv(3)*weights(1) * [0, 0,  0,  1,   0,   u1,   v1,    u1*v1]
    //        hs_inv(3)*weights(2) * [0, 0,  0,  1,   0,   u2,   v2,    u2*v2]
    //
    //  See also gen_vander_1d, gen_vander_2d
    //  Handle input arguments
    flag = (npoints <= us.size(0));

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if ((order >= 0) || (order == -1) || (order == -2) || (order == -4)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Order must be 0, 1, 2, -1, -2, or -4");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    stride = us.size(0);
    nrblks = (order + 1) * (order + 2) * (order + 3) / 6;
    switch (order) {
     case -1:
      nrblks = 4;
      break;

     case -2:
      nrblks = 7;
      break;

     case -4:
      if (degree > 0) {
        nrblks = 13;
      } else {
        nrblks = 19;
      }
      break;
    }

    if (degree >= 0) {
      ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
    } else {
      ncols = (1 - degree) * (1 - degree) * (1 - degree);
    }

    //  Allocate storage for V
    nrows = us.size(0) * nrblks;
    if ((V.size(1) != nrows) || (V.size(0) != ncols)) {
      V.set_size(ncols, nrows);
    }

    //  compute 0th order generalized Vandermonde matrix
    //  Compute generalized Vandermonde matrix of order-0
    if (weights.size(0) == 0) {
      if (degree != 0) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
          V[iPnt + V.size(1)] = us[us.size(1) * iPnt];
          V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1];
          V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2];
        }
      } else {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt] = 1.0;
        }
      }
    } else if (degree != 0) {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
        V[iPnt + V.size(1)] = us[us.size(1) * iPnt] * weights[iPnt];
        V[iPnt + V.size(1) * 2] = us[us.size(1) * iPnt + 1] * weights[iPnt];
        V[iPnt + V.size(1) * 3] = us[us.size(1) * iPnt + 2] * weights[iPnt];
      }
    } else {
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt] = weights[iPnt];
      }
    }

    c = 4;
    d = 3;
    if (0 > order) {
      i = order;
    } else {
      i = 0;
    }

    if (-degree > 0) {
      b_degree = 1;
    } else if (-degree < 0) {
      b_degree = -1;
    } else {
      b_degree = 0;
    }

    x = b_degree * i;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    if (x < 0) {
      x = 0;
    }

    i1 = b_degree - x;
    for (deg = 2; deg <= i1; deg++) {
      //  Within each level, use convention of Pascal triangle with x^deg at peak
      for (j = 0; j < deg; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)] * us[us.size(1)
            * iPnt];
        }

        c++;
      }

      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * ((c - d) - 1)] *
          us[us.size(1) * iPnt + 1];
      }

      c++;
      for (j = 0; j < d; j++) {
        for (iPnt = 0; iPnt < npoints; iPnt++) {
          V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (((c - d) - deg) - 1)] *
            us[us.size(1) * iPnt + 2];
        }

        c++;
      }

      d = (d + deg) + 1;
    }

    //  Compute the tri-degree terms if degree<0
    if (degree < 0) {
      deg = -degree;
      maxLayers = -degree * 3 + i;

      // max number of layers needed in the Pascal tetrahedron
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      nTermsInLayer = d;

      // initializing number of elements in layer
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      i = 1 - degree;
      for (p = i; p <= maxLayers; p++) {
        int gap;

        //  Within each level, x^deg is at the peak of Pascal triangle
        //  implicitly calculating number of elements in corner Pascal triangles
        cornerTriangle = (cornerTriangle + p) + degree;
        counterBottomRow = 1;

        // counter for the bottom row to be subtracted later
        for (k = 0; k < deg; k++) {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - nTermsInLayer)] *
              us[us.size(1) * iPnt + 1];
          }

          c++;
          counterBottomRow++;
        }

        deg--;
        x = ((degree + degree) + p) - 1;
        if (x < 0) {
          x = 0;
        }

        excess += x;
        d = (d + p) + 1;

        // number of terms in Pascal tetrahedron
        nTermsInPrevLayer = nTermsInLayer;
        nTermsInLayer = d + 3 * (excess - cornerTriangle);
        gap = (nTermsInPrevLayer + counterBottomRow) - 1;
        i1 = nTermsInLayer - counterBottomRow;
        for (j = 0; j <= i1; j++) {
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[iPnt + V.size(1) * c] = V[iPnt + V.size(1) * (c - gap)] *
              us[us.size(1) * iPnt + 2];
          }

          c++;
        }
      }
    }

    flag = (order <= 2);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (order != 0) {
      //     %% compute higher order confluent Vandermonde matrix blocks incrementally
      switch (order) {
       case 1:
        {
          int balance;
          int kdegree;
          int offset;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt];
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            double scaleu;
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt];
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            double scalev;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt];
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (int b_i{0}; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * (static_cast<double>(k) + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              int degg;

              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                int partition;
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case 2:
        {
          double scaleu;
          double scalev;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt];
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt];
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt];
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * (static_cast<double>(k) + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = 2.0 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dv
            offset += us.size(0);
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = 1.0;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev++;
                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * static_cast<double>(deg);
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree);
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = 2.0 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * static_cast<
                        double>(k + 1);
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * static_cast<double>(k
                        + 1);
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = 2.0 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>(k)
                    + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -1:
        {
          int balance;
          int kdegree;
          int offset;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt];
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            double scaleu;
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt];
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            double scalev;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt];
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (int b_i{0}; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * (static_cast<double>(k) + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              int degg;

              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                int partition;
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -2:
        {
          double scaleu;
          double scalev;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt];
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt];
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt];
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * (static_cast<double>(k) + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = 2.0 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dv
            offset += us.size(0);
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = 1.0;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev++;
                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * static_cast<double>(deg);
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree);
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = 2.0 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * static_cast<
                        double>(k + 1);
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * static_cast<double>(k
                        + 1);
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = 2.0 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>(k)
                    + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }
        }
        break;

       case -4:
        {
          double scaleu;
          double scalev;
          double uu4_tmp;
          int b_i;
          int balance;
          int degg;
          int kdegree;
          int offset;
          int partition;
          flag = (degree > 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag,
                   "Biharnomic is only supported for Pascal-tetrahedral monomials in 3D.");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr,
                    "Biharnomic is only supported for Pascal-tetrahedral monomials in 3D.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          //  Compute order-1 CVM row blocks from order-0 GVM.
          flag = (degree != 0);

          //  Throw error if condition false
          //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

          mxAssert(flag, "");

#else //MATLAB_MEX_FILE

          if (!flag) {
            fprintf(stderr, "Runtime assertion error.\n");
            fflush(stderr);
          }

          assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

          // compute derivatives with respect to u
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = stride + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = V[iPnt];
            V[i + V.size(1) * 2] = 0.0;
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 3;
          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i1 = b_degree - x;
          for (deg = 2; deg <= i1; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(stride + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = stride + iPnt;
                V[b_degree + V.size(1) * (c + 1)] = V[b_degree + V.size(1) * ((c
                  - d) - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(stride + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = stride + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = stride + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to v
          offset = us.size(0) + us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            i = offset + iPnt;
            V[i] = 0.0;
            V[i + V.size(1)] = 0.0;
            V[i + V.size(1) * 2] = V[iPnt];
            V[i + V.size(1) * 3] = 0.0;
          }

          c = 4;
          d = 4;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                  * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (c - d)]
                * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i1 = offset + iPnt;
                  V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (c - nTermsInLayer)]
                    * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              i1 = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= i1; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    balance)] * us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          // compute derivatives with respect to w
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = V[iPnt];
          }

          c = 4;
          d = 3;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 1) {
            i = order + 1;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 2; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[iPnt + V.size(1) * (((c
                    - d) - deg) - 1)] * (static_cast<double>(k) + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 1) {
              i = order + 1;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[iPnt + V.size(1) * (c
                      - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^2
          offset = us.size(0) << 2;
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            V[x] = 0.0;
            V[x + V.size(1)] = 0.0;
            V[x + V.size(1) * 2] = 0.0;
            V[x + V.size(1) * 3] = 0.0;
            V[x + V.size(1) * 4] = 2.0 * V[iPnt];
            for (i = 0; i < 5; i++) {
              V[x + V.size(1) * (i + 5)] = 0.0;
            }
          }

          c = 10;
          d = 6;
          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          if (degree < 0) {
            i1 = -degree;
          } else {
            i1 = degree;
          }

          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          x = b_degree * i;
          if (x < 0) {
            x = 0;
          }

          b_degree = i1 - x;
          for (deg = 3; deg <= b_degree; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size(1)
                * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * scaleu;
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                x = offset + iPnt;
                V[x + V.size(1) * (c + 1)] = V[x + V.size(1) * ((c - d) - deg)] *
                  us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          // compute tri degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt + 1];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = (nTermsInPrevLayer + counterBottomRow) - 1;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dv
            offset += us.size(0);
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 5; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 5] = V[iPnt];
              V[x + V.size(1) * 6] = 0.0;
              V[x + V.size(1) * 7] = 0.0;
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 7;
            for (deg = 3; deg <= i1; deg++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              scalev = 1.0;
              for (j = 0; j <= deg - 2; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) +
                    V.size(1) * (c - d)] * scalev;
                }

                scalev++;
                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + stride) + V.size
                  (1) * (c - d)] * static_cast<double>(deg);
              }

              c++;
              for (kdegree = 0; kdegree <= d - 3; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  i = offset + iPnt;
                  V[i + V.size(1) * c] = V[i + V.size(1) * ((c - d) - deg)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }

              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
              d = (d + deg) + 1;
            }

            // compute the tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d - 2;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (kdegree = 0; kdegree < deg; kdegree++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt)
                      + V.size(1) * (c - nTermsInLayer)] * static_cast<double>
                      (-degree - kdegree);
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                x = (((p + degree) << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer + 1;
                nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
                balance = nTermsInPrevLayer + counterBottomRow;
                b_degree = nTermsInLayer - counterBottomRow;
                for (j = 0; j <= b_degree; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    x = offset + iPnt;
                    V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                      us[us.size(1) * iPnt + 2];
                  }

                  c++;
                }
              }
            }
          }

          //  compute dv^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 6; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 6] = 2.0 * V[iPnt];
            V[x + V.size(1) * 7] = 0.0;
            V[x + V.size(1) * 8] = 0.0;
            V[x + V.size(1) * 9] = 0.0;
          }

          c = 10;
          d = 7;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1)) +
                  V.size(1) * (c - d)] * scalev;
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[((stride << 1) + iPnt) +
                V.size(1) * (c - d)] * static_cast<double>(deg);
            }

            c++;
            for (kdegree = 0; kdegree <= d - 3; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                b_degree = offset + iPnt;
                V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * ((c - d)
                  - deg)] * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            c++;
            d = (d + deg) + 1;
          }

          // compute the tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d - 2;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  b_degree = offset + iPnt;
                  V[b_degree + V.size(1) * c] = V[b_degree + V.size(1) * (c -
                    nTermsInLayer)] * us[us.size(1) * iPnt];
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              x = (((p + degree) << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer + 1;
              nTermsInLayer = (d + 3 * (excess - cornerTriangle)) - 2;
              balance = nTermsInPrevLayer + counterBottomRow;
              b_degree = nTermsInLayer - counterBottomRow;
              for (j = 0; j <= b_degree; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  x = offset + iPnt;
                  V[x + V.size(1) * c] = V[x + V.size(1) * (c - balance)] *
                    us[us.size(1) * iPnt + 2];
                }

                c++;
              }
            }
          }

          if (order > 0) {
            //     %% compute du*dw
            offset = (offset + us.size(0)) - 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = (offset + iPnt) + 1;
              for (i = 0; i < 7; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 7] = V[iPnt];
              V[x + V.size(1) * 8] = 0.0;
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = V[(iPnt + stride)
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[((offset + iPnt) + V.size(1) * c) + 1] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[((iPnt + offset) + V.size(1) * c) + 1] = V[(iPnt +
                        stride) + V.size(1) * (c - balance)] * static_cast<
                        double>(k + 1);
                    }

                    c++;
                  }
                }
              }
            }

            //     %% compute dv*dw
            offset = (offset + us.size(0)) + 1;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              x = offset + iPnt;
              for (i = 0; i < 8; i++) {
                V[x + V.size(1) * i] = 0.0;
              }

              V[x + V.size(1) * 8] = V[iPnt];
              V[x + V.size(1) * 9] = 0.0;
            }

            c = 10;
            d = 6;
            for (deg = 3; deg <= i1; deg++) {
              for (j = 0; j <= deg; j++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
              }

              for (k = 0; k < deg; k++) {
                i = (deg - k) - 1;
                for (b_i = 0; b_i <= i; b_i++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 1))
                      + V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>
                      (k) + 1.0);
                  }

                  c++;
                }
              }

              d = (d + deg) + 1;
            }

            // compute tri-degree terms if degree < 0
            if (degree < 0) {
              deg = -degree;
              maxLayers = -degree * 3 + 1;

              // max number of layers needed in the Pascal tetrahedron
              cornerTriangle = 0;

              // number of elements subtracted in each corner Pascal triangle
              nTermsInLayer = d;

              // initializing number of elements in layer
              excess = 0;

              // excess based on overlapping of growing Pascal triangles
              i = 1 - degree;
              for (p = i; p <= maxLayers; p++) {
                //  Within each level, x^deg is at the peak of Pascal triangle
                //  implicitly calculating number of elements in corner Pascal triangles
                cornerTriangle = (cornerTriangle + p) + degree;
                counterBottomRow = 0;

                // counter for the bottom row to be subtracted later
                for (k = 0; k < deg; k++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(offset + iPnt) + V.size(1) * c] = 0.0;
                  }

                  c++;
                  counterBottomRow++;
                }

                deg--;
                b_degree = p + degree;
                x = ((b_degree << 1) - p) - 1;
                if (x < 0) {
                  x = 0;
                }

                excess += x;
                d = (d + p) + 1;

                // number of terms in Pascal tetrahedron
                nTermsInPrevLayer = nTermsInLayer;
                nTermsInLayer = d + 3 * (excess - cornerTriangle);
                balance = nTermsInPrevLayer + counterBottomRow;
                degg = -degree;
                for (k = 0; k < degg; k++) {
                  x = (b_degree - k) - 1;
                  if (x < 0) {
                    x = -x;
                  }

                  partition = -degree - x;
                  for (j = 0; j <= partition; j++) {
                    for (iPnt = 0; iPnt < npoints; iPnt++) {
                      V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + (stride <<
                        1)) + V.size(1) * (c - balance)] * static_cast<double>(k
                        + 1);
                    }

                    c++;
                  }
                }
              }
            }
          }

          //  compute dw^2
          offset += us.size(0);
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            x = offset + iPnt;
            for (i = 0; i < 9; i++) {
              V[x + V.size(1) * i] = 0.0;
            }

            V[x + V.size(1) * 9] = 2.0 * V[iPnt];
          }

          c = 10;
          d = 6;
          if (-degree > 0) {
            b_degree = 1;
          } else if (-degree < 0) {
            b_degree = -1;
          } else {
            b_degree = 0;
          }

          if (0 > order + 2) {
            i = order + 2;
          } else {
            i = 0;
          }

          x = b_degree * i;
          if (degree < 0) {
            b_degree = -degree;
          } else {
            b_degree = degree;
          }

          if (x < 0) {
            x = 0;
          }

          i = b_degree - x;
          for (deg = 3; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                    V.size(1) * (((c - d) - deg) - 1)] * (static_cast<double>(k)
                    + 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          // compute tri-degree terms if degree < 0
          if (degree < 0) {
            deg = -degree;
            if (0 > order + 2) {
              i = order + 2;
            } else {
              i = 0;
            }

            maxLayers = -degree * 3 + i;

            // max number of layers needed in the Pascal tetrahedron
            cornerTriangle = 0;

            // number of elements subtracted in each corner Pascal triangle
            nTermsInLayer = d;

            // initializing number of elements in layer
            excess = 0;

            // excess based on overlapping of growing Pascal triangles
            i = 1 - degree;
            for (p = i; p <= maxLayers; p++) {
              //  Within each level, x^deg is at the peak of Pascal triangle
              //  implicitly calculating number of elements in corner Pascal triangles
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 0;

              // counter for the bottom row to be subtracted later
              for (k = 0; k < deg; k++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = 0.0;
                }

                c++;
                counterBottomRow++;
              }

              deg--;
              b_degree = p + degree;
              x = ((b_degree << 1) - p) - 1;
              if (x < 0) {
                x = 0;
              }

              excess += x;
              d = (d + p) + 1;

              // number of terms in Pascal tetrahedron
              nTermsInPrevLayer = nTermsInLayer;
              nTermsInLayer = d + 3 * (excess - cornerTriangle);
              balance = nTermsInPrevLayer + counterBottomRow;
              degg = -degree;
              for (k = 0; k < degg; k++) {
                x = (b_degree - k) - 1;
                if (x < 0) {
                  x = -x;
                }

                partition = -degree - x;
                for (j = 0; j <= partition; j++) {
                  for (iPnt = 0; iPnt < npoints; iPnt++) {
                    V[(iPnt + offset) + V.size(1) * c] = V[(iPnt + 3 * stride) +
                      V.size(1) * (c - balance)] * static_cast<double>(k + 1);
                  }

                  c++;
                }
              }
            }
          }

          //  compute du^4
          offset = us.size(0) * 7;
          uu4_tmp = 24.0 * std::pow(1.0, 4.0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 20] = uu4_tmp * V[iPnt];
          }

          c = 35;
          d = 15;
          if (degree < 0) {
            i = -degree;
          } else {
            i = degree;
          }

          for (deg = 5; deg <= i; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2)) +
                V.size(1) * ((c - (d << 1)) + deg)] * static_cast<double>(deg) *
                (static_cast<double>(deg) - 1.0);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2)) +
                  V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu - 1.0);
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * (c + 1)] = V[i1 + V.size(1) * ((c - d) - deg)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          //  compute du^2dv^2
          offset = us.size(0) << 3;
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 22] = 4.0 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            scaleu = deg;
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                V.size(1) * ((c - (d << 1)) + deg)] * static_cast<double>(deg) *
                (static_cast<double>(deg) - 1.0);
            }

            c++;
            for (j = 0; j <= deg - 2; j++) {
              scaleu--;
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                  V.size(1) * ((c - (d << 1)) + deg)] * scaleu * (scaleu - 1.0);
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * (c + 1)] = V[i1 + V.size(1) * ((c - d) - deg)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            c += 2;
            d = (d + deg) + 1;
          }

          //  compute dv^4
          offset = us.size(0) * 9;
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 24] = uu4_tmp * V[iPnt];
          }

          c = 34;
          d = 15;
          for (deg = 5; deg <= degree; deg++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 1)] = 0.0;
            }

            scalev = 1.0;
            for (j = 0; j <= deg - 2; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * (c + 2)] = V[(iPnt + 5 * stride)
                  + V.size(1) * ((c - (d << 1)) + deg)] * scalev * (scalev - 1.0);
              }

              scalev++;
              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * (c + 2)] = V[(5 * stride + iPnt) +
                V.size(1) * ((c - (d << 1)) + deg)] * static_cast<double>(deg) *
                (static_cast<double>(deg) - 1.0);
            }

            c += 3;
            for (kdegree = 0; kdegree <= d - 2; kdegree++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                i1 = offset + iPnt;
                V[i1 + V.size(1) * c] = V[i1 + V.size(1) * (((c - d) - deg) - 1)]
                  * us[us.size(1) * iPnt + 2];
              }

              c++;
            }

            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * c] = 0.0;
            }

            d = (d + deg) + 1;
          }

          //  compute du^2*dw^2
          offset += us.size(0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 29] = 4.0 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + (stride << 2))
                    + V.size(1) * (((c - (d << 1)) - deg) - 1)] * (static_cast<
                    double>(k) + 1.0) * ((static_cast<double>(k) + 1.0) - 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          //  compute dv^2*dw^2
          offset += us.size(0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 31] = 4.0 * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 5 * stride) +
                    V.size(1) * (((c - (d << 1)) - deg) - 1)] * (static_cast<
                    double>(k) + 1.0) * ((static_cast<double>(k) + 1.0) - 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }

          //  compute dw^4
          offset += us.size(0);
          for (b_i = 0; b_i < 35; b_i++) {
            for (iPnt = 0; iPnt < npoints; iPnt++) {
              V[(offset + iPnt) + V.size(1) * b_i] = 0.0;
            }
          }

          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[(offset + iPnt) + V.size(1) * 34] = uu4_tmp * V[iPnt];
          }

          c = 35;
          d = 15;
          for (deg = 5; deg <= i; deg++) {
            for (j = 0; j <= deg; j++) {
              for (iPnt = 0; iPnt < npoints; iPnt++) {
                V[(offset + iPnt) + V.size(1) * c] = 0.0;
              }

              c++;
            }

            for (k = 0; k < deg; k++) {
              i1 = (deg - k) - 1;
              for (b_i = 0; b_i <= i1; b_i++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(offset + iPnt) + V.size(1) * c] = V[(iPnt + 6 * stride) +
                    V.size(1) * (((c - (d << 1)) - deg) - 1)] * (static_cast<
                    double>(k) + 1.0) * ((static_cast<double>(k) + 1.0) - 1.0);
                }

                c++;
              }
            }

            d = (d + deg) + 1;
          }
        }
        break;

       default:
        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(false, "Order must be 0, 1, 2, -1, -2, or -4.");

#else //MATLAB_MEX_FILE

        fprintf(stderr, "Order must be 0, 1, 2, -1, -2, or -4.\n");
        fflush(stderr);
        assert(false);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        break;
      }
    }
  }

  static inline
  void gen_vander_3d_dag(int degree, ::coder::array<unsigned char, 2U>
    &dag)
  {
    int b_i;
    int c;
    int d;
    int i;
    int j;
    int maxterms;
    int x;

    //  Build a DAG for Vandermonde matrix in 3D.
    //
    //     dag = gen_vander_3d_dag(degree)
    //     dag = gen_vander_3d_dag(degree, dag)
    //
    //  Parameters
    //  ----------
    //     degree:  Maximum degree of monomials. Use negative for tensor-product monomials.
    //
    //  Returns
    //  -------
    //     dag:     A direct acyclic graph stored in a 3-by-(#monomials+1) M-array
    //        (in column major, or a (#monomials+1)-by-3 M-array if row major).
    //        Each column represents a monomial (x^i*y^j*z^k), and its column
    //        stores the offsets to the indices of its "child" monomials
    //        (i.e., x^(i+1)*y^j*z^k, x^i*y^(j+1)*z^k, and x^i*y^j+*z^(k+1)).
    //        The last entry is used store a signature, so that it can be
    //       recomputed when needed.
    //
    //     dag_colMajor is an optional output in column-major for testing.
    //
    //  See also gen_vander_3d, gen_vander_1d_dag, gen_vander_2d_dag, rrqr_trunc
    //  Handle input arguments
    if (degree != 0) {
      dag[0] = 1U;

      //  x-child
      dag[1] = 2U;

      //  y-child
      dag[2] = 3U;

      //  z-child
    } else {
      dag[0] = 0U;
      dag[1] = 0U;
      dag[2] = 0U;

      //  No children
    }

    c = 1;
    d = 3;
    if (degree < 0) {
      i = -degree;
    } else {
      i = degree;
    }

    for (int deg{2}; deg <= i; deg++) {
      maxterms = (d + deg) + 1;
      for (j = deg; j >= 1; j--) {
        for (b_i = 0; b_i < j; b_i++) {
          x = c + b_i;
          dag[dag.size(1) * x] = static_cast<unsigned char>(d);
          dag[dag.size(1) * x + 1] = static_cast<unsigned char>(d + 1);
          dag[dag.size(1) * x + 2] = static_cast<unsigned char>(maxterms);
        }

        c += j;
        d++;
      }

      d = maxterms;
    }

    if (degree > 0) {
      i = dag.size(0);
      for (b_i = c + 1; b_i <= i; b_i++) {
        dag[dag.size(1) * (b_i - 1)] = 0U;
        dag[dag.size(1) * (b_i - 1) + 1] = 0U;
        dag[dag.size(1) * (b_i - 1) + 2] = 0U;
      }
    } else if (degree < 0) {
      int cornerTriangle;
      int excess;
      int maxlayers;
      int num_elem_group;
      maxlayers = -3 * degree + 1;
      cornerTriangle = 0;

      // number of elements subtracted in each corner Pascal triangle
      excess = 0;

      // excess based on overlapping of growing Pascal triangles
      num_elem_group = -1;
      i = -degree;
      for (int p{i}; p <= maxlayers; p++) {
        int ntermsinlayer;
        int x_tmp;
        int y;
        cornerTriangle = (cornerTriangle + p) + degree;
        y = excess << 1;
        x_tmp = p + degree;
        x = (x_tmp << 1) - p;
        if (x < 0) {
          x = 0;
        }

        excess += x;
        maxterms = ((((d + (-degree << 1)) - 3 * cornerTriangle) - p) + y) +
          excess;
        ntermsinlayer = d + 3 * (excess - cornerTriangle);
        for (int group{0}; group <= i; group++) {
          y = x_tmp - group;
          if (y < 0) {
            x = -y;
          } else {
            x = y;
          }

          num_elem_group = -degree - x;
          if (num_elem_group + 1 < 1) {
            ntermsinlayer -= 3;
          } else if (((-degree - p) + group) - 1 < 0) {
            dag[dag.size(1) * c] = 0U;
            for (b_i = 0; b_i < num_elem_group; b_i++) {
              x = c + b_i;
              dag[dag.size(1) * (x + 1)] = static_cast<unsigned char>
                (ntermsinlayer - 1);
              dag[dag.size(1) * x + 1] = static_cast<unsigned char>
                (ntermsinlayer);
              dag[dag.size(1) * x + 2] = static_cast<unsigned char>(maxterms);
            }

            c += num_elem_group;
            dag[dag.size(1) * c + 1] = 0U;
            dag[dag.size(1) * c + 2] = static_cast<unsigned char>(maxterms);
            c++;
            if (y > 1) {
              y = 1;
            }

            ntermsinlayer -= y;
          } else {
            for (b_i = 0; b_i <= num_elem_group; b_i++) {
              x = c + b_i;
              dag[dag.size(1) * x] = static_cast<unsigned char>(ntermsinlayer -
                1);
              dag[dag.size(1) * x + 1] = static_cast<unsigned char>
                (ntermsinlayer);
              dag[dag.size(1) * x + 2] = static_cast<unsigned char>(maxterms);
            }

            c = (c + num_elem_group) + 1;
            ntermsinlayer++;
          }
        }

        for (j = 0; j <= num_elem_group; j++) {
          dag[dag.size(1) * (((c + j) - num_elem_group) - 1) + 2] = 0U;
        }

        d = (d + p) + 2;
      }
    }

    //  Use last entry as signature
    i = dag.size(1) * dag.size(0) - 1;
    x = dag.size(0);
    dag[i % x * dag.size(1) + i / x] = static_cast<unsigned char>(degree + 127);
  }

  static inline
  int rrqr_factor(const ::coder::array<double, 2U> &A, double thres, int
    m, int n, ::coder::array<double, 2U> &QR, ::coder::array<int, 1U> &p, ::
    coder::array<double, 1U> &work, const ::coder::array<unsigned char, 2U> &dag)
  {
    int rank;
    boolean_T flag;

    //  Compute rank-revealing QR with column pivoting
    //
    //  [rank, QR, p] = rrqr_factor(A)
    //  [rank, QR, p] = rrqr_factor(A, thres)
    //  [rank, QR, p] = rrqr_factor(A, thres, m, n)
    //  [rank, QR, p, work] = rrqr_factor(A, thres, m, n, QR, p, work)
    //  [rank, QR, p, work] = rrqr_factor(A, thres, m, n, QR, p, work, dag)
    //
    //  Parameters
    //  ----------
    //      A:      Weighted generalized Vandermonde matrix
    //      thres:  Threshold for 2-norm condition number
    //      m:      Number of rows (use 0 for the default value of nrows(A))
    //      n:      Number of columns (use 0 for for the default value of ncolumns(A))
    //      QR:     Preallocated buffer for output QR (nrows(A)-by-(n+1))
    //      p:      Preallocated buffer for column pivoting vector (length >= n)
    //      work:   Work space; you can start with zeros(0,1) and its allocation
    //              will grow automatically as needed
    //      dag:    Direct acyclic graph for the Vandermonde matrix
    //
    //  Output
    //  ------
    //      rank:   The estimated numerical rank of matrix A.
    //      QR:     Internal representation of QRCP (including tau in (n+1)st column)
    //      p:      Column pivoting vector (1-based index) (length >= n)
    //      work:   Work space returned for reuse.
    //
    //  See also rrqr_rsolve, rrqr_rtsolve, rrqr_qmulti, rrqr_qtmulti, gen_vander_dag
    //  Must not inline to prevent buffer allocation
    //  Obtain input arguments
    if (m == 0) {
      m = A.size(1);
    }

    if (n == 0) {
      n = A.size(0);
    }

    //  Preallocate output arguments
    flag = (QR.size(1) == A.size(1));

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "The number of rows in QR must be equal to that of A.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "The number of rows in QR must be equal to that of A.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    flag = (QR.size(0) >= n + 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "The number of columns in QR must be greater than that of A.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr,
              "The number of columns in QR must be greater than that of A.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    flag = (p.size(0) >= n);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag,
             "Length of permutation vector must be no smaller than the number of columns.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr,
              "Length of permutation vector must be no smaller than the number of columns.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  Allocate work space if needed
    //  Invoke C++ function
    p[0] = 0;
    if ((0 == n) || ((dag.size(0) == 0) || (dag.size(1) == 0))) {
      rank = wls::rrqr_factor_nodag(&A[0], thres, m, n, &QR[0], &(p.data())[0],
        &(work.data())[0], work.size(0), A.size(1));
    } else {
      rank = wls::rrqr_factor(&A[0], thres, m, n, &QR[0], &(p.data())[0],
        &(work.data())[0], work.size(0), A.size(1), &dag[0], dag.size(1),
        dag.size(0));
    }

    return rank;
  }

  static inline
  void rrqr_qmulti(const ::coder::array<double, 2U> &QR, int m, int n,
    int rank, ::coder::array<double, 2U> &bs, int nrhs, ::coder::array<double,
    1U> &work)
  {
    int stride_bs;
    int u1;
    int wsize;
    boolean_T flag;

    //  Perform Q*bs, where Q is stored implicitly in QR
    //
    //     bs = rrqr_qmulti(QR, m, n, rank, bs)
    //     bs = rrqr_qmulti(QR, m, n, rank, bs, nrhs)
    //     [bs, work] = rrqr_qmulti(QR, m, n, rank, bs, nrhs, work)
    //
    //  Parameters
    //  ----------
    //     QR:      Output structure of `rrqr_factor`
    //     m:       Number of rows in Q (use 0 for default, nrows(QR))
    //     n:       Number of columns in Q (use 0 for default, ncolumns(QR)-1)
    //     rank:    Numerical rank of R (must be <= n)
    //     bs:      Right-hand side vectors of size n-by-nrhs, preallocated to
    //              max(m,n)-by-nrhs.
    //     nrhs:    Number of columns in bs (use 0 for default, ncolumns(bs))
    //     work:    Work space; you can start with zeros(1,0) and its allocation
    //              will grow automatically as needed.
    //
    //  Output
    //     bs:      Solution vectors, which overwrite bs from input
    //     work:    Work space returned for reuse.
    //
    //  See also rrqr_factor, rrqr_rsolve, rrqr_rtsolve, rrqr_qtmulti
    //  Obtain stride
    stride_bs = bs.size(1);

    //  Obtain input arguments
    if (m == 0) {
      m = QR.size(1);
    }

    if (n == 0) {
      n = QR.size(0) - 1;
    }

    if (rank == 0) {
      rank = n;
    }

    u1 = n;
    if (m <= n) {
      u1 = m;
    }

    if ((rank <= u1) && (rank >= 1)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Rank must be a positive value no greater than min(m,n).");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr,
              "Rank must be a positive value no greater than min(m,n).\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Resize work space if needed
    wsize = wls::query_work_size(m, n);
    if (work.size(0) < wsize) {
      work.set_size(wsize);
    }

    //  zero out extra rows in bs to avoid errors in LAPACK
    u1 = n + 1;
    for (int i{u1}; i <= m; i++) {
      for (int j{0}; j < nrhs; j++) {
        bs[(i + bs.size(1) * j) - 1] = 0.0;
      }
    }

    //  Invoke C++ function
    wls::rrqr_qmulti(&QR[0], m, n, rank, QR.size(1), nrhs, &bs[0], stride_bs,
                     &(work.data())[0], wsize);
  }

  static inline
  void rrqr_rtsolve(const ::coder::array<double, 2U> &QR, int n, int rank,
    ::coder::array<double, 2U> &bs, int nrhs)
  {
    int i;
    boolean_T flag;

    //  Perform forward substitution to compute bs=R'\bs, where R is stored in QR
    //
    //     bs = rrqr_rtsolve(QR, n, rank, bs)
    //     bs = rrqr_rtsolve(QR, n, rank, bs, nrhs)
    //
    //  Parameters
    //  ----------
    //    QR:      Output structure of `rrqr_factor`
    //    n:       Number of columns in R (use 0 for default, ncolumns(QR)-1)
    //    rank:    Numerical rank of R (must be >0 and <= n)
    //    bs:      Right-hand side vectors
    //    nrhs:    Number of columns in bs (use 0 for default, ncolumns(bs))
    //
    //  Output
    //    bs:      Solution vectors, which overwrite bs from input
    //
    //  See also rrqr_factor, rrqr_rsolve, rrqr_qmulti, rrqr_qtmulti
    //  Obtain input arguments
    if (n == 0) {
      n = QR.size(0) - 1;
    }

    if (rank == 0) {
      rank = n;
    }

    if (QR.size(1) > n) {
      i = n;
    } else {
      i = QR.size(1);
    }

    if ((rank <= i) && (rank >= 1)) {
      flag = true;
    } else {
      flag = false;
    }

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "Rank must be a positive value no greater than n.");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Rank must be a positive value no greater than n.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Obtain stride
    //  Invoke C++ function
    wls::rrqr_rtsolve(&QR[0], n, rank, QR.size(1), nrhs, &bs[0], bs.size(1));
  }

  static inline
  void wls_buhmann_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const ::coder::array<double, 1U> &params_sh, const ::
    coder::array<double, 2U> &params_pw, ::coder::array<double, 1U> &ws)
  {
    static const double b_dv[7]{ 2.6, 2.0, 1.6, 1.6, 1.6, 1.5, 1.4 };

    double dist_k;
    double rho;
    double sigma;
    int abs_degree;
    int i;
    boolean_T flag;

    //  Weights based on Buhmann's radial basis function
    //
    //     ws = wls_buhmann_weights(us, npoints, degree, params_sh, params_pw, ws)
    //
    //     The weight for the ith point is computed as
    //       w_i = gamma_i *
    //       (1/9+r_i*r_i*(-14/15+r_i*sqrt(r_i)*(16/3+sqrt(r_i)*(-7+sqrt(r_i)*112/45))))
    //     where r_i = norm(us(i,:), 2)/rho   if degree>=0, and
    //           r_i = norm(us(i,:), inf)/rho if degree <0,
    //     for rho=sigma*dist_k and gamma_i>=0.
    //
    //  Parameters
    //  ----------
    //     us:  Shifted and rescaled local coordinates in d dimensions (size m-by-d)
    //     npoints:  Number of points (use 0 for the default value size(us, 1))
    //     degree:   Degree of polynomials (default is 2)
    //     params_sh: Shared parameters of size k-by-1 (use empty for default).
    //           params_sh(1): sigma (use 0 for the degree-based default)
    //     params_pw: Pointwise parameters of size npoints-by-1 (use empty for default).
    //           params_pw(:,1): gamma value for each point for safeguards
    //     ws:       Preallocated storage for the weights
    //
    //  Returns
    //  -------
    //     ws:       Computed weights for each point.
    //
    //  See also
    //     wls_init, wls_invdist_weights, wls_eno_weights, WlsWeight
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "npoints must be nonnegative");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (degree == 0) {
      degree = 2;
    }

    if (degree < 0) {
      abs_degree = 1 - degree;
    } else {
      abs_degree = degree + 1;
    }

    if ((params_sh.size(0) != 0) && (params_sh[0] != 0.0)) {
      sigma = params_sh[0];

      //  Assign default rho
    } else if (abs_degree - 1 >= 7) {
      sigma = 1.4;
    } else {
      sigma = b_dv[abs_degree - 2];
    }

    if (ws.size(0) == 0) {
      ws.set_size(npoints);
    } else {
      flag = (ws.size(0) >= npoints);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "length of ws cannot be smaller than npoints");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "length of ws cannot be smaller than npoints\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    //  Compute rho to be sigma times the kth distance for k=ceil(1.5*ncoff)
    compute_distances(us, npoints, degree < 0, ws);
    if (us.size(1) == 1) {
      i = abs_degree;
    } else if (us.size(1) == 2) {
      if (degree < 0) {
        i = abs_degree * abs_degree;
      } else {
        i = (abs_degree + 1) * abs_degree / 2;
      }
    } else if (degree < 0) {
      i = abs_degree * abs_degree * abs_degree;
    } else {
      i = (abs_degree + 2) * (abs_degree + 1) * abs_degree / 6;
    }

    dist_k = find_kth_shortest_dist(ws, (i * 3 + 1) / 2, 1, npoints);
    rho = sigma * dist_k;
    if ((params_pw.size(0) == 0) || (params_pw.size(1) == 0)) {
      for (int b_i{0}; b_i < npoints; b_i++) {
        double r;
        if (degree > 0) {
          double d;
          double r2;

          //  Compute 2-norm
          d = us[us.size(1) * b_i];
          r2 = d * d;
          i = us.size(1);
          for (int j{2}; j <= i; j++) {
            d = us[(j + us.size(1) * b_i) - 1];
            r2 += d * d;
          }

          r = std::sqrt(r2);
        } else {
          //  Compute inf-norm for tensor-product
          r = std::abs(us[us.size(1) * b_i]);
          i = us.size(1);
          for (int j{2}; j <= i; j++) {
            double r1;
            r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
            if (r1 > r) {
              r = r1;
            }
          }
        }

        if (r > rho) {
          ws[b_i] = 0.0;
        } else {
          double r_sqrt;
          r /= rho;
          r_sqrt = std::sqrt(r);
          ws[b_i] = r * r * (r * r_sqrt * (r_sqrt * (r_sqrt * 112.0 / 45.0 +
            -7.0) + 5.333333333333333) + -0.93333333333333335) +
            0.1111111111111111;
        }
      }
    } else {
      flag = (params_pw.size(0) >= npoints);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "size(params_pw,1) should be >=npoints");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "size(params_pw,1) should be >=npoints\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      for (int b_i{0}; b_i < npoints; b_i++) {
        double b_gamma;
        b_gamma = params_pw[params_pw.size(1) * b_i];
        if (b_gamma <= 0.0) {
          ws[b_i] = 0.0;
        } else {
          double r;
          if (degree > 0) {
            double d;
            double r2;

            //  Compute 2-norm
            d = us[us.size(1) * b_i];
            r2 = d * d;
            i = us.size(1);
            for (int j{2}; j <= i; j++) {
              d = us[(j + us.size(1) * b_i) - 1];
              r2 += d * d;
            }

            r = std::sqrt(r2);
          } else {
            //  Compute inf-norm for tensor-product
            r = std::abs(us[us.size(1) * b_i]);
            i = us.size(1);
            for (int j{2}; j <= i; j++) {
              double r1;
              r1 = std::abs(us[(j + us.size(1) * b_i) - 1]);
              if (r1 > r) {
                r = r1;
              }
            }
          }

          if (r > rho) {
            ws[b_i] = 0.0;
          } else {
            double r_sqrt;
            r /= rho;
            r_sqrt = std::sqrt(r);
            ws[b_i] = b_gamma * (r * r * (r * r_sqrt * (r_sqrt * (r_sqrt * 112.0
              / 45.0 + -7.0) + 5.333333333333333) + -0.93333333333333335) +
                                 0.1111111111111111);
          }
        }
      }
    }
  }

  static inline
  void wls_eno_weights(const ::coder::array<double, 2U> &us, int npoints,
    int degree, const ::coder::array<double, 2U> &us_unscaled, const ::coder::
    array<double, 1U> &params_sh, const ::coder::array<double, 2U> &params_pw, ::
    coder::array<double, 1U> &ws)
  {
    double alpha;
    double b_degree;
    double c0;
    double c1;
    double c1dfg;
    double epsilon;
    double epsilon_ENO;
    double h2bar;
    double h2bar_tmp;
    double safegauard;
    int i;
    boolean_T flag;

    //  WLS-ENO weights based on function values
    //
    //     ws = wls_eno_weights(us, npoints, degree, us_unscaled, params_sh, params_pw, ws)
    //
    //     The weight for the ith point is computed as
    //       w_i = gamma*(r_i^2 + epsilon_ID) .^ -0.25 / ...
    //           (c0*|f_i-g0|^2 + c1*df*alpha_i + epsilon_ENO*df^2*hbar^2)
    //     where r_i = norm(us(i,:), 2)   if degree>=0 and
    //           r_i = norm(us(i,:), inf) if degree <0,
    //
    //  Parameters
    //  ----------
    //     us:  Shifted and rescaled local coordinates in d dimensions (size m-by-d)
    //     npoints:  Number of points (use 0 for the default value size(us, 1))
    //     degree:   Degree of polynomials (negative for tensor-product)
    //     us_unscaled:    Unscaled local coordinates for computing hbar_squared
    //     params_sh: Shared parameters of size k-by-1, where k>=3
    //           params_sh(1): g0, reference solution at center node
    //           params_sh(2): df_global, i.e., global range of f
    //           params_sh(3): c0 (use 0 for default 1)
    //           params_sh(4): c1 (use 0 for default 0.05)
    //           params_sh(5): epsilon_ID (use 0 for default 0.01)
    //           params_sh(6): epsilon_ENO (use 0 for default 0.001)
    //     params_pw: Pointwise parameters of size npoints-by-2 or npoints-by-3
    //           params_pw(:,1): function value for each point
    //           params_pw(:,2): max |alpha| in 1-ring for each node
    //           params_pw(:,3): gamma value for each point (optional)
    //     ws:       Preallocated storage for the weights
    //
    //  Returns
    //  -------
    //     ws:       Computed weight for each point
    //
    //  See also
    //     wls_init, wls_invdist_weights, wls_buhmann_weights, WlsWeight
    flag = (params_sh.size(0) >= 2);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "first two shared parameters are required");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "first two shared parameters are required\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    flag = (params_pw.size(0) >= npoints);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "size(params_pw,1) should be >=npoints");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "size(params_pw,1) should be >=npoints\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    flag = (params_pw.size(1) >= 2);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "size(params_pw,2) should be >=2");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "size(params_pw,2) should be >=2\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "npoints must be nonnegative");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    if (ws.size(0) == 0) {
      ws.set_size(npoints);
    } else {
      flag = (ws.size(0) >= npoints);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "length of ws cannot be smaller than npoints");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "length of ws cannot be smaller than npoints\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    //  Compute hbar using ws as buffer space
    compute_distances(us_unscaled, npoints, degree < 0, ws);
    h2bar = ws[0] * ws[0];
    for (i = 2; i <= npoints; i++) {
      h2bar_tmp = ws[i - 1];
      h2bar += h2bar_tmp * h2bar_tmp;
    }

    h2bar /= static_cast<double>(npoints);

    //  Evaluate the inverse-distance weights as base
    b_degree = 0.5 - static_cast<double>(degree < 0);

    //  use 0.5 or -0.5
    if ((params_sh.size(0) >= 5) && (params_sh[4] != 0.0)) {
      h2bar_tmp = params_sh[4];
    } else {
      h2bar_tmp = 0.01;
    }

    //  Weights based on inverse distance
    //
    //     ws = wls_invdist_weights(us, npoints, degree, params_sh, params_pw, ws)
    //
    //     The weight for the ith point is computed as
    //       w_i = gamma_i * (r_i + epsilon) .^ -alpha
    //     where r_i = norm(us(i,:), 2)   if degree>=0 and
    //           r_i = norm(us(i,:), inf) if degree <0,
    //     alpha=degree/2 by default, and gamma_i>=0;
    //
    //  Parameters
    //  ----------
    //     us:  Shifted and rescaled local coordinates in d dimensions (size m-by-d)
    //     npoints:  Number of points (use 0 for the default value size(us, 1))
    //     degree:   Degree of polynomials (negative for tensor-product)
    //     params_sh: Shared parameters of size k-by-1, where k>=0.
    //           params_sh(1): epsilon (use 0 for default 0.01)
    //           params_sh(2): alpha (use 0 for default abs(degree)/2)
    //     params_pw: Pointwise parameters of size npoints-by-1 (use empty for default).
    //           params_pw(:,1): gamma value for each point for safeguards
    //     ws:       Preallocated storage for the weights
    //
    //  Returns
    //  -------
    //     ws:       Computed weight for each point
    //
    //  See also
    //     wls_init, wls_buhmann_weights, wls_eno_weights, WlsWeight
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "npoints must be nonnegative");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    epsilon = 0.01;
    alpha = std::abs(b_degree) / 2.0;
    if (h2bar_tmp != 0.0) {
      epsilon = h2bar_tmp;
    }

    flag = (ws.size(0) >= npoints);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "length of ws cannot be smaller than npoints");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "length of ws cannot be smaller than npoints\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    for (i = 0; i < npoints; i++) {
      double r;
      double r2;
      r = std::abs(us[us.size(1) * i]);
      if (us.size(1) > 1) {
        if (b_degree > 0.0) {
          int b_i;

          //  Compute 2-norm
          r2 = r * r;
          b_i = us.size(1);
          for (int c_i{2}; c_i <= b_i; c_i++) {
            h2bar_tmp = us[(c_i + us.size(1) * i) - 1];
            r2 += h2bar_tmp * h2bar_tmp;
          }
        } else {
          int b_i;

          //  Compute inf-norm for tensor-product
          b_i = us.size(1);
          for (int c_i{2}; c_i <= b_i; c_i++) {
            double r1;
            r1 = std::abs(us[(c_i + us.size(1) * i) - 1]);
            if (r1 > r) {
              r = r1;
            }
          }

          r2 = r * r;
        }
      } else {
        r2 = r * r;
      }

      //  Compute weight
      ws[i] = std::pow(r2 + epsilon, -alpha);
    }

    if ((params_sh.size(0) >= 3) && (params_sh[2] != 0.0)) {
      c0 = params_sh[2];
    } else {
      c0 = 1.0;
    }

    if ((params_sh.size(0) >= 4) && (params_sh[3] != 0.0)) {
      c1 = params_sh[3];
    } else {
      c1 = 0.05;
    }

    c1dfg = c1 * params_sh[1];
    if ((params_sh.size(0) >= 6) && (params_sh[5] != 0.0)) {
      epsilon_ENO = params_sh[5];
    } else {
      epsilon_ENO = 0.001;
    }

    safegauard = epsilon_ENO * (params_sh[1] * params_sh[1]) * h2bar;
    if (params_pw.size(1) > 2) {
      for (i = 0; i < npoints; i++) {
        h2bar_tmp = params_pw[params_pw.size(1) * i] - params_sh[0];
        ws[i] = ws[i] / ((c0 * (h2bar_tmp * h2bar_tmp) + c1dfg *
                          params_pw[params_pw.size(1) * i + 1]) + safegauard);
      }
    }
  }

  static inline
  void wls_invdist_weights(const ::coder::array<double, 2U> &us, int
    npoints, int degree, const ::coder::array<double, 1U> &params_sh, const ::
    coder::array<double, 2U> &params_pw, ::coder::array<double, 1U> &ws)
  {
    double alpha;
    double epsilon;
    int b_degree;
    boolean_T flag;

    //  Weights based on inverse distance
    //
    //     ws = wls_invdist_weights(us, npoints, degree, params_sh, params_pw, ws)
    //
    //     The weight for the ith point is computed as
    //       w_i = gamma_i * (r_i + epsilon) .^ -alpha
    //     where r_i = norm(us(i,:), 2)   if degree>=0 and
    //           r_i = norm(us(i,:), inf) if degree <0,
    //     alpha=degree/2 by default, and gamma_i>=0;
    //
    //  Parameters
    //  ----------
    //     us:  Shifted and rescaled local coordinates in d dimensions (size m-by-d)
    //     npoints:  Number of points (use 0 for the default value size(us, 1))
    //     degree:   Degree of polynomials (negative for tensor-product)
    //     params_sh: Shared parameters of size k-by-1, where k>=0.
    //           params_sh(1): epsilon (use 0 for default 0.01)
    //           params_sh(2): alpha (use 0 for default abs(degree)/2)
    //     params_pw: Pointwise parameters of size npoints-by-1 (use empty for default).
    //           params_pw(:,1): gamma value for each point for safeguards
    //     ws:       Preallocated storage for the weights
    //
    //  Returns
    //  -------
    //     ws:       Computed weight for each point
    //
    //  See also
    //     wls_init, wls_buhmann_weights, wls_eno_weights, WlsWeight
    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "npoints must be nonnegative");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    epsilon = 0.01;
    if (degree < 0) {
      b_degree = -degree;
    } else {
      b_degree = degree;
    }

    alpha = static_cast<double>(b_degree) / 2.0;
    if ((params_sh.size(0) != 0) && (params_sh[0] != 0.0)) {
      epsilon = params_sh[0];
    }

    if ((params_sh.size(0) > 1) && (params_sh[1] != 0.0)) {
      alpha = params_sh[1];
    }

    if (ws.size(0) == 0) {
      ws.set_size(npoints);
    } else {
      flag = (ws.size(0) >= npoints);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "length of ws cannot be smaller than npoints");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "length of ws cannot be smaller than npoints\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    if ((params_pw.size(0) == 0) || (params_pw.size(1) == 0)) {
      for (int i{0}; i < npoints; i++) {
        double r;
        double r2;
        r = std::abs(us[us.size(1) * i]);
        if (us.size(1) > 1) {
          if (degree > 0) {
            //  Compute 2-norm
            r2 = r * r;
            b_degree = us.size(1);
            for (int b_i{2}; b_i <= b_degree; b_i++) {
              double d;
              d = us[(b_i + us.size(1) * i) - 1];
              r2 += d * d;
            }
          } else {
            //  Compute inf-norm for tensor-product
            b_degree = us.size(1);
            for (int b_i{2}; b_i <= b_degree; b_i++) {
              double r1;
              r1 = std::abs(us[(b_i + us.size(1) * i) - 1]);
              if (r1 > r) {
                r = r1;
              }
            }

            r2 = r * r;
          }
        } else {
          r2 = r * r;
        }

        //  Compute weight
        ws[i] = std::pow(r2 + epsilon, -alpha);
      }
    } else {
      flag = (params_pw.size(0) >= npoints);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "size(params_pw,1) should be >=npoints");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "size(params_pw,1) should be >=npoints\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      for (int i{0}; i < npoints; i++) {
        double b_gamma;
        b_gamma = params_pw[params_pw.size(1) * i];
        if (b_gamma <= 0.0) {
          ws[i] = 0.0;
        } else {
          double r;
          double r2;
          r = std::abs(us[us.size(1) * i]);
          if (us.size(1) > 1) {
            if (degree > 0) {
              //  Compute 2-norm
              r2 = r * r;
              b_degree = us.size(1);
              for (int b_i{2}; b_i <= b_degree; b_i++) {
                double d;
                d = us[(b_i + us.size(1) * i) - 1];
                r2 += d * d;
              }
            } else {
              //  Compute inf-norm for tensor-product
              b_degree = us.size(1);
              for (int b_i{2}; b_i <= b_degree; b_i++) {
                double r1;
                r1 = std::abs(us[(b_i + us.size(1) * i) - 1]);
                if (r1 > r) {
                  r = r1;
                }
              }

              r2 = r * r;
            }
          } else {
            r2 = r * r;
          }

          //  Compute weight
          ws[i] = b_gamma * std::pow(r2 + epsilon, -alpha);
        }
      }
    }
  }

  static inline
  void wls_resize(WlsObject *b_wls, int dim, int npoints, int degree, int
    order, boolean_T use_dag)
  {
    int b_degree;
    int ncols;
    int nrows;
    int stride;
    int u1;

    //  Reinitialize the buffers of WlsObject
    //
    //     wls = wls_resize(wls, dim, npoints)
    //     wls = wls_resize(wls, dim, npoints, degree)
    //     wls = wls_resize(wls, dim, npoints, degree, order)
    //     wls = wls_resize(wls, dim, npoints, degree, order, use_dag)
    //
    //  Parameters
    //  ----------
    //     wls:       Previously allocated WlsObject object
    //     dim:       Dimension
    //     npoints:   Number of points.
    //     degree:    Degree of polynomials (negative for tensor product monomials).
    //     order:     Order of the confluent Vandermonde system (default is 0).
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //
    //  Returns
    //  -------
    //     wls:       a MATLAB/C struct that contains the buffer spaces for
    //                computing WLS and its differentiation.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, or use_dag is missing, then its
    //  corresponding value in the input wls object will be preserved. If
    //  present, they will overwrite the values in `wls`.
    //
    //   See also WlsObject, wls_init, wls_func, wls_var_diff
    b_wls->degree = degree;
    b_wls->order = order;
    b_wls->use_dag = use_dag;

    //  make stride a multiple of four
    stride = ((npoints + 3) / 4) << 2;
    b_wls->stride = stride;
    if (order == 0) {
      nrows = stride;
    } else if (order == 1) {
      nrows = (dim + 1) * stride;
    } else {
      int nrblks;
      switch (dim) {
       case 1:
        nrblks = order + 1;
        break;

       case 2:
        nrblks = (order + 1) * (order + 2) / 2;
        break;

       default:
        //  Throw error if condition false
        //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

        mxAssert(true, "Dimension must be 1, 2, or 3.");

#else //MATLAB_MEX_FILE

        assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

        nrblks = (order + 1) * (order + 2) * (order + 3) / 6;
        break;
      }

      nrows = nrblks * stride;
    }

    b_wls->us.set_size(stride, dim);
    b_wls->rweights.set_size(stride);
    b_wls->npoints = npoints;

    //  determine number of columns and allocate V and QR
    switch (dim) {
     case 1:
      if (degree < 0) {
        b_degree = -degree;
      } else {
        b_degree = degree;
      }

      ncols = b_degree + 1;
      break;

     case 2:
      if (degree > 0) {
        ncols = (degree + 1) * (degree + 2) / 2;
      } else {
        ncols = (1 - degree) * (1 - degree);
      }
      break;

     default:
      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "Dimension must be 1, 2, or 3.");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      if (degree > 0) {
        ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;
      } else {
        ncols = (1 - degree) * (1 - degree) * (1 - degree);
      }
      break;
    }

    b_wls->hs_inv.size[1] = dim;
    b_wls->hs_inv.size[0] = 1;
    if (use_dag) {
      if ((dim != b_wls->dag.size(1)) || (ncols + 1 != b_wls->dag.size(0))) {
        //  Reset DAG if dimension or degree has changed.
        //  Add one additional row for storing the signature.
        b_wls->dag.set_size(ncols + 1, dim);
        u1 = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_degree = b_wls->dag.size(0);
        b_wls->dag[u1 % b_degree * b_wls->dag.size(1) + u1 / b_degree] =
          MAX_uint8_T;
      }
    } else {
      b_wls->dag.set_size(0, dim);
    }

    b_wls->jpvt.set_size(ncols);
    b_wls->V.set_size(ncols, nrows);
    b_wls->QR.set_size(ncols + 1, nrows);
    b_wls->rank = 0;

    //  work space
    b_degree = ncols << 2;
    u1 = ncols + 1;
    if (nrows >= u1) {
      u1 = nrows;
    }

    if (b_degree < 4160) {
      b_degree = 4160;
    }

    b_wls->work.set_size((u1 << 5) + b_degree);
  }

  static inline
  void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, const
                 ::coder::array<double, 2U> &fs, int npoints, ::coder::array<
                 double, 2U> &vdops, ::coder::array<double, 2U> &result)
  {
    int iPoint;
    int iRow;
    int j;
    int nDims;
    int ncols;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    //
    //  [wls, vdops] = wls_func(wls, pnts) computes fitting at given points.
    //  [wls, vdops, result] = wls_func(wls, pnts, fs)
    //  [wls, vdops, result] = wls_func(wls, pnts, fs, npoints) specifies number
    //    of points in pnts
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     pnts:      Local coordinates of the quadrature points (n-by-d)
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //     npoints:   Number of points in pnts
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-n. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size n-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_func
    ncols = b_wls->ncols;
    nDims = pnts.size(1);

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((npoints + 3) / 4) << 2, pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint < npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
          iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the generalized Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, npoints, b_wls->degree, 0, b_wls->hs_inv.data,
               b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(npoints, nrows_vdops);

    //  Extract vopts from Vandermonde matrix
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint < npoints; iPoint++) {
        b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iPoint] = b_wls->
          V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, npoints);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                npoints, b_wls->work);

    //  Transpose the operator for row-major
    vdops.set_size(nrows_vdops, npoints);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j < npoints; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k < npoints; k++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(npoints, fs.size(1));
      u0 = fs.size(1) * npoints;
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        for (int iPnt{0}; iPnt < npoints; iPnt++) {
          for (iRow = 0; iRow <= nrows; iRow++) {
            result[iFunc + result.size(1) * iPnt] = result[iFunc + result.size(1)
              * iPnt] + fs[iFunc + fs.size(1) * iRow] * vdops[iPnt + vdops.size
              (1) * iRow];
          }
        }
      }
    }
  }

  static inline
  void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, const
                 ::coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                 &vdops, ::coder::array<double, 2U> &result)
  {
    int iPoint;
    int iRow;
    int j;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    //
    //  [wls, vdops] = wls_func(wls, pnts) computes fitting at given points.
    //  [wls, vdops, result] = wls_func(wls, pnts, fs)
    //  [wls, vdops, result] = wls_func(wls, pnts, fs, npoints) specifies number
    //    of points in pnts
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     pnts:      Local coordinates of the quadrature points (n-by-d)
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //     npoints:   Number of points in pnts
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-n. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size n-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_func
    npoints = pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = pnts.size(1);

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((pnts.size(0) + 3) / 4) << 2, pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
          iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the generalized Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, pnts.size(0), b_wls->degree, 0, b_wls->hs_inv.data,
               b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(pnts.size(0), nrows_vdops);

    //  Extract vopts from Vandermonde matrix
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iPoint] = b_wls->
          V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, pnts.size(0));
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                pnts.size(0), b_wls->work);

    //  Transpose the operator for row-major
    vdops.set_size(nrows_vdops, pnts.size(0));
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= npoints; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= npoints; k++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(pnts.size(0), fs.size(1));
      u0 = fs.size(1) * pnts.size(0);
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        for (int iPnt{0}; iPnt <= npoints; iPnt++) {
          for (iRow = 0; iRow <= nrows; iRow++) {
            result[iFunc + result.size(1) * iPnt] = result[iFunc + result.size(1)
              * iPnt] + fs[iFunc + fs.size(1) * iRow] * vdops[iPnt + vdops.size
              (1) * iRow];
          }
        }
      }
    }
  }

  static inline
  void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, ::
                 coder::array<double, 2U> &vdops)
  {
    int iPoint;
    int j;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    //
    //  [wls, vdops] = wls_func(wls, pnts) computes fitting at given points.
    //  [wls, vdops, result] = wls_func(wls, pnts, fs)
    //  [wls, vdops, result] = wls_func(wls, pnts, fs, npoints) specifies number
    //    of points in pnts
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     pnts:      Local coordinates of the quadrature points (n-by-d)
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //     npoints:   Number of points in pnts
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-n. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size n-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_func
    npoints = pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = pnts.size(1);

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((pnts.size(0) + 3) / 4) << 2, pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
          iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the generalized Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, pnts.size(0), b_wls->degree, 0, b_wls->hs_inv.data,
               b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(pnts.size(0), nrows_vdops);

    //  Extract vopts from Vandermonde matrix
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iPoint] = b_wls->
          V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, pnts.size(0));
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                pnts.size(0), b_wls->work);

    //  Transpose the operator for row-major
    vdops.set_size(nrows_vdops, pnts.size(0));
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= npoints; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= npoints; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, boolean_T interp0,
                 boolean_T use_dag, int npoints)
  {
    int dim;
    boolean_T flag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    b_wls->interp0 = interp0;
    b_wls->use_dag = use_dag;
    if (npoints <= 0) {
      npoints = us.size(0);
    } else {
      flag = (npoints <= us.size(0));

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag,
               "Number of points cannot be greater than the first dimension of `us`.");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr,
                "Number of points cannot be greater than the first dimension of `us`.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    }

    //  Resize buffers
    wls_resize(b_wls, us.size(1), npoints, degree, order, use_dag);
    if (npoints != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i < npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i < npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i < npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i < npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i < npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i < npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, npoints, degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, npoints, degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, npoints, degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, npoints, degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (npoints != b_wls->stride) && (npoints != b_wls->stride))
      {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = npoints;
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i < npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * npoints;
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j < npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j < npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j < npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (npoints != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = npoints;
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i < npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg += npoints;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us)
  {
    int degree;
    int dim;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;

    //  Default is to use unit weight
    degree = b_wls->degree;
    order = b_wls->order;
    use_dag = b_wls->use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), b_wls->degree, b_wls->order,
               b_wls->use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int i1;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            i1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > i1) {
                wls_idx_0 -= i1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            i1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > i1) {
                wls_idx_0 -= i1;
              }
            }
          }

          i1 = b_wls->us.size(0);
          b_wls->us[b_i % i1 * b_wls->us.size(1) + b_i / i1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (order == 0) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));

        //  unit weights
        wls_idx_0 = b_wls->rweights.size(0);
        b_wls->rweights.set_size(wls_idx_0);
        for (i = 0; i < wls_idx_0; i++) {
          b_wls->rweights[i] = 1.0;
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight)
  {
    int degree;
    int dim;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    degree = b_wls->degree;
    order = b_wls->order;
    use_dag = b_wls->use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), b_wls->degree, b_wls->order,
               b_wls->use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            boolean_T flag;
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, us.size(0), degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree)
  {
    int dim;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    order = b_wls->order;
    use_dag = b_wls->use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), degree, b_wls->order,
               b_wls->use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            boolean_T flag;
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, us.size(0), degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order)
  {
    int dim;
    int npoints;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    use_dag = b_wls->use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), degree, order, b_wls->use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            boolean_T flag;
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, us.size(0), degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, boolean_T interp0)
  {
    int dim;
    int npoints;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    b_wls->interp0 = interp0;
    use_dag = b_wls->use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), degree, order, b_wls->use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            boolean_T flag;
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, us.size(0), degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, boolean_T interp0,
                 boolean_T use_dag)
  {
    int dim;
    int npoints;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    //
    //     wls = wls_init(wls, us)
    //     wls = wls_init(wls, us, weight)
    //     wls = wls_init(wls, us, weight, degree)
    //     wls = wls_init(wls, us, weight, degree, order)
    //     wls = wls_init(wls, us, weight, degree, order, interp0)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag)
    //     wls = wls_init(wls, us, weight, degree, order, interp0, use_dag, npoints)
    //
    //  Parameters
    //  ----------
    //     wls:     An WlsObject with (partially) preallocated work space.
    //     us:      input points in local coordinate system (m-by-dim)
    //     weight:  An WlsWeight object specifying the weighting scheme.
    //     degree:  degree of polynomials (default is 2)
    //     order:   degree of derivatives in the input (default is 0)
    //     interp0: Whether to enforce WLS to pass through the first point.
    //     use_dag:   Whether to use DAG if truncation is needed for QRCP,
    //                so that high-degree polynomials will be truncated first.
    //                (Default is true. Use false for better efficiency.)
    //     npoints: number of points (default is size(us,1))
    //
    //  Returns
    //  -------
    //     wls:   The WlsObject containing the V matrix, its QR factorization,
    //            permutation vector, rank, weights, dag, and work space.
    //
    //  Notes
    //  -----
    //  If one of the arguments degree, order, interp0, or use_dag is missing,
    //  then its corresponding value in the input wls object will be preserved.
    //
    //  See also
    //     WlsWeight, WlsObject, wls_var_kernel, wls_var_diff, wls_var_func, wls_var_grad,
    //  Process input arguments
    dim = us.size(1) - 1;
    b_wls->interp0 = interp0;
    b_wls->use_dag = use_dag;
    npoints = us.size(0) - 1;

    //  Resize buffers
    wls_resize(b_wls, us.size(1), us.size(0), degree, order, use_dag);
    if (us.size(0) != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;
      int wls_idx_0;
      boolean_T b;
      boolean_T b1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        wls_idx_0 = b_wls->dag.size(0);
        if (b_wls->dag[i % wls_idx_0 * b_wls->dag.size(1) + i / wls_idx_0] !=
            degree + 127) {
          //  Wrapper function for building DAG in nD.
          //
          //     dag = gen_vander_dag(dim, degree)
          //     dag = gen_vander_dag(dim, degree, dag)
          //
          //  Parameters
          //  ----------
          //     dim:     Dimension
          //     degree:  Maximum degree of monomials. Use negative
          //              for tensor-product monomials.
          //
          //  Returns
          //  -------
          //     dag: Dependency graph of Vandermonde system
          //
          //  See also gen_vander, rrqr_trunc
          //  Compute the dag
          switch (us.size(1)) {
           case 1:
            gen_vander_1d_dag(degree, b_wls->dag);
            break;

           case 2:
            gen_vander_2d_dag(degree, b_wls->dag);
            break;

           default:
            gen_vander_3d_dag(degree, b_wls->dag);
            break;
          }
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          maxx = std::fmax(maxx, std::abs(us[wls_idx_0]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(us[us.size(1) * b_i]), std::
            abs(us[us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(us[us.size(1) *
            b_i]), std::abs(us[us.size(1) * b_i + 1])), std::abs(us[us.size(1) *
            b_i + 2])));
        }
        break;
      }

      if (maxx == 0.0) {
        maxx_inv = 1.0;
      } else {
        maxx_inv = 1.0 / maxx;
      }

      for (b_i = 0; b_i <= dim; b_i++) {
        b_wls->hs_inv.data[b_i] = maxx_inv;
      }

      //  scale us and save solution into wls.us
      switch (us.size(1)) {
       case 1:
        b = true;
        b1 = (us.size(1) <= 0);
        i = us.size(1) * us.size(0);
        wls_idx_0 = 0;
        for (b_i = 0; b_i <= npoints; b_i++) {
          if (b1 || (b_i >= i)) {
            wls_idx_0 = 0;
            b = true;
          } else if (b) {
            b = false;
            wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
          } else {
            u1 = us.size(1) * us.size(0) - 1;
            if (wls_idx_0 > MAX_int32_T - us.size(1)) {
              wls_idx_0 = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
            } else {
              wls_idx_0 += us.size(1);
              if (wls_idx_0 > u1) {
                wls_idx_0 -= u1;
              }
            }
          }

          u1 = b_wls->us.size(0);
          b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[wls_idx_0] *
            maxx_inv;
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] * maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] *
            maxx_inv;
          b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] *
            maxx_inv;
        }
        break;
      }

      //  Compute point-wise weights
      if (((weight->name.size(1) == 0) || (weight->name[0] == 'U')) && (order ==
           0)) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));
        if ((weight->name.size(1) == 0) || (weight->name[0] == 'U')) {
          //  unit weights
          wls_idx_0 = b_wls->rweights.size(0);
          b_wls->rweights.set_size(wls_idx_0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->rweights[i] = 1.0;
          }
        } else {
          char c;
          c = cv[static_cast<unsigned char>(weight->name[0]) & 127];
          if (c == 'I') {
            //  inverse distance
            wls_invdist_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else if (c == 'B') {
            //  Buhmann weights. All points share same parameters
            wls_buhmann_weights(b_wls->us, us.size(0), degree,
                                weight->params_shared, weight->params_pointwise,
                                b_wls->rweights);
          } else {
            boolean_T flag;
            flag = (c == 'E');

            //  Throw error if condition false
            //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

            mxAssert(flag,
                     "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.");

#else //MATLAB_MEX_FILE

            if (!flag) {
              fprintf(stderr,
                      "Weighting scheme must be unit, inverse, Buhmann, or WLS-ENO.\n");
              fflush(stderr);
            }

            assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

            //  WLS-ENO
            wls_eno_weights(b_wls->us, us.size(0), degree, us,
                            weight->params_shared, weight->params_pointwise,
                            b_wls->rweights);
          }
        }
      }

      //  Compute Vandermonde system and recompute DAG if needed
      gen_vander(b_wls->us, us.size(0), degree, order, b_wls->rweights, b_wls->V);
      nrblks = b_wls->V.size(1) / b_wls->stride;
      ncols = b_wls->V.size(0);

      //  Compact CVM if needed
      if ((order > 0) && (us.size(0) != b_wls->stride) && (us.size(0) !=
           b_wls->stride)) {
        //  Compact the storage of Vandermonde matrix
        //
        //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
        //
        //  Parameters
        //  ----------
        //     V:        confluent Vandermonde matrix
        //     npoints:  number of points
        //     nrblks:   number of row blocks
        //     stride:   length of each row block
        //     rowmajor: Array is in row-major
        //
        //  Returns
        //  -------
        //     V:       V with data in first npoints*nrblks rows
        trg = us.size(0);
        for (b_b = 2; b_b <= nrblks; b_b++) {
          src = (b_b - 1) * b_wls->stride;
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->V.size(0);
            for (j = 0; j < i; j++) {
              b_wls->V[trg + b_wls->V.size(1) * j] = b_wls->V[src +
                b_wls->V.size(1) * j];
            }

            src++;
            trg++;
          }
        }
      }

      b_wls->nrows = nrblks * us.size(0);
      b_wls->ncols = ncols;

      //  Omit rows in CVM if needed
      wls_idx_0 = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (wls_idx_0 <= u1) {
        u1 = wls_idx_0;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          wls_idx_0 = b_wls->V.size(0);
          for (i = 0; i < wls_idx_0; i++) {
            b_wls->V[b_i + b_wls->V.size(1) * i] = 0.0;
          }
        }
      }

      //  Perform QR with column pivoting
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      b_wls->rank = rrqr_factor(b_wls->V, thres, b_wls->nrows, ncols, b_wls->QR,
        b_wls->jpvt, b_wls->work, b_wls->dag);
      if ((b_wls->rweights.size(0) != 0) && (order > 0)) {
        //  Compute weights for derivatives
        if (order <= 2) {
          double s;
          int J;
          int blk;
          s = 1.0 / maxx_inv;
          for (blk = 0; blk <= dim; blk++) {
            J = (blk + 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] * s;
            }
          }

          if (order == 2) {
            s = 1.0 / (maxx_inv * maxx_inv);
            i = us.size(1) + 2;
            for (blk = i; blk <= nrblks; blk++) {
              J = (blk - 1) * b_wls->stride;
              for (j = 0; j <= npoints; j++) {
                b_wls->rweights[J + j] = b_wls->rweights[j] * s;
              }
            }
          }
        } else {
          //  Compute scaling factors for each block. Use wls.vdops as work space.
          //  It is important to `coder.ignoreConst` on fourth argument so that
          //  V is not optimized to be a 1-D array by MATLAB Coder.
          gen_vander(b_wls->hs_inv.data, b_wls->hs_inv.size, order, b_wls->V);
          for (int blk{2}; blk <= nrblks; blk++) {
            int J;
            J = (blk - 1) * b_wls->stride;
            for (j = 0; j <= npoints; j++) {
              b_wls->rweights[J + j] = b_wls->rweights[j] / b_wls->V
                [b_wls->V.size(1) * (blk - 1)];
            }
          }
        }

        if (us.size(0) != b_wls->stride) {
          //  Compact the storage of Vandermonde matrix
          //
          //   V = compact_vander(V, npoints, nrblks, stride, rowmajor)
          //
          //  Parameters
          //  ----------
          //     V:        confluent Vandermonde matrix
          //     npoints:  number of points
          //     nrblks:   number of row blocks
          //     stride:   length of each row block
          //     rowmajor: Array is in row-major
          //
          //  Returns
          //  -------
          //     V:       V with data in first npoints*nrblks rows
          trg = us.size(0);
          for (b_b = 2; b_b <= nrblks; b_b++) {
            src = (b_b - 1) * b_wls->stride;
            for (b_i = 0; b_i <= npoints; b_i++) {
              b_wls->rweights[trg + b_i] = b_wls->rweights[src + b_i];
            }

            trg = (trg + npoints) + 1;
          }
        }
      }
    }
  }

  static inline
  void wls_var_bilap(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                      coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                      &vdops, ::coder::array<double, 2U> &result)
  {
    int bilap_size_idx_1;
    int iPoint;
    int iRow;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];
    boolean_T flag;

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_bilap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      bilap_size_idx_1 = 1;
      bilap_data[0] = 4;
      break;

     case 2:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 6;
        bilap_data[1] = 7;
        bilap_data[2] = 7;
        bilap_data[3] = 8;
      } else {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 9;
        bilap_data[1] = 10;
        bilap_data[2] = 10;
        bilap_data[3] = 11;
      }
      break;

     default:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 9;
        for (u1 = 0; u1 < 9; u1++) {
          bilap_data[u1] = iv[u1];
        }
      } else {
        bilap_size_idx_1 = 9;
        for (u1 = 0; u1 < 9; u1++) {
          bilap_data[u1] = iv1[u1];
        }
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u1 = 0; u1 < nrows_vdops; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset
              + iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
              iPoint] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow <= nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(1, fs.size(1));
      u0 = fs.size(1);
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_bilap(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, ::coder::array<double, 1U> &vdops)
  {
    int bilap_size_idx_1;
    int iPoint;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_bilap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      bilap_size_idx_1 = 1;
      bilap_data[0] = 4;
      break;

     case 2:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 6;
        bilap_data[1] = 7;
        bilap_data[2] = 7;
        bilap_data[3] = 8;
      } else {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 9;
        bilap_data[1] = 10;
        bilap_data[2] = 10;
        bilap_data[3] = 11;
      }
      break;

     default:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 9;
        for (u0 = 0; u0 < 9; u0++) {
          bilap_data[u0] = iv[u0];
        }
      } else {
        bilap_size_idx_1 = 9;
        for (u0 = 0; u0 < 9; u0++) {
          bilap_data[u0] = iv1[u0];
        }
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset +
            iPoint) + b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }
  }

  static inline
  void wls_var_bilap(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                      array<double, 1U> &vdops)
  {
    int bilap_size_idx_1;
    int iPoint;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];
    boolean_T flag;

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_bilap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_bilap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      bilap_size_idx_1 = 1;
      bilap_data[0] = 4;
      break;

     case 2:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 6;
        bilap_data[1] = 7;
        bilap_data[2] = 7;
        bilap_data[3] = 8;
      } else {
        bilap_size_idx_1 = 4;
        bilap_data[0] = 9;
        bilap_data[1] = 10;
        bilap_data[2] = 10;
        bilap_data[3] = 11;
      }
      break;

     default:
      if (b_wls->degree > 0) {
        bilap_size_idx_1 = 9;
        for (u0 = 0; u0 < 9; u0++) {
          bilap_data[u0] = iv[u0];
        }
      } else {
        bilap_size_idx_1 = 9;
        for (u0 = 0; u0 < 9; u0++) {
          bilap_data[u0] = iv1[u0];
        }
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset
              + iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
              iPoint] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }
  }

  static inline
  void wls_var_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                     &vdops, double result_data[], int result_size[2])
  {
    double c_vdops;
    double e_vdops;
    double f_vdops;
    int i;
    int u1;

    //  Variational curl operators as weighted sum at quadrature points in 3D
    //
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-3)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or 3. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-3 will be used
    //                to scale the three components of the curl, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-3).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-9.
    //                To apply the operator, use
    //                   [sum(fs(:,2:3) .* vdops(:,2:3), 'all');
    //                    sum(fs(:,[1,3]) .* vdops(:,[4,6]), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,7:8), 'all')];
    //     result:    Computed solution of size 3-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  compute vdops and reorganize
    if (ws.size(1) <= 1) {
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int dim{0}; dim < nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(9, nrows_vdops);
      u0 = nrows_vdops * 9;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = 9;
      i = 9;
      while ((i > 0) && (iv4[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv4[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, 9);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + 9 * iRow] = vdops[k + 9 * iRow] * b_wls->rweights[iRow];
          }
        }
      }

      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double d_vdops;
        b_vdops = vdops[9 * i + 2];
        c_vdops = vdops[9 * i + 1];
        d_vdops = vdops[9 * i];
        vdops[9 * i] = 0.0;
        vdops[9 * i + 1] = -b_vdops;
        vdops[9 * i + 2] = c_vdops;
        vdops[9 * i + 3] = b_vdops;
        vdops[9 * i + 4] = 0.0;
        vdops[9 * i + 5] = -d_vdops;
        vdops[9 * i + 6] = -c_vdops;
        vdops[9 * i + 7] = d_vdops;
        vdops[9 * i + 8] = 0.0;
      }
    } else {
      int iPoint;
      int iWeight;
      int j;
      int lenWs;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      boolean_T flag;

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (9 / lenWs * lenWs == 9);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int dim{0}; dim < nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(9, nrows_vdops);
      u0 = nrows_vdops * 9;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = 9;
      i = 9;
      while ((i > 0) && (iv5[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv5[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if (ws.size(0) == 0) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                  [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                  iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }

        if (iWeight == lenWs) {
          iWeight = 1;
        } else {
          iWeight++;
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, 9);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + 9 * iRow] = vdops[k + 9 * iRow] * b_wls->rweights[iRow];
          }
        }
      }

      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double d_vdops;
        double g_vdops;
        b_vdops = vdops[9 * i + 3];
        d_vdops = vdops[9 * i];
        e_vdops = vdops[9 * i + 1];
        f_vdops = vdops[9 * i + 4];
        c_vdops = vdops[9 * i + 5];
        g_vdops = vdops[9 * i + 2];
        vdops[9 * i] = 0.0;
        vdops[9 * i + 1] = -b_vdops;
        vdops[9 * i + 2] = d_vdops;
        vdops[9 * i + 3] = e_vdops;
        vdops[9 * i + 4] = 0.0;
        vdops[9 * i + 5] = -f_vdops;
        vdops[9 * i + 6] = -c_vdops;
        vdops[9 * i + 7] = g_vdops;
        vdops[9 * i + 8] = 0.0;
      }
    }

    //  compute output value
    if ((fs.size(0) != 0) && (fs.size(1) != 0)) {
      result_size[1] = 1;
      result_size[0] = 3;
      result_data[0] = 0.0;
      result_data[1] = 0.0;
      result_data[2] = 0.0;
      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        c_vdops = fs[fs.size(1) * i + 1];
        e_vdops = fs[fs.size(1) * i + 2];
        result_data[0] = (result_data[0] + vdops[9 * i + 1] * c_vdops) + vdops[9
          * i + 2] * e_vdops;
        f_vdops = fs[fs.size(1) * i];
        result_data[1] = (result_data[1] + f_vdops * vdops[9 * i + 3]) + e_vdops
          * vdops[9 * i + 5];
        result_data[2] = (result_data[2] + f_vdops * vdops[9 * i + 6]) + c_vdops
          * vdops[9 * i + 7];
      }
    } else {
      result_size[1] = 0;
      result_size[0] = 3;
    }
  }

  static inline
  void wls_var_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int i;
    int iPoint;
    int j;
    int nDims;
    int nOps;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;

    //  Variational curl operators as weighted sum at quadrature points in 3D
    //
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-3)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or 3. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-3 will be used
    //                to scale the three components of the curl, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-3).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-9.
    //                To apply the operator, use
    //                   [sum(fs(:,2:3) .* vdops(:,2:3), 'all');
    //                    sum(fs(:,[1,3]) .* vdops(:,[4,6]), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,7:8), 'all')];
    //     result:    Computed solution of size 3-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  compute vdops and reorganize
    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(9, nrows_vdops);
    u0 = nrows_vdops * 9;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    nOps = 9;
    i = 9;
    while ((i > 0) && (iv4[i - 1] == 0)) {
      nOps--;
      i--;
    }

    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = iv4[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, 9);
    for (i = 0; i < nrows_vdops; i++) {
      for (j = 0; j < 9; j++) {
        vdops[j + 9 * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k < nOps; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + 9 * iRow] = vdops[k + 9 * iRow] * b_wls->rweights[iRow];
        }
      }
    }

    u1 = b_wls->nrows;
    for (i = 0; i < u1; i++) {
      double b_vdops;
      double c_vdops;
      double d;
      b_vdops = vdops[9 * i + 2];
      d = vdops[9 * i + 1];
      c_vdops = vdops[9 * i];
      vdops[9 * i] = 0.0;
      vdops[9 * i + 1] = -b_vdops;
      vdops[9 * i + 2] = d;
      vdops[9 * i + 3] = b_vdops;
      vdops[9 * i + 4] = 0.0;
      vdops[9 * i + 5] = -c_vdops;
      vdops[9 * i + 6] = -d;
      vdops[9 * i + 7] = c_vdops;
      vdops[9 * i + 8] = 0.0;
    }

    //  compute output value
  }

  static inline
  void wls_var_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    //  Variational curl operators as weighted sum at quadrature points in 3D
    //
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-3)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or 3. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-3 will be used
    //                to scale the three components of the curl, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-3).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-9.
    //                To apply the operator, use
    //                   [sum(fs(:,2:3) .* vdops(:,2:3), 'all');
    //                    sum(fs(:,[1,3]) .* vdops(:,[4,6]), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,7:8), 'all')];
    //     result:    Computed solution of size 3-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  compute vdops and reorganize
    if (ws.size(1) <= 1) {
      int i;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int dim{0}; dim < nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(9, nrows_vdops);
      u0 = nrows_vdops * 9;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = 9;
      i = 9;
      while ((i > 0) && (iv4[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv4[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, 9);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + 9 * iRow] = vdops[k + 9 * iRow] * b_wls->rweights[iRow];
          }
        }
      }

      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        double d_vdops;
        b_vdops = vdops[9 * i + 2];
        c_vdops = vdops[9 * i + 1];
        d_vdops = vdops[9 * i];
        vdops[9 * i] = 0.0;
        vdops[9 * i + 1] = -b_vdops;
        vdops[9 * i + 2] = c_vdops;
        vdops[9 * i + 3] = b_vdops;
        vdops[9 * i + 4] = 0.0;
        vdops[9 * i + 5] = -d_vdops;
        vdops[9 * i + 6] = -c_vdops;
        vdops[9 * i + 7] = d_vdops;
        vdops[9 * i + 8] = 0.0;
      }
    } else {
      int i;
      int iPoint;
      int iWeight;
      int j;
      int lenWs;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;
      boolean_T flag;

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (9 / lenWs * lenWs == 9);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int dim{0}; dim < nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(9, nrows_vdops);
      u0 = nrows_vdops * 9;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = 9;
      i = 9;
      while ((i > 0) && (iv5[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv5[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if (ws.size(0) == 0) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                  [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                  iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }

        if (iWeight == lenWs) {
          iWeight = 1;
        } else {
          iWeight++;
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, 9);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + 9 * iRow] = vdops[k + 9 * iRow] * b_wls->rweights[iRow];
          }
        }
      }

      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        double d_vdops;
        double e_vdops;
        double f_vdops;
        double g_vdops;
        b_vdops = vdops[9 * i + 3];
        d_vdops = vdops[9 * i];
        c_vdops = vdops[9 * i + 1];
        e_vdops = vdops[9 * i + 4];
        f_vdops = vdops[9 * i + 5];
        g_vdops = vdops[9 * i + 2];
        vdops[9 * i] = 0.0;
        vdops[9 * i + 1] = -b_vdops;
        vdops[9 * i + 2] = d_vdops;
        vdops[9 * i + 3] = c_vdops;
        vdops[9 * i + 4] = 0.0;
        vdops[9 * i + 5] = -e_vdops;
        vdops[9 * i + 6] = -f_vdops;
        vdops[9 * i + 7] = g_vdops;
        vdops[9 * i + 8] = 0.0;
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::coder::array<
    double, 2U> &fs, ::coder::array<double, 2U> &vdops, double result_data[],
    int result_size[2])
  {
    double b_vdops;
    double d;
    double d1;
    int dim;
    int i;
    int u0;
    int u1;
    signed char grad_div_data[18];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_curl_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size_idx_1;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;

      //  All components share the same weight, we need to compute Hessian
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size_idx_1 = 1;
        hess_data[0] = 0;
        break;

       case 2:
        grad_div_size_idx_1 = 4;
        hess_data[0] = 4;
        hess_data[1] = 5;
        hess_data[2] = 6;
        hess_data[3] = 0;
        break;

       default:
        grad_div_size_idx_1 = 9;
        for (u1 = 0; u1 < 9; u1++) {
          hess_data[u1] = iv2[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size_idx_1, nrows_vdops);
      u0 = nrows_vdops * grad_div_size_idx_1;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size_idx_1;
      i = grad_div_size_idx_1;
      while ((i > 0) && (hess_data[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size_idx_1);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size_idx_1; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double c_vdops;
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i + 2];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = b_vdops;
          vdops[vdops.size(1) * i + 3] = -c_vdops;
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double c_vdops;
          double d_vdops;
          double e_vdops;
          b_vdops = vdops[vdops.size(1) * i + 2];
          c_vdops = vdops[vdops.size(1) * i + 5];
          d_vdops = vdops[vdops.size(1) * i + 1];
          e_vdops = vdops[vdops.size(1) * i + 3];
          d = vdops[vdops.size(1) * i];
          d1 = vdops[vdops.size(1) * i + 4];
          vdops[vdops.size(1) * i] = -b_vdops - c_vdops;
          vdops[vdops.size(1) * i + 1] = d_vdops;
          vdops[vdops.size(1) * i + 2] = e_vdops;
          vdops[vdops.size(1) * i + 3] = d_vdops;
          vdops[vdops.size(1) * i + 4] = -d - c_vdops;
          vdops[vdops.size(1) * i + 5] = d1;
          vdops[vdops.size(1) * i + 6] = e_vdops;
          vdops[vdops.size(1) * i + 7] = d1;
          vdops[vdops.size(1) * i + 8] = -d - b_vdops;
        }
      }
    } else {
      int grad_div_size_idx_0;
      int grad_div_size_idx_1;
      int iPoint;
      int j;
      int lenWs;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      boolean_T flag;

      //  Each component has its own weight, so we need to compute all components
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size_idx_1 = 1;
        grad_div_size_idx_0 = 1;
        grad_div_data[0] = 0;
        break;

       case 2:
        grad_div_size_idx_1 = 1;
        grad_div_size_idx_0 = 4;
        grad_div_data[0] = 6;
        grad_div_data[1] = 5;
        grad_div_data[2] = 5;
        grad_div_data[3] = 4;
        break;

       default:
        grad_div_size_idx_1 = 2;
        grad_div_size_idx_0 = 9;
        for (u1 = 0; u1 < 18; u1++) {
          grad_div_data[u1] = iv3[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (grad_div_size_idx_0 / lenWs * lenWs == grad_div_size_idx_0);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size_idx_0, nrows_vdops);
      u0 = nrows_vdops * grad_div_size_idx_0;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size_idx_0;
      i = grad_div_size_idx_0;
      while ((i > 0) && (grad_div_data[grad_div_size_idx_1 * (i - 1)] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      for (int jDiff{0}; jDiff < grad_div_size_idx_1; jDiff++) {
        int iWeight;
        iWeight = 1;

        //  Loop through the operators
        for (int iOp{0}; iOp < nOps; iOp++) {
          signed char b_i;

          //  Skip padded zeros in the differential operator
          b_i = grad_div_data[jDiff + grad_div_size_idx_1 * iOp];
          if (b_i > 0) {
            int offset;
            offset = (b_i - 1) * stride;

            //  Sum up monomials weighted by weights for each component
            for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
              j = b_wls->jpvt[iMonomial] - 1;
              if (ws.size(0) == 0) {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                    b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                    b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
                }
              } else {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                    b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                    [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                    iPoint) + b_wls->V.size(1) * j];
                }
              }
            }
          }

          if (iWeight == lenWs) {
            iWeight = 1;
          } else {
            iWeight++;
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size_idx_0);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size_idx_0; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      //  Flip the signs of diagonal entries
      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i + 3] = -vdops[vdops.size(1) * i + 3];
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i + 4] = -vdops[vdops.size(1) * i + 4];
          vdops[vdops.size(1) * i + 8] = -vdops[vdops.size(1) * i + 8];
        }
      }
    }

    //  compute output value
    if ((fs.size(0) != 0) && (fs.size(1) != 0)) {
      result_size[1] = 1;
      result_size[0] = 3;
      u0 = static_cast<signed char>(dim);
      if (0 <= u0 - 1) {
        std::memset(&result_data[0], 0, u0 * sizeof(double));
      }

      switch (dim) {
       case 1:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          result_data[0] += vdops[vdops.size(1) * i] * fs[fs.size(1) * i];
        }
        break;

       case 2:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          d = fs[fs.size(1) * i];
          d1 = fs[fs.size(1) * i + 1];
          result_data[0] = (result_data[0] + vdops[vdops.size(1) * i] * d) +
            vdops[vdops.size(1) * i + 1] * d1;
          result_data[1] = (result_data[1] + d * vdops[vdops.size(1) * i + 2]) +
            d1 * vdops[vdops.size(1) * i + 3];
        }
        break;

       case 3:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          d = fs[fs.size(1) * i];
          d1 = fs[fs.size(1) * i + 1];
          b_vdops = fs[fs.size(1) * i + 2];
          result_data[0] = ((result_data[0] + vdops[vdops.size(1) * i] * d) +
                            vdops[vdops.size(1) * i + 1] * d1) +
            vdops[vdops.size(1) * i + 2] * b_vdops;
          result_data[1] = ((result_data[1] + d * vdops[vdops.size(1) * i + 3])
                            + d1 * vdops[vdops.size(1) * i + 4]) + b_vdops *
            vdops[vdops.size(1) * i + 5];
          result_data[2] = ((result_data[2] + d * vdops[vdops.size(1) * i + 6])
                            + d1 * vdops[vdops.size(1) * i + 7]) + b_vdops *
            vdops[vdops.size(1) * i + 8];
        }
        break;
      }
    } else {
      result_size[1] = 0;
      result_size[0] = 3;
    }
  }

  static inline
  void wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int dim;
    int hess_size;
    int i;
    int iPoint;
    int j;
    int nDims;
    int nOps;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_curl_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    //  All components share the same weight, we need to compute Hessian
    switch (b_wls->us.size(1)) {
     case 1:
      hess_size = 1;
      hess_data[0] = 0;
      break;

     case 2:
      hess_size = 4;
      hess_data[0] = 4;
      hess_data[1] = 5;
      hess_data[2] = 6;
      hess_data[3] = 0;
      break;

     default:
      hess_size = 9;
      for (u1 = 0; u1 < 9; u1++) {
        hess_data[u1] = iv2[u1];
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int b_dim{0}; b_dim < nDims; b_dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(hess_size, nrows_vdops);
    u0 = nrows_vdops * hess_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    nOps = hess_size;
    i = hess_size;
    while ((i > 0) && (hess_data[i - 1] == 0)) {
      nOps--;
      i--;
    }

    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = hess_data[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, hess_size);
    for (i = 0; i < nrows_vdops; i++) {
      for (j = 0; j < hess_size; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k < nOps; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (dim == 2) {
      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        b_vdops = vdops[vdops.size(1) * i + 1];
        c_vdops = vdops[vdops.size(1) * i];
        vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i + 2];
        vdops[vdops.size(1) * i + 1] = b_vdops;
        vdops[vdops.size(1) * i + 2] = b_vdops;
        vdops[vdops.size(1) * i + 3] = -c_vdops;
      }
    } else if (dim == 3) {
      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        double d;
        double d1;
        double d_vdops;
        double e_vdops;
        b_vdops = vdops[vdops.size(1) * i + 2];
        c_vdops = vdops[vdops.size(1) * i + 5];
        d_vdops = vdops[vdops.size(1) * i + 1];
        e_vdops = vdops[vdops.size(1) * i + 3];
        d = vdops[vdops.size(1) * i];
        d1 = vdops[vdops.size(1) * i + 4];
        vdops[vdops.size(1) * i] = -b_vdops - c_vdops;
        vdops[vdops.size(1) * i + 1] = d_vdops;
        vdops[vdops.size(1) * i + 2] = e_vdops;
        vdops[vdops.size(1) * i + 3] = d_vdops;
        vdops[vdops.size(1) * i + 4] = -d - c_vdops;
        vdops[vdops.size(1) * i + 5] = d1;
        vdops[vdops.size(1) * i + 6] = e_vdops;
        vdops[vdops.size(1) * i + 7] = d1;
        vdops[vdops.size(1) * i + 8] = -d - b_vdops;
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::array<double, 2U>
    &vdops)
  {
    int dim;
    signed char grad_div_data[18];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts)
    //  [wls, vdops] = wls_var_curl_curl(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_curl_curl(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_grad_div, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size_idx_1;
      int i;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

      //  All components share the same weight, we need to compute Hessian
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size_idx_1 = 1;
        hess_data[0] = 0;
        break;

       case 2:
        grad_div_size_idx_1 = 4;
        hess_data[0] = 4;
        hess_data[1] = 5;
        hess_data[2] = 6;
        hess_data[3] = 0;
        break;

       default:
        grad_div_size_idx_1 = 9;
        for (u1 = 0; u1 < 9; u1++) {
          hess_data[u1] = iv2[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size_idx_1, nrows_vdops);
      u0 = nrows_vdops * grad_div_size_idx_1;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size_idx_1;
      i = grad_div_size_idx_1;
      while ((i > 0) && (hess_data[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size_idx_1);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size_idx_1; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double b_vdops;
          double c_vdops;
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i + 2];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = b_vdops;
          vdops[vdops.size(1) * i + 3] = -c_vdops;
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double b_vdops;
          double c_vdops;
          double d;
          double d1;
          double d_vdops;
          double e_vdops;
          b_vdops = vdops[vdops.size(1) * i + 2];
          c_vdops = vdops[vdops.size(1) * i + 5];
          d_vdops = vdops[vdops.size(1) * i + 1];
          e_vdops = vdops[vdops.size(1) * i + 3];
          d = vdops[vdops.size(1) * i];
          d1 = vdops[vdops.size(1) * i + 4];
          vdops[vdops.size(1) * i] = -b_vdops - c_vdops;
          vdops[vdops.size(1) * i + 1] = d_vdops;
          vdops[vdops.size(1) * i + 2] = e_vdops;
          vdops[vdops.size(1) * i + 3] = d_vdops;
          vdops[vdops.size(1) * i + 4] = -d - c_vdops;
          vdops[vdops.size(1) * i + 5] = d1;
          vdops[vdops.size(1) * i + 6] = e_vdops;
          vdops[vdops.size(1) * i + 7] = d1;
          vdops[vdops.size(1) * i + 8] = -d - b_vdops;
        }
      }
    } else {
      int grad_div_size_idx_0;
      int grad_div_size_idx_1;
      int i;
      int iPoint;
      int j;
      int lenWs;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;
      boolean_T flag;

      //  Each component has its own weight, so we need to compute all components
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size_idx_1 = 1;
        grad_div_size_idx_0 = 1;
        grad_div_data[0] = 0;
        break;

       case 2:
        grad_div_size_idx_1 = 1;
        grad_div_size_idx_0 = 4;
        grad_div_data[0] = 6;
        grad_div_data[1] = 5;
        grad_div_data[2] = 5;
        grad_div_data[3] = 4;
        break;

       default:
        grad_div_size_idx_1 = 2;
        grad_div_size_idx_0 = 9;
        for (u1 = 0; u1 < 18; u1++) {
          grad_div_data[u1] = iv3[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (grad_div_size_idx_0 / lenWs * lenWs == grad_div_size_idx_0);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size_idx_0, nrows_vdops);
      u0 = nrows_vdops * grad_div_size_idx_0;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size_idx_0;
      i = grad_div_size_idx_0;
      while ((i > 0) && (grad_div_data[grad_div_size_idx_1 * (i - 1)] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      for (int jDiff{0}; jDiff < grad_div_size_idx_1; jDiff++) {
        int iWeight;
        iWeight = 1;

        //  Loop through the operators
        for (int iOp{0}; iOp < nOps; iOp++) {
          signed char b_i;

          //  Skip padded zeros in the differential operator
          b_i = grad_div_data[jDiff + grad_div_size_idx_1 * iOp];
          if (b_i > 0) {
            int offset;
            offset = (b_i - 1) * stride;

            //  Sum up monomials weighted by weights for each component
            for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
              j = b_wls->jpvt[iMonomial] - 1;
              if (ws.size(0) == 0) {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                    b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                    b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
                }
              } else {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                    b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                    [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                    iPoint) + b_wls->V.size(1) * j];
                }
              }
            }
          }

          if (iWeight == lenWs) {
            iWeight = 1;
          } else {
            iWeight++;
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size_idx_0);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size_idx_0; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      //  Flip the signs of diagonal entries
      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i + 3] = -vdops[vdops.size(1) * i + 3];
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          vdops[vdops.size(1) * i] = -vdops[vdops.size(1) * i];
          vdops[vdops.size(1) * i + 4] = -vdops[vdops.size(1) * i + 4];
          vdops[vdops.size(1) * i + 8] = -vdops[vdops.size(1) * i + 8];
        }
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_div(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                    coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                    &vdops, double result_data[], int result_size[2])
  {
    int grad_size;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];
    boolean_T flag;

    //  Compute variational divergence operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_div(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use sum(fs .* vdops, 'all').
    //     result:    Computed solution of size 1-by-1 (scalar).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (grad_size / lenWs * lenWs == grad_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result_size[1] = 0;
      result_size[0] = 0;
    } else {
      int iFunc;
      u0 = grad_size / quad_pnts.size(1);
      result_size[1] = 1;
      result_size[0] = static_cast<signed char>(u0);
      u0 = static_cast<signed char>(u0);
      if (0 <= u0 - 1) {
        std::memset(&result_data[0], 0, u0 * sizeof(double));
      }

      //  Compute solution
      iFunc = 1;
      iOp = 0;
      for (int iDiff{0}; iDiff <= nDiff; iDiff++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          result_data[iOp] += fs[(iFunc + fs.size(1) * iRow) - 1] * vdops[iDiff
            + vdops.size(1) * iRow];
        }

        if (iOp + 1 == result_size[0]) {
          iOp = 0;
        } else {
          iOp++;
        }

        if (iFunc == fs.size(1)) {
          iFunc = 1;
        } else {
          iFunc++;
        }
      }
    }
  }

  static inline
  void wls_var_div(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int grad_size;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational divergence operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_div(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use sum(fs .* vdops, 'all').
    //     result:    Computed solution of size 1-by-1 (scalar).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
            vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
            iPoint) + b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_div(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                    array<double, 2U> &vdops)
  {
    int grad_size;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];
    boolean_T flag;

    //  Compute variational divergence operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_div(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use sum(fs .* vdops, 'all').
    //     result:    Computed solution of size 1-by-1 (scalar).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (grad_size / lenWs * lenWs == grad_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_func(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int iPoint;
    int iRow;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;
    boolean_T flag;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    //
    //  [wls, vdops] = wls_var_func(wls, pnts) computes sum of fittings at quadrature
    //               points (i.e., with unit weights). Can also pass a single point.
    //  [wls, vdops] = wls_var_func(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_func(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should be the product
    //                of the quadrature weight, Jacobian determinant, and value of a
    //                weighting function (e.g., a test function or the derivative
    //                of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_grad, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u1 = 0; u1 < nrows_vdops; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    //  Skip padded zeros in the differential operator
    //  Sum up monomials weighted by weights for each component
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial] - 1;
      if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[iPoint
            + b_wls->V.size(1) * j];
        }
      } else {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
            iPoint] * b_wls->V[iPoint + b_wls->V.size(1) * j];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow <= nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(1, fs.size(1));
      u0 = fs.size(1);
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_func(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 1U> &vdops)
  {
    int iPoint;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    //
    //  [wls, vdops] = wls_var_func(wls, pnts) computes sum of fittings at quadrature
    //               points (i.e., with unit weights). Can also pass a single point.
    //  [wls, vdops] = wls_var_func(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_func(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should be the product
    //                of the quadrature weight, Jacobian determinant, and value of a
    //                weighting function (e.g., a test function or the derivative
    //                of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_grad, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    //  Skip padded zeros in the differential operator
    //  Sum up monomials weighted by weights for each component
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[iPoint
          + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }
  }

  static inline
  void wls_var_func(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 1U> &vdops)
  {
    int iPoint;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;
    boolean_T flag;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    //
    //  [wls, vdops] = wls_var_func(wls, pnts) computes sum of fittings at quadrature
    //               points (i.e., with unit weights). Can also pass a single point.
    //  [wls, vdops] = wls_var_func(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_func(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should be the product
    //                of the quadrature weight, Jacobian determinant, and value of a
    //                weighting function (e.g., a test function or the derivative
    //                of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_grad, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    //  Skip padded zeros in the differential operator
    //  Sum up monomials weighted by weights for each component
    for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial] - 1;
      if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[iPoint
            + b_wls->V.size(1) * j];
        }
      } else {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
            iPoint] * b_wls->V[iPoint + b_wls->V.size(1) * j];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }
  }

  static inline
  void wls_var_grad(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int grad_size;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];
    boolean_T flag;

    //  Compute variational gradient operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_grad(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (eg., a test function).
    //                Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (grad_size / lenWs * lenWs == grad_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(grad_size, fs.size(1));
      u0 = fs.size(1) * grad_size;
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        iOp = 0;
        for (int iDiff{0}; iDiff <= nDiff; iDiff++) {
          for (iRow = 0; iRow <= nrows; iRow++) {
            result[iFunc + result.size(1) * iOp] = result[iFunc + result.size(1)
              * iOp] + fs[iFunc + fs.size(1) * iRow] * vdops[iDiff + vdops.size
              (1) * iRow];
          }

          if (iOp + 1 == result.size(0)) {
            iOp = 0;
          } else {
            iOp++;
          }
        }
      }
    }
  }

  static inline
  void wls_var_grad(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int grad_size;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational gradient operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_grad(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (eg., a test function).
    //                Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
            vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
            iPoint) + b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_grad(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    int grad_size;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];
    boolean_T flag;

    //  Compute variational gradient operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_grad(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (eg., a test function).
    //                Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-1. To apply the
    //                operator, use vdops' * fs
    //     result:    Computed solution of size d-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_div, wls_var_curl, wls_var_hess,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      grad_size = 1;
      grad_data[0] = 2;
      break;

     case 2:
      grad_size = 2;
      grad_data[0] = 2;
      grad_data[1] = 3;
      break;

     default:
      grad_size = 3;
      grad_data[0] = 2;
      grad_data[1] = 3;
      grad_data[2] = 4;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (grad_size / lenWs * lenWs == grad_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(grad_size, nrows_vdops);
    u0 = nrows_vdops * grad_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, grad_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_grad_div(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::coder::array<
    double, 2U> &fs, ::coder::array<double, 2U> &vdops, double result_data[],
    int result_size[2])
  {
    double b_vdops;
    double c_vdops;
    double d;
    int dim;
    int i;
    int u0;
    int u1;
    signed char grad_div_data[9];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_grad_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;

      //  All components share the same weight, we need to compute Hessian
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size = 1;
        hess_data[0] = 3;
        break;

       case 2:
        grad_div_size = 4;
        hess_data[0] = 4;
        hess_data[1] = 5;
        hess_data[2] = 6;
        hess_data[3] = 0;
        break;

       default:
        grad_div_size = 9;
        for (u1 = 0; u1 < 9; u1++) {
          hess_data[u1] = iv2[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size, nrows_vdops);
      u0 = nrows_vdops * grad_div_size;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size;
      i = grad_div_size;
      while ((i > 0) && (hess_data[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i + 2];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = b_vdops;
          vdops[vdops.size(1) * i + 3] = c_vdops;
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double d_vdops;
          double e_vdops;
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i + 3];
          d_vdops = vdops[vdops.size(1) * i + 2];
          d = vdops[vdops.size(1) * i + 4];
          e_vdops = vdops[vdops.size(1) * i + 5];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = c_vdops;
          vdops[vdops.size(1) * i + 3] = b_vdops;
          vdops[vdops.size(1) * i + 4] = d_vdops;
          vdops[vdops.size(1) * i + 5] = d;
          vdops[vdops.size(1) * i + 6] = c_vdops;
          vdops[vdops.size(1) * i + 7] = d;
          vdops[vdops.size(1) * i + 8] = e_vdops;
        }
      }
    } else {
      int grad_div_size;
      int iPoint;
      int iWeight;
      int j;
      int lenWs;
      int nDiff;
      int nDims;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      boolean_T flag;

      //  Each component has its own weight, we need to compute all components
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size = 1;
        grad_div_data[0] = 3;
        break;

       case 2:
        grad_div_size = 4;
        grad_div_data[0] = 4;
        grad_div_data[1] = 5;
        grad_div_data[2] = 5;
        grad_div_data[3] = 6;
        break;

       default:
        grad_div_size = 9;
        for (u1 = 0; u1 < 9; u1++) {
          grad_div_data[u1] = iv6[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      nDiff = grad_div_size - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (grad_div_size / lenWs * lenWs == grad_div_size);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size, nrows_vdops);
      u0 = nrows_vdops * grad_div_size;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      //  Summing up rows in the differential operator
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp <= nDiff; iOp++) {
        int offset;

        //  Skip padded zeros in the differential operator
        offset = (grad_div_data[iOp] - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
          j = b_wls->jpvt[iMonomial] - 1;
          if (ws.size(0) == 0) {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V
                [(offset + iPoint) + b_wls->V.size(1) * j];
            }
          } else {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                iPoint) + b_wls->V.size(1) * j];
            }
          }
        }

        if (iWeight == lenWs) {
          iWeight = 1;
        } else {
          iWeight++;
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops,
                   grad_div_size);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, grad_div_size, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k <= nDiff; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }
    }

    //  compute output value
    if ((fs.size(0) != 0) && (fs.size(1) != 0)) {
      result_size[1] = 1;
      result_size[0] = 3;
      u0 = static_cast<signed char>(dim);
      if (0 <= u0 - 1) {
        std::memset(&result_data[0], 0, u0 * sizeof(double));
      }

      switch (dim) {
       case 1:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          result_data[0] += vdops[vdops.size(1) * i] * fs[fs.size(1) * i];
        }
        break;

       case 2:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          d = fs[fs.size(1) * i];
          b_vdops = fs[fs.size(1) * i + 1];
          result_data[0] = (result_data[0] + vdops[vdops.size(1) * i] * d) +
            vdops[vdops.size(1) * i + 1] * b_vdops;
          result_data[1] = (result_data[1] + d * vdops[vdops.size(1) * i + 2]) +
            b_vdops * vdops[vdops.size(1) * i + 3];
        }
        break;

       case 3:
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          d = fs[fs.size(1) * i];
          b_vdops = fs[fs.size(1) * i + 1];
          c_vdops = fs[fs.size(1) * i + 2];
          result_data[0] = ((result_data[0] + vdops[vdops.size(1) * i] * d) +
                            vdops[vdops.size(1) * i + 1] * b_vdops) +
            vdops[vdops.size(1) * i + 2] * c_vdops;
          result_data[1] = ((result_data[1] + d * vdops[vdops.size(1) * i + 3])
                            + b_vdops * vdops[vdops.size(1) * i + 4]) + c_vdops *
            vdops[vdops.size(1) * i + 5];
          result_data[2] = ((result_data[2] + d * vdops[vdops.size(1) * i + 6])
                            + b_vdops * vdops[vdops.size(1) * i + 7]) + c_vdops *
            vdops[vdops.size(1) * i + 8];
        }
        break;
      }
    } else {
      result_size[1] = 0;
      result_size[0] = 3;
    }
  }

  static inline
  void wls_var_grad_div(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int dim;
    int hess_size;
    int i;
    int iPoint;
    int j;
    int nDims;
    int nOps;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_grad_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    //  All components share the same weight, we need to compute Hessian
    switch (b_wls->us.size(1)) {
     case 1:
      hess_size = 1;
      hess_data[0] = 3;
      break;

     case 2:
      hess_size = 4;
      hess_data[0] = 4;
      hess_data[1] = 5;
      hess_data[2] = 6;
      hess_data[3] = 0;
      break;

     default:
      hess_size = 9;
      for (u1 = 0; u1 < 9; u1++) {
        hess_data[u1] = iv2[u1];
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int b_dim{0}; b_dim < nDims; b_dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(hess_size, nrows_vdops);
    u0 = nrows_vdops * hess_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    nOps = hess_size;
    i = hess_size;
    while ((i > 0) && (hess_data[i - 1] == 0)) {
      nOps--;
      i--;
    }

    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = hess_data[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, hess_size);
    for (i = 0; i < nrows_vdops; i++) {
      for (j = 0; j < hess_size; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k < nOps; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (dim == 2) {
      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        b_vdops = vdops[vdops.size(1) * i + 1];
        c_vdops = vdops[vdops.size(1) * i + 2];
        vdops[vdops.size(1) * i + 1] = b_vdops;
        vdops[vdops.size(1) * i + 2] = b_vdops;
        vdops[vdops.size(1) * i + 3] = c_vdops;
      }
    } else if (dim == 3) {
      u1 = b_wls->nrows;
      for (i = 0; i < u1; i++) {
        double b_vdops;
        double c_vdops;
        double d;
        double d_vdops;
        double e_vdops;
        b_vdops = vdops[vdops.size(1) * i + 1];
        c_vdops = vdops[vdops.size(1) * i + 3];
        d_vdops = vdops[vdops.size(1) * i + 2];
        d = vdops[vdops.size(1) * i + 4];
        e_vdops = vdops[vdops.size(1) * i + 5];
        vdops[vdops.size(1) * i + 1] = b_vdops;
        vdops[vdops.size(1) * i + 2] = c_vdops;
        vdops[vdops.size(1) * i + 3] = b_vdops;
        vdops[vdops.size(1) * i + 4] = d_vdops;
        vdops[vdops.size(1) * i + 5] = d;
        vdops[vdops.size(1) * i + 6] = c_vdops;
        vdops[vdops.size(1) * i + 7] = d;
        vdops[vdops.size(1) * i + 8] = e_vdops;
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_grad_div(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::array<double, 2U>
    &vdops)
  {
    int dim;
    signed char grad_div_data[9];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts)
    //  [wls, vdops] = wls_var_grad_div(wls, quad_pnts, ws)
    //  [wls, vdops, reslut] = wls_var_grad_div(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each quadrature point
    //                of size n-by-c, where c is 0, 1 or d. In general, each
    //                weight should be the product of the quadrature weight,
    //                Jacobian determinant, and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //                Use n-by-0 for unit weight; n-by-1 indicates a single
    //                weight for all entries at each point; n-by-d will be used
    //                to scale each component of the operator, respectively.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-d).
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-(d^2).
    //                To apply the operator in 2D, use
    //                   [sum(fs(:,1:2) .* vdops(:,1:2), 'all');
    //                    sum(fs(:,1:2) .* vdops(:,3:4), 'all')];
    //                To apply the operator in 3D, use
    //                   [sum(fs(:,1:3) .* vdops(:,1:3), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,4:6), 'all');
    //                    sum(fs(:,1:3) .* vdops(:,7:9), 'all')];
    //     result:    Computed solution of size d-by-1.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_lap, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size;
      int i;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

      //  All components share the same weight, we need to compute Hessian
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size = 1;
        hess_data[0] = 3;
        break;

       case 2:
        grad_div_size = 4;
        hess_data[0] = 4;
        hess_data[1] = 5;
        hess_data[2] = 6;
        hess_data[3] = 0;
        break;

       default:
        grad_div_size = 9;
        for (u1 = 0; u1 < 9; u1++) {
          hess_data[u1] = iv2[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(true, "");

#else //MATLAB_MEX_FILE

      assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (int b_dim{0}; b_dim < nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size, nrows_vdops);
      u0 = nrows_vdops * grad_div_size;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      nOps = grad_div_size;
      i = grad_div_size;
      while ((i > 0) && (hess_data[i - 1] == 0)) {
        nOps--;
        i--;
      }

      //  Summing up rows in the differential operator
      //  Loop through the operators
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                  b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] +
                  ws[ws.size(1) * iPoint] * b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size);
      for (i = 0; i < nrows_vdops; i++) {
        for (j = 0; j < grad_div_size; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k < nOps; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }

      if (dim == 2) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double b_vdops;
          double c_vdops;
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i + 2];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = b_vdops;
          vdops[vdops.size(1) * i + 3] = c_vdops;
        }
      } else if (dim == 3) {
        u1 = b_wls->nrows;
        for (i = 0; i < u1; i++) {
          double b_vdops;
          double c_vdops;
          double d;
          double d_vdops;
          double e_vdops;
          b_vdops = vdops[vdops.size(1) * i + 1];
          c_vdops = vdops[vdops.size(1) * i + 3];
          d_vdops = vdops[vdops.size(1) * i + 2];
          d = vdops[vdops.size(1) * i + 4];
          e_vdops = vdops[vdops.size(1) * i + 5];
          vdops[vdops.size(1) * i + 1] = b_vdops;
          vdops[vdops.size(1) * i + 2] = c_vdops;
          vdops[vdops.size(1) * i + 3] = b_vdops;
          vdops[vdops.size(1) * i + 4] = d_vdops;
          vdops[vdops.size(1) * i + 5] = d;
          vdops[vdops.size(1) * i + 6] = c_vdops;
          vdops[vdops.size(1) * i + 7] = d;
          vdops[vdops.size(1) * i + 8] = e_vdops;
        }
      }
    } else {
      int grad_div_size;
      int iPoint;
      int iWeight;
      int j;
      int lenWs;
      int nDiff;
      int nDims;
      int ncols;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;
      boolean_T flag;

      //  Each component has its own weight, we need to compute all components
      switch (b_wls->us.size(1)) {
       case 1:
        grad_div_size = 1;
        grad_div_data[0] = 3;
        break;

       case 2:
        grad_div_size = 4;
        grad_div_data[0] = 4;
        grad_div_data[1] = 5;
        grad_div_data[2] = 5;
        grad_div_data[3] = 6;
        break;

       default:
        grad_div_size = 9;
        for (u1 = 0; u1 < 9; u1++) {
          grad_div_data[u1] = iv6[u1];
        }
        break;
      }

      //  Compute variational differential operators as weighted sum at quadrature points
      //
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
      //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
      //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
      //
      //  Parameters
      //  ----------
      //     wls:       A d-dimensional WlsObject, including work spaces.
      //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
      //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
      //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
      //                (for Laplacian). In each row, diff_idx may be padded
      //                with zeros for memory preallocation. The operators in a row
      //                will be summed up before applying QR factoriztaion.
      //     ws:        Weight at each quadrature point for each differential operator
      //                (n-by-k), where k is 0, 1 or d. Each weight should be the
      //                product of the quadrature weight, Jacobian determinant,
      //                and value of a weighting function (e.g., a test function
      //                or the derivative of test function). Use empty for unit weight.
      //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
      //                where s is 1 or d for scalar and vector-valued functions.
      //
      //  Returns
      //  -------
      //     wls:       Updated WlsObject
      //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
      //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
      //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
      //
      //  See also
      //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
      //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
      npoints = quad_pnts.size(0) - 1;
      ncols = b_wls->ncols;
      nDims = quad_pnts.size(1);
      nDiff = grad_div_size - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;
      flag = (grad_div_size / lenWs * lenWs == grad_div_size);

      //  Throw error if condition false
      //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

      mxAssert(flag, "");

#else //MATLAB_MEX_FILE

      if (!flag) {
        fprintf(stderr, "Runtime assertion error.\n");
        fflush(stderr);
      }

      assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

      //  scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      for (dim = 0; dim < nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }

      //  compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  force each operator (rhs) to be stored contiguously
      b_wls->vdops.set_size(grad_div_size, nrows_vdops);
      u0 = nrows_vdops * grad_div_size;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      //  Summing up rows in the differential operator
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp <= nDiff; iOp++) {
        int offset;

        //  Skip padded zeros in the differential operator
        offset = (grad_div_data[iOp] - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
          j = b_wls->jpvt[iMonomial] - 1;
          if (ws.size(0) == 0) {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V
                [(offset + iPoint) + b_wls->V.size(1) * j];
            }
          } else {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] =
                b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws
                [(iWeight + ws.size(1) * iPoint) - 1] * b_wls->V[(offset +
                iPoint) + b_wls->V.size(1) * j];
            }
          }
        }

        if (iWeight == lenWs) {
          iWeight = 1;
        } else {
          iWeight++;
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops,
                   grad_div_size);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, grad_div_size, b_wls->work);

      //  Transpose the operators to column major
      vdops.set_size(nrows_vdops, grad_div_size);
      for (int i{0}; i < nrows_vdops; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) *
            j];
        }
      }

      nrows = b_wls->nrows;
      if (b_wls->rweights.size(0) != 0) {
        for (int k{0}; k <= nDiff; k++) {
          for (int iRow{0}; iRow < nrows; iRow++) {
            vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
              b_wls->rweights[iRow];
          }
        }
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_hess(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int hess_size;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];
    boolean_T flag;

    //  Compute variational Hessian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts)
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_hess(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or length of operator (d*(d+1)/2.
    //                Each weight should be the product of the quadrature weight,
    //                Jacobian determinant and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size D-by-size(fs, 2), where D=d*(d+1)/2,
    //                the size of tril(H), in the order of [dx^2; dxdy; dy^2] in
    //                2D and [dx^2; dxdy; dy^2; dxdz; dydz; dz^2] in 3D.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      hess_size = 1;
      hess_data[0] = 3;
      break;

     case 2:
      hess_size = 3;
      hess_data[0] = 4;
      hess_data[1] = 5;
      hess_data[2] = 6;
      break;

     default:
      hess_size = 6;
      for (u1 = 0; u1 < 6; u1++) {
        hess_data[u1] = static_cast<signed char>(u1 + 5);
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = hess_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (hess_size / lenWs * lenWs == hess_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(hess_size, nrows_vdops);
    u0 = nrows_vdops * hess_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, hess_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                hess_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, hess_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(hess_size, fs.size(1));
      u0 = fs.size(1) * hess_size;
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        iOp = 0;
        for (int iDiff{0}; iDiff <= nDiff; iDiff++) {
          for (iRow = 0; iRow <= nrows; iRow++) {
            result[iFunc + result.size(1) * iOp] = result[iFunc + result.size(1)
              * iOp] + fs[iFunc + fs.size(1) * iRow] * vdops[iDiff + vdops.size
              (1) * iRow];
          }

          if (iOp + 1 == result.size(0)) {
            iOp = 0;
          } else {
            iOp++;
          }
        }
      }
    }
  }

  static inline
  void wls_var_hess(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int hess_size;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];

    //  Compute variational Hessian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts)
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_hess(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or length of operator (d*(d+1)/2.
    //                Each weight should be the product of the quadrature weight,
    //                Jacobian determinant and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size D-by-size(fs, 2), where D=d*(d+1)/2,
    //                the size of tril(H), in the order of [dx^2; dxdy; dy^2] in
    //                2D and [dx^2; dxdy; dy^2; dxdz; dydz; dz^2] in 3D.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      hess_size = 1;
      hess_data[0] = 3;
      break;

     case 2:
      hess_size = 3;
      hess_data[0] = 4;
      hess_data[1] = 5;
      hess_data[2] = 6;
      break;

     default:
      hess_size = 6;
      for (u1 = 0; u1 < 6; u1++) {
        hess_data[u1] = static_cast<signed char>(u1 + 5);
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = hess_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(hess_size, nrows_vdops);
    u0 = nrows_vdops * hess_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
            vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
            iPoint) + b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, hess_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                hess_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, hess_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_hess(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    int hess_size;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];
    boolean_T flag;

    //  Compute variational Hessian operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts)
    //  [wls, vdops] = wls_var_hess(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_hess(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1, or length of operator (d*(d+1)/2.
    //                Each weight should be the product of the quadrature weight,
    //                Jacobian determinant and value of a weighting function
    //                (eg., a test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator, use vdops' * fs.
    //     result:    Computed solution of size D-by-size(fs, 2), where D=d*(d+1)/2,
    //                the size of tril(H), in the order of [dx^2; dxdy; dy^2] in
    //                2D and [dx^2; dxdy; dy^2; dxdz; dydz; dz^2] in 3D.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    switch (b_wls->us.size(1)) {
     case 1:
      hess_size = 1;
      hess_data[0] = 3;
      break;

     case 2:
      hess_size = 3;
      hess_data[0] = 4;
      hess_data[1] = 5;
      hess_data[2] = 6;
      break;

     default:
      hess_size = 6;
      for (u1 = 0; u1 < 6; u1++) {
        hess_data[u1] = static_cast<signed char>(u1 + 5);
      }
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    nDiff = hess_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (hess_size / lenWs * lenWs == hess_size);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(hess_size, nrows_vdops);
    u0 = nrows_vdops * hess_size;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
              iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial + b_wls->vdops.size(1) * iOp] = b_wls->
              vdops[iMonomial + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
              ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        }
      }

      if (iWeight == lenWs) {
        iWeight = 1;
      } else {
        iWeight++;
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, hess_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                hess_size, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops, hess_size);
    for (int i{0}; i < nrows_vdops; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * i] = b_wls->vdops[i + b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (int iRow{0}; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_lap(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                    coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                    &vdops, ::coder::array<double, 2U> &result)
  {
    int iPoint;
    int iRow;
    int lap_size_idx_1;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];
    boolean_T flag;

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_lap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      lap_size_idx_1 = 1;
      lap_data[0] = 3;
      break;

     case 2:
      lap_size_idx_1 = 2;
      lap_data[0] = 4;
      lap_data[1] = 5;
      break;

     default:
      lap_size_idx_1 = 3;
      lap_data[0] = 5;
      lap_data[1] = 6;
      lap_data[2] = 7;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u1 = 0; u1 < nrows_vdops; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset
              + iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
              iPoint] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows - 1;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow <= nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if ((fs.size(0) == 0) || (fs.size(1) == 0)) {
      result.set_size(0, 0);
    } else {
      result.set_size(1, fs.size(1));
      u0 = fs.size(1);
      for (u1 = 0; u1 < u0; u1++) {
        result[u1] = 0.0;
      }

      //  Compute solution
      u1 = fs.size(1);
      for (int iFunc{0}; iFunc < u1; iFunc++) {
        for (iRow = 0; iRow <= nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_lap(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, ::coder::array<double, 1U> &vdops, int
                    result_size[2])
  {
    int iPoint;
    int lap_size_idx_1;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_lap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      lap_size_idx_1 = 1;
      lap_data[0] = 3;
      break;

     case 2:
      lap_size_idx_1 = 2;
      lap_data[0] = 4;
      lap_data[1] = 5;
      break;

     default:
      lap_size_idx_1 = 3;
      lap_data[0] = 5;
      lap_data[1] = 6;
      lap_data[2] = 7;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(true, "");

#else //MATLAB_MEX_FILE

    assert(true);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset +
            iPoint) + b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    result_size[1] = 0;
    result_size[0] = 0;
  }

  static inline
  void wls_var_lap(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                    array<double, 1U> &vdops, int result_size[2])
  {
    int iPoint;
    int lap_size_idx_1;
    int lenWs;
    int nDims;
    int ncols;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];
    boolean_T flag;

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts)
    //  [wls, vdops] = wls_var_lap(wls, quad_pnts, ws)
    //  [wls, vdops, result] = wls_var_lap(wls, quad_pnts, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d)
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0 or 1. Each weight should
    //                be the product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (m-by-s),
    //                where s is the number of scalar functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-d. To apply the
    //                operator to compute vector Laplacian, use vdops' * fs.
    //     result:    Computed solution of size 1-by-size(fs, 2).
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_vec_hess, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap, wls_var_kernel
    //  The operators are row vectors, so they will be summed up before solve
    switch (b_wls->us.size(1)) {
     case 1:
      lap_size_idx_1 = 1;
      lap_data[0] = 3;
      break;

     case 2:
      lap_size_idx_1 = 2;
      lap_data[0] = 4;
      lap_data[1] = 5;
      break;

     default:
      lap_size_idx_1 = 3;
      lap_data[0] = 5;
      lap_data[1] = 6;
      lap_data[2] = 7;
      break;
    }

    //  Compute variational differential operators as weighted sum at quadrature points
    //
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx)
    //  [wls, vdops] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws)
    //  [wls, vdops, result] = wls_var_kernel(wls, quad_pnts, order, diff_idx, ws, fs)
    //
    //  Parameters
    //  ----------
    //     wls:       A d-dimensional WlsObject, including work spaces.
    //     quad_pnts: Local coordinates of the quadrature points (n-by-d).
    //     order:     order of confluent Vandermonde matrix. Use -1 to indicate div.
    //     diff_idx:  Differential operators of size nDiff-by-1 or 1-by-njDiff
    //                (for Laplacian). In each row, diff_idx may be padded
    //                with zeros for memory preallocation. The operators in a row
    //                will be summed up before applying QR factoriztaion.
    //     ws:        Weight at each quadrature point for each differential operator
    //                (n-by-k), where k is 0, 1 or d. Each weight should be the
    //                product of the quadrature weight, Jacobian determinant,
    //                and value of a weighting function (e.g., a test function
    //                or the derivative of test function). Use empty for unit weight.
    //     fs:        Values corresponding to rows in CVM in wls object (n-by-s),
    //                where s is 1 or d for scalar and vector-valued functions.
    //
    //  Returns
    //  -------
    //     wls:       Updated WlsObject
    //     vdops:     The computed operator of size wls.nrows-by-size(diff_idx, 1).
    //     result:    Computed solution. Its size is size(diff_idx, 1)-by-size(fs, 2)
    //                if order is positive or size(diff_idx, 1)/d-by-size(fs, 2) otherwise.
    //
    //  See also
    //     wls_init, wls_var_diff, wls_var_func, wls_var_grad, wls_var_div, wls_var_curl,
    //     wls_var_hess, wls_var_lap, wls_var_grad_div, wls_var_curl_curl, wls_var_bilap
    npoints = quad_pnts.size(0) - 1;
    ncols = b_wls->ncols;
    nDims = quad_pnts.size(1);
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;
    flag = (1 / lenWs * lenWs == 1);

    //  Throw error if condition false
    //  C++
#ifndef NDEBUG
#ifdef MATLAB_MEX_FILE

    mxAssert(flag, "");

#else //MATLAB_MEX_FILE

    if (!flag) {
      fprintf(stderr, "Runtime assertion error.\n");
      fflush(stderr);
    }

    assert(flag);

#endif //MATLAB_MEX_FILE
#endif //NDEBUG

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    for (int dim{0}; dim < nDims; dim++) {
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
          quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
      }
    }

    //  compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  force each operator (rhs) to be stored contiguously
    b_wls->vdops.set_size(1, nrows_vdops);
    for (u0 = 0; u0 < nrows_vdops; u0++) {
      b_wls->vdops[u0] = 0.0;
    }

    //  Omit zeros in the diff operators
    //  Summing up rows in the differential operator
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      //  Sum up monomials weighted by weights for each component
      for (int iMonomial{0}; iMonomial < ncols; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + b_wls->V[(offset
              + iPoint) + b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->vdops[iMonomial] = b_wls->vdops[iMonomial] + ws[ws.size(1) *
              iPoint] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);

    //  Transpose the operators to column major
    vdops.set_size(nrows_vdops);
    for (int i{0}; i < nrows_vdops; i++) {
      vdops[i] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    result_size[1] = 0;
    result_size[0] = 0;
  }
}

// End of code generation (wls_internal.cpp)
