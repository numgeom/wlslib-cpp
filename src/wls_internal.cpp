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
    rowoffset, int coloffset, int m, int n, ::coder::array<double, 2U> &QR, ::
    coder::array<int, 1U> &p, ::coder::array<double, 1U> &work, const ::coder::
    array<unsigned char, 2U> &dag);
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
        b_npoints = npoints - 1;

        //  Generate (confluent) Vandermonde matrix in 1D.
        if (npoints == 0) {
          b_npoints = us.size(0) - 1;
        }

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
    switch (us_size[1]) {
     case 1:
      {
        int i;

        //  Generate (confluent) Vandermonde matrix in 1D.
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

        //  Generate (confluent) Vandermonde matrix in 1D.
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
    int b_i;

    //  Build a dag for Vandermonde matrix in 1D.
    for (int i{0}; i < degree; i++) {
      b_i = dag.size(0);
      dag[i % b_i * dag.size(1) + i / b_i] = 1U;
    }

    int i1;
    b_i = dag.size(0);
    dag[degree % b_i * dag.size(1) + degree / b_i] = 0U;

    //  a leaf has no child
    b_i = dag.size(1) * dag.size(0) - 1;
    i1 = dag.size(0);
    dag[b_i % i1 * dag.size(1) + b_i / i1] = static_cast<unsigned char>(degree +
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
    ncols = (degree + 1) * (degree + 2) / 2;

    //  Allocate storage for V
    if ((V.size(1) != 1) || (V.size(0) != ncols)) {
      V.set_size(ncols, 1);
    }

    //  compute 0th order generalized Vandermonde matrix
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

    //  Generate generalized/confluent Vandermonde matrix in 2D.
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
      for (iPnt = 0; iPnt < npoints; iPnt++) {
        V[stride + iPnt] = 0.0;
        V[(stride + iPnt) + V.size(1)] = V[iPnt];
        V[(stride + iPnt) + V.size(1) * 2] = 0.0;
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
              offset = 3 * us.size(0);

              //  Compute derivative with respect to u
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
              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);

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

    //  Generate generalized/confluent Vandermonde matrix in 2D.
    if (npoints == 0) {
      npoints = us.size(0);
    }

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
              offset = 3 * us.size(0);

              //  Compute derivative with respect to u
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
              offset = 5 * us.size(0);
              offset_prev = 3 * us.size(0);

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
    ncols = (degree + 1) * (degree + 2) * (degree + 3) / 6;

    //  Allocate storage for V
    if ((V.size(1) != 1) || (V.size(0) != ncols)) {
      V.set_size(ncols, 1);
    }

    //  compute 0th order generalized Vandermonde matrix
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

    //  Generate generalized/confluent Vandermonde matrix in 3D.
    if (npoints == 0) {
      npoints = us.size(0);
    }

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

          //  Compute order-1 CVM row blocks from order-0 GVM.
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

    //  Generate generalized/confluent Vandermonde matrix in 3D.
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

    if (order != 0) {
      //     %% compute higher order confluent Vandermonde matrix blocks incrementally
      switch (order) {
       case 1:
        {
          int balance;
          int kdegree;
          int offset;

          //  Compute order-1 CVM row blocks from order-0 GVM.
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[stride + iPnt] = 0.0;
            V[(stride + iPnt) + V.size(1)] = V[iPnt];
            V[(stride + iPnt) + V.size(1) * 2] = 0.0;
            V[(stride + iPnt) + V.size(1) * 3] = 0.0;
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
                V[(stride + iPnt) + V.size(1) * (c + 1)] = V[(stride + iPnt) +
                  V.size(1) * ((c - d) - deg)] * us[us.size(1) * iPnt + 2];
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
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - nTermsInLayer)] * us[us.size(1) * iPnt + 1];
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
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - balance)] * us[us.size(1) * iPnt + 2];
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
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[stride + iPnt] = 0.0;
            V[(stride + iPnt) + V.size(1)] = V[iPnt];
            V[(stride + iPnt) + V.size(1) * 2] = 0.0;
            V[(stride + iPnt) + V.size(1) * 3] = 0.0;
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
                V[(stride + iPnt) + V.size(1) * (c + 1)] = V[(stride + iPnt) +
                  V.size(1) * ((c - d) - deg)] * us[us.size(1) * iPnt + 2];
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
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - nTermsInLayer)] * us[us.size(1) * iPnt + 1];
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
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - balance)] * us[us.size(1) * iPnt + 2];
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
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[stride + iPnt] = 0.0;
            V[(stride + iPnt) + V.size(1)] = V[iPnt];
            V[(stride + iPnt) + V.size(1) * 2] = 0.0;
            V[(stride + iPnt) + V.size(1) * 3] = 0.0;
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
                V[(stride + iPnt) + V.size(1) * (c + 1)] = V[(stride + iPnt) +
                  V.size(1) * ((c - d) - deg)] * us[us.size(1) * iPnt + 2];
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
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - nTermsInLayer)] * us[us.size(1) * iPnt + 1];
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
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - balance)] * us[us.size(1) * iPnt + 2];
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
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[stride + iPnt] = 0.0;
            V[(stride + iPnt) + V.size(1)] = V[iPnt];
            V[(stride + iPnt) + V.size(1) * 2] = 0.0;
            V[(stride + iPnt) + V.size(1) * 3] = 0.0;
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
                V[(stride + iPnt) + V.size(1) * (c + 1)] = V[(stride + iPnt) +
                  V.size(1) * ((c - d) - deg)] * us[us.size(1) * iPnt + 2];
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
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - nTermsInLayer)] * us[us.size(1) * iPnt + 1];
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
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - balance)] * us[us.size(1) * iPnt + 2];
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

          //  Compute order-1 CVM row blocks from order-0 GVM.
          for (iPnt = 0; iPnt < npoints; iPnt++) {
            V[stride + iPnt] = 0.0;
            V[(stride + iPnt) + V.size(1)] = V[iPnt];
            V[(stride + iPnt) + V.size(1) * 2] = 0.0;
            V[(stride + iPnt) + V.size(1) * 3] = 0.0;
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
                V[(stride + iPnt) + V.size(1) * (c + 1)] = V[(stride + iPnt) +
                  V.size(1) * ((c - d) - deg)] * us[us.size(1) * iPnt + 2];
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
              cornerTriangle = (cornerTriangle + p) + degree;
              counterBottomRow = 1;

              // counter for the bottom row to be subtracted later
              for (kdegree = 0; kdegree < deg; kdegree++) {
                for (iPnt = 0; iPnt < npoints; iPnt++) {
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - nTermsInLayer)] * us[us.size(1) * iPnt + 1];
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
                  V[(stride + iPnt) + V.size(1) * c] = V[(stride + iPnt) +
                    V.size(1) * (c - balance)] * us[us.size(1) * iPnt + 2];
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
    rowoffset, int coloffset, int m, int n, ::coder::array<double, 2U> &QR, ::
    coder::array<int, 1U> &p, ::coder::array<double, 1U> &work, const ::coder::
    array<unsigned char, 2U> &dag)
  {
    int rank;

    //  rrqr_factor  Compute rank-revealing QR with column pivoting
    if (m == 0) {
      m = A.size(1) - rowoffset;
    }

    if (n == 0) {
      n = A.size(0) - coloffset;
    }

    //  Preallocate output arguments
    p[0] = 0;
    if ((0 == n) || ((dag.size(0) == 0) || (dag.size(1) == 0))) {
      int i;
      int i1;
      int i2;

      //  Note: A and Q are always stored in column major
      i = coloffset * A.size(1) + rowoffset;
      i1 = A.size(1);
      i2 = A.size(0);
      rank = wls::rrqr_factor_nodag(&A[i % i2 * i1 + i / i2], thres, m, n, &QR[0],
        &(p.data())[0], &(work.data())[0], work.size(0), A.size(1));
    } else {
      int i;
      int i1;
      int i2;
      int i3;
      int i4;
      int i5;

      //  Note: A, Q, and dag are always stored in column major
      i = coloffset * A.size(1) + rowoffset;
      i1 = A.size(1);
      i2 = A.size(0);
      i3 = static_cast<int>(static_cast<double>(coloffset) * static_cast<double>
                            (dag.size(1)));
      i4 = dag.size(1);
      i5 = dag.size(0);
      rank = wls::rrqr_factor(&A[i % i2 * i1 + i / i2], thres, m, n, &QR[0],
        &(p.data())[0], &(work.data())[0], work.size(0), A.size(1), &dag[i3 % i5
        * i4 + i3 / i5], dag.size(1), dag.size(0) - coloffset);
    }

    return rank;
  }

  static inline
  void rrqr_qmulti(const ::coder::array<double, 2U> &QR, int m, int n,
    int rank, ::coder::array<double, 2U> &bs, int nrhs, ::coder::array<double,
    1U> &work)
  {
    int i;
    int stride_bs;
    int wsize;

    //  Perform Q*bs, where Q is stored implicitly in QR
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

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Resize work space if needed
    wsize = wls::query_work_size(m, n);
    if (work.size(0) < wsize) {
      work.set_size(wsize);
    }

    //  zero out extra rows in bs to avoid errors in LAPACK
    i = n + 1;
    for (int b_i{i}; b_i <= m; b_i++) {
      for (int j{0}; j < nrhs; j++) {
        bs[(b_i + bs.size(1) * j) - 1] = 0.0;
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
    //  Perform forward substitution to compute bs=R'\bs, where R is stored in QR
    if (n == 0) {
      n = QR.size(0) - 1;
    }

    if (rank == 0) {
      rank = n;
    }

    if (nrhs == 0) {
      nrhs = bs.size(0);
    }

    //  Obtain stride
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

    //  Weights based on Buhmann's radial basis function
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

    //  WLS-ENO weights based on function values
    if (ws.size(0) == 0) {
      ws.set_size(npoints);
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
    epsilon = 0.01;
    alpha = std::abs(b_degree) / 2.0;
    if (h2bar_tmp != 0.0) {
      epsilon = h2bar_tmp;
    }

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

    //  Weights based on inverse distance
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
    int unnamed_idx_0;

    //  Reinitialize the buffers of WlsObject
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
        b_wls->dag.set_size(ncols + 1, dim);
        unnamed_idx_0 = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_degree = b_wls->dag.size(0);
        b_wls->dag[unnamed_idx_0 % b_degree * b_wls->dag.size(1) + unnamed_idx_0
          / b_degree] = MAX_uint8_T;
      }
    } else {
      b_wls->dag.set_size(0, dim);
    }

    b_wls->jpvt.set_size(ncols);

    //  V is always full, but QR has one fewer row and column in interp0 mode
    b_wls->V.set_size(ncols, nrows);
    unnamed_idx_0 = (ncols - b_wls->interp0) + 1;
    b_wls->QR.set_size(unnamed_idx_0, nrows);
    b_wls->rank = 0;

    //  work space
    b_degree = ncols << 2;
    unnamed_idx_0 = ncols + 1;
    if (nrows >= unnamed_idx_0) {
      unnamed_idx_0 = nrows;
    }

    if (b_degree < 4160) {
      b_degree = 4160;
    }

    b_wls->work.set_size((unnamed_idx_0 << 5) + b_degree);
  }

  static inline
  void wls_func1(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, const
                 ::coder::array<double, 2U> &fs, int npoints, ::coder::array<
                 double, 2U> &vdops, ::coder::array<double, 2U> &result)
  {
    int i;
    int iPoint;
    int iRow;
    int j;
    int nDims;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    nDims = pnts.size(1) - 1;

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((npoints + 3) / 4) << 2, pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint < npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (pnts[dim + pnts.size(1)
            * iPoint] - b_wls->origin.data[dim]) * b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint < npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
            iPoint] * b_wls->hs_inv.data[dim];
        }
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
    u1 = nrows_vdops - b_wls->interp0;
    b_wls->vdops.set_size(npoints, u1);

    //  Extract vopts from Vandermonde matrix
    u1 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint < npoints; iPoint++) {
        b_wls->vdops[(iMonomial + b_wls->interp0) + b_wls->vdops.size(1) *
          iPoint] = b_wls->V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, npoints);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                npoints, b_wls->work);
    vdops.set_size(nrows_vdops, npoints);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j < npoints; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k < npoints; k++) {
        for (iRow = 0; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j < npoints; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 1.0 - s;
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
          for (iRow = 0; iRow < nrows; iRow++) {
            result[iFunc + result.size(1) * iPnt] = result[iFunc + result.size(1)
              * iPnt] + fs[iFunc + fs.size(1) * iRow] * vdops[iPnt + vdops.size
              (1) * iRow];
          }
        }
      }
    }
  }

  static inline
  void wls_func2(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, const
                 ::coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                 &vdops, ::coder::array<double, 2U> &result)
  {
    int i;
    int iPoint;
    int iRow;
    int j;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    npoints = pnts.size(0) - 1;
    nDims = pnts.size(1) - 1;

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((pnts.size(0) + 3) / 4) << 2, pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (pnts[dim + pnts.size(1)
            * iPoint] - b_wls->origin.data[dim]) * b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
            iPoint] * b_wls->hs_inv.data[dim];
        }
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
    u1 = nrows_vdops - b_wls->interp0;
    b_wls->vdops.set_size(pnts.size(0), u1);

    //  Extract vopts from Vandermonde matrix
    u1 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->vdops[(iMonomial + b_wls->interp0) + b_wls->vdops.size(1) *
          iPoint] = b_wls->V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, pnts.size(0));
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                pnts.size(0), b_wls->work);
    vdops.set_size(nrows_vdops, pnts.size(0));

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= npoints; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= npoints; k++) {
        for (iRow = 0; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= npoints; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 1.0 - s;
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
          for (iRow = 0; iRow < nrows; iRow++) {
            result[iFunc + result.size(1) * iPnt] = result[iFunc + result.size(1)
              * iPnt] + fs[iFunc + fs.size(1) * iRow] * vdops[iPnt + vdops.size
              (1) * iRow];
          }
        }
      }
    }
  }

  static inline
  void wls_func3(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts, ::
                 coder::array<double, 2U> &vdops)
  {
    int i;
    int iPoint;
    int j;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute WLS-fitting at one or more points.
    npoints = pnts.size(0) - 1;
    nDims = pnts.size(1) - 1;

    //  scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((pnts.size(0) + 3) / 4) << 2, pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (pnts[dim + pnts.size(1)
            * iPoint] - b_wls->origin.data[dim]) * b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = pnts[dim + pnts.size(1) *
            iPoint] * b_wls->hs_inv.data[dim];
        }
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
    u0 = nrows_vdops - b_wls->interp0;
    b_wls->vdops.set_size(pnts.size(0), u0);

    //  Extract vopts from Vandermonde matrix
    u0 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u0; iMonomial++) {
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        b_wls->vdops[(iMonomial + b_wls->interp0) + b_wls->vdops.size(1) *
          iPoint] = b_wls->V[iPoint + b_wls->V.size(1) * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, pnts.size(0));
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                pnts.size(0), b_wls->work);
    vdops.set_size(nrows_vdops, pnts.size(0));

    //  Transpose the operator for row-major
    u0 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u0; i++) {
      for (j = 0; j <= npoints; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= npoints; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 1.0 - s;
      }
    }
  }

  static inline
  void wls_init1(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, int interp0,
                 boolean_T use_dag, int npoints)
  {
    int dim;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    b_wls->interp0 = (interp0 != 0);
    interp0 = b_wls->interp0;
    b_wls->use_dag = use_dag;
    if (npoints <= 0) {
      npoints = us.size(0);
    }

    //  Resize buffers
    wls_resize(b_wls, us.size(1), npoints, degree, order, use_dag);
    if (npoints != 0) {
      double maxx;
      double maxx_inv;
      double thres;
      int b_b;
      int b_i;
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i < npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = ((us.size(1) <= 0) || (us.size(0) <= 0));
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i < npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i < npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i < npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i < npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i < npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i < npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init2(WlsObject *b_wls, const ::coder::array<double, 2U> &us)
  {
    int degree;
    int dim;
    int interp0;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;

    //  Default is to use unit weight
    degree = b_wls->degree;
    order = b_wls->order;
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              int i1;
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                i1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > i1) {
                    b_us -= i1;
                  }
                }
              }

              i1 = b_wls->us.size(0);
              b_wls->us[b_i % i1 * b_wls->us.size(1) + b_i / i1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              int i1;
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                i1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > i1) {
                    b_us -= i1;
                  }
                }
              }

              i1 = b_wls->us.size(0);
              b_wls->us[b_i % i1 * b_wls->us.size(1) + b_i / i1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
      }

      //  Compute point-wise weights
      if (order == 0) {
        //  Unit weights
        b_wls->rweights.set_size(0);
      } else {
        b_wls->rweights.set_size(b_wls->V.size(1));

        //  unit weights
        b_us = b_wls->rweights.size(0);
        b_wls->rweights.set_size(b_us);
        for (i = 0; i < b_us; i++) {
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
      if ((degree > 1) && (degree < 7)) {
        thres = dv[degree - 1];
      } else {
        thres = 1.0E+8;
      }

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init3(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight)
  {
    int degree;
    int dim;
    int interp0;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    degree = b_wls->degree;
    order = b_wls->order;
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init4(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree)
  {
    int dim;
    int interp0;
    int npoints;
    int order;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    order = b_wls->order;
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init5(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order)
  {
    int dim;
    int interp0;
    int npoints;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init6(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, int interp0)
  {
    int dim;
    int npoints;
    boolean_T use_dag;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    b_wls->interp0 = (interp0 != 0);
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_init7(WlsObject *b_wls, const ::coder::array<double, 2U> &us, const
                 WlsWeight *weight, int degree, int order, int interp0,
                 boolean_T use_dag)
  {
    int dim;
    int npoints;

    //  Initialize WlsObject in 1D, 2D, or 3D.
    dim = us.size(1) - 1;
    b_wls->interp0 = (interp0 != 0);
    interp0 = b_wls->interp0;
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
      int b_us;
      int i;
      int j;
      int ncols;
      int nrblks;
      int src;
      int trg;
      int u1;

      //  Recompute DAG if use_dag and its signature does not match
      if (use_dag) {
        i = b_wls->dag.size(1) * b_wls->dag.size(0) - 1;
        b_us = b_wls->dag.size(0);
        if (b_wls->dag[i % b_us * b_wls->dag.size(1) + i / b_us] != degree + 127)
        {
          //  Wrapper function for building DAG in nD.
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

      if (b_wls->interp0 != 0) {
        //  Make the first node the origin in interp0 mode
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = us[0];
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us] -
                us[0];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = us[0];
          b_wls->origin.data[1] = us[1];
          b_wls->origin.data[2] = us[2];
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i] - us[0];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1] -
              us[1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2] -
              us[2];
          }
          break;
        }
      } else {
        switch (us.size(1)) {
         case 1:
          {
            boolean_T b;
            boolean_T b1;
            b_wls->origin.size[1] = 1;
            b_wls->origin.size[0] = 1;
            b_wls->origin.data[0] = 0.0;
            b = true;
            b1 = (us.size(1) <= 0);
            i = us.size(1) * us.size(0);
            b_us = 0;
            for (b_i = 0; b_i <= npoints; b_i++) {
              if (b1 || (b_i >= i)) {
                b_us = 0;
                b = true;
              } else if (b) {
                b = false;
                b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
              } else {
                u1 = us.size(1) * us.size(0) - 1;
                if (b_us > MAX_int32_T - us.size(1)) {
                  b_us = b_i % us.size(0) * us.size(1) + b_i / us.size(0);
                } else {
                  b_us += us.size(1);
                  if (b_us > u1) {
                    b_us -= u1;
                  }
                }
              }

              u1 = b_wls->us.size(0);
              b_wls->us[b_i % u1 * b_wls->us.size(1) + b_i / u1] = us[b_us];
            }
          }
          break;

         case 2:
          b_wls->origin.size[1] = 2;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
          }
          break;

         default:
          b_wls->origin.size[1] = 3;
          b_wls->origin.size[0] = 1;
          b_wls->origin.data[0] = 0.0;
          b_wls->origin.data[1] = 0.0;
          b_wls->origin.data[2] = 0.0;
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = us[us.size(1) * b_i];
            b_wls->us[b_wls->us.size(1) * b_i + 1] = us[us.size(1) * b_i + 1];
            b_wls->us[b_wls->us.size(1) * b_i + 2] = us[us.size(1) * b_i + 2];
          }
          break;
        }
      }

      //  Scale us to be between -1 and 1
      maxx = 0.0;
      switch (us.size(1)) {
       case 1:
        for (b_i = 0; b_i <= npoints; b_i++) {
          i = b_wls->us.size(0);
          maxx = std::fmax(maxx, std::abs(b_wls->us[b_i % i * b_wls->us.size(1)
            + b_i / i]));
        }
        break;

       case 2:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::abs(b_wls->us[b_wls->us.size(1) *
            b_i]), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 1])));
        }
        break;

       default:
        for (b_i = 0; b_i <= npoints; b_i++) {
          maxx = std::fmax(maxx, std::fmax(std::fmax(std::abs(b_wls->us
            [b_wls->us.size(1) * b_i]), std::abs(b_wls->us[b_wls->us.size(1) *
            b_i + 1])), std::abs(b_wls->us[b_wls->us.size(1) * b_i + 2])));
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

      //  scale wls.us
      if (maxx_inv != 1.0) {
        switch (us.size(1)) {
         case 1:
          for (b_i = 0; b_i <= npoints; b_i++) {
            i = b_wls->us.size(0);
            b_us = b_wls->us.size(0);
            b_wls->us[b_i % i * b_wls->us.size(1) + b_i / i] = b_wls->us[b_i %
              b_us * b_wls->us.size(1) + b_i / b_us] * maxx_inv;
          }
          break;

         case 2:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
          }
          break;

         default:
          for (b_i = 0; b_i <= npoints; b_i++) {
            b_wls->us[b_wls->us.size(1) * b_i] = b_wls->us[b_wls->us.size(1) *
              b_i] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 1] = b_wls->us[b_wls->us.size(1)
              * b_i + 1] * maxx_inv;
            b_wls->us[b_wls->us.size(1) * b_i + 2] = b_wls->us[b_wls->us.size(1)
              * b_i + 2] * maxx_inv;
          }
          break;
        }
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
          b_us = b_wls->rweights.size(0);
          b_wls->rweights.set_size(b_us);
          for (i = 0; i < b_us; i++) {
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
      b_us = weight->omit_rows.size(0);
      u1 = b_wls->nrows;
      if (b_us <= u1) {
        u1 = b_us;
      }

      for (b_i = 0; b_i < u1; b_i++) {
        if (weight->omit_rows[b_i]) {
          b_us = b_wls->V.size(0);
          for (i = 0; i < b_us; i++) {
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

      //  In interp0 mode, we trim off the first row and first column.
      b_wls->rank = rrqr_factor(b_wls->V, thres, interp0, interp0, b_wls->nrows
        - interp0, ncols - interp0, b_wls->QR, b_wls->jpvt, b_wls->work,
        b_wls->dag);
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
  void wls_var_bilap1(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                      coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                      &vdops, ::coder::array<double, 2U> &result)
  {
    int bilap_size_idx_1;
    int i;
    int iPoint;
    int iRow;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
              b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
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
        for (iRow = 0; iRow < nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_bilap2(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, ::coder::array<double, 1U> &vdops)
  {
    int bilap_size_idx_1;
    int i;
    int iPoint;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
            b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
    }
  }

  static inline
  void wls_var_bilap3(WlsObject *b_wls, const ::coder::array<double, 2U>
                      &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                      array<double, 1U> &vdops)
  {
    int bilap_size_idx_1;
    int i;
    int iPoint;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char bilap_data[9];

    //  Compute variational (vector) bi-Laplacian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -4,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < bilap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (bilap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
              b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
    }
  }

  static inline
  void wls_var_curl1(WlsObject *b_wls, const ::coder::array<double, 2U>
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
    if (ws.size(1) <= 1) {
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;

      //  Compute variational differential operators as weighted sum at quadrature points
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
              b_wls->hs_inv.data[dim];
          }
        }
      } else {
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = b_wls->izero_ + 9;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv4[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, 9);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j < 9; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + 9 * (i - 1)];
          }

          vdops[j] = static_cast<double>(iv4[j] == 1) - s;
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
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;

      //  Compute variational differential operators as weighted sum at quadrature points
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
              b_wls->hs_inv.data[dim];
          }
        }
      } else {
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = b_wls->izero_ + 9;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if (ws.size(0) == 0) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) *
                  iPoint) - 1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) *
                  j];
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
      vdops.set_size(nrows_vdops, 9);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j < 9; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + 9 * (i - 1)];
          }

          vdops[j] = static_cast<double>(iv5[j] == 1) - s;
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
  void wls_var_curl2(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int i;
    int iPoint;
    int j;
    int nDims;
    int nOps;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;

    //  Variational curl operators as weighted sum at quadrature points in 3D
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 9;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
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
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = iv4[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        u1 = b_wls->ncols - b_wls->interp0;
        for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);
    vdops.set_size(nrows_vdops, 9);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j < 9; j++) {
        vdops[j + 9 * (i + b_wls->interp0)] = b_wls->vdops[i + b_wls->vdops.size
          (1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j < 9; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + 9 * (i - 1)];
        }

        vdops[j] = static_cast<double>(iv4[j] == 1) - s;
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
  void wls_var_curl3(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    //  Variational curl operators as weighted sum at quadrature points in 3D
    if (ws.size(1) <= 1) {
      int i;
      int iPoint;
      int j;
      int nDims;
      int nOps;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

      //  Compute variational differential operators as weighted sum at quadrature points
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
              b_wls->hs_inv.data[dim];
          }
        }
      } else {
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = b_wls->izero_ + 9;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = iv4[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, 9);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j < 9; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + 9 * (i - 1)];
          }

          vdops[j] = static_cast<double>(iv4[j] == 1) - s;
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
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

      //  Compute variational differential operators as weighted sum at quadrature points
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
              b_wls->hs_inv.data[dim];
          }
        }
      } else {
        for (int dim{0}; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = b_wls->izero_ + 9;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if (ws.size(0) == 0) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) *
                  iPoint) - 1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) *
                  j];
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
      vdops.set_size(nrows_vdops, 9);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j < 9; j++) {
          vdops[j + 9 * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j < 9; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + 9 * (i - 1)];
          }

          vdops[j] = static_cast<double>(iv5[j] == 1) - s;
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
  void wls_var_curl_curl1(WlsObject *b_wls, const ::coder::array<double, 2U>
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
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size_idx_1;
      int iPoint;
      int j;
      int nDiff;
      int nDims;
      int nOps;
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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size_idx_1 - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size_idx_1 + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, grad_div_size_idx_1);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
      int nDiff;
      int nDims;
      int nOps;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;

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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size_idx_0 - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size_idx_0 + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
            u1 = b_wls->ncols - b_wls->interp0;
            for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
              j = b_wls->jpvt[iMonomial] - 1;
              if (ws.size(0) == 0) {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  u0 = iMonomial + b_wls->interp0;
                  b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->
                    vdops[u0 + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
                    iPoint) + b_wls->V.size(1) * j];
                }
              } else {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  u0 = iMonomial + b_wls->interp0;
                  b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->
                    vdops[u0 + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
                    ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
                    b_wls->V.size(1) * j];
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
      vdops.set_size(nrows_vdops, grad_div_size_idx_0);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        boolean_T b;

        //  In interp0 mode, we set the first entry based on partition of unity
        b = true;
        u1 = grad_div_size_idx_0 * grad_div_size_idx_1;
        u0 = 0;
        for (j = 0; j <= nDiff; j++) {
          double s;
          if (j >= u1) {
            u0 = 0;
            b = true;
          } else if (b) {
            b = false;
            u0 = j % grad_div_size_idx_0 * grad_div_size_idx_1 + j /
              grad_div_size_idx_0;
          } else if (u0 > MAX_int32_T - grad_div_size_idx_1) {
            u0 = j % grad_div_size_idx_0 * grad_div_size_idx_1 + j /
              grad_div_size_idx_0;
          } else {
            u0 += grad_div_size_idx_1;
            if (u0 > u1 - 1) {
              u0 = (u0 - u1) + 1;
            }
          }

          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(grad_div_data[u0] == 1) - s;
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
  void wls_var_curl_curl2(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int dim;
    int hess_size;
    int i;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int nOps;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = hess_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int b_dim{0}; b_dim <= nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
            b_wls->hs_inv.data[b_dim];
        }
      }
    } else {
      for (int b_dim{0}; b_dim <= nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = hess_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
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
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = hess_data[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        u1 = b_wls->ncols - b_wls->interp0;
        for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);
    vdops.set_size(nrows_vdops, hess_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
  void wls_var_curl_curl3(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::array<double, 2U>
    &vdops)
  {
    int dim;
    signed char grad_div_data[18];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size_idx_1;
      int i;
      int iPoint;
      int j;
      int nDiff;
      int nDims;
      int nOps;
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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size_idx_1 - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size_idx_1 + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, grad_div_size_idx_1);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
      int nDiff;
      int nDims;
      int nOps;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size_idx_0 - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size_idx_0 + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
            u1 = b_wls->ncols - b_wls->interp0;
            for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
              j = b_wls->jpvt[iMonomial] - 1;
              if (ws.size(0) == 0) {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  u0 = iMonomial + b_wls->interp0;
                  b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->
                    vdops[u0 + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset +
                    iPoint) + b_wls->V.size(1) * j];
                }
              } else {
                for (iPoint = 0; iPoint <= npoints; iPoint++) {
                  u0 = iMonomial + b_wls->interp0;
                  b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->
                    vdops[u0 + b_wls->vdops.size(1) * iOp] + ws[(iWeight +
                    ws.size(1) * iPoint) - 1] * b_wls->V[(offset + iPoint) +
                    b_wls->V.size(1) * j];
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
      vdops.set_size(nrows_vdops, grad_div_size_idx_0);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        boolean_T b;

        //  In interp0 mode, we set the first entry based on partition of unity
        b = true;
        u1 = grad_div_size_idx_0 * grad_div_size_idx_1;
        u0 = 0;
        for (j = 0; j <= nDiff; j++) {
          double s;
          if (j >= u1) {
            u0 = 0;
            b = true;
          } else if (b) {
            b = false;
            u0 = j % grad_div_size_idx_0 * grad_div_size_idx_1 + j /
              grad_div_size_idx_0;
          } else if (u0 > MAX_int32_T - grad_div_size_idx_1) {
            u0 = j % grad_div_size_idx_0 * grad_div_size_idx_1 + j /
              grad_div_size_idx_0;
          } else {
            u0 += grad_div_size_idx_1;
            if (u0 > u1 - 1) {
              u0 = (u0 - u1) + 1;
            }
          }

          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(grad_div_data[u0] == 1) - s;
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
  void wls_var_div1(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                    coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                    &vdops, double result_data[], int result_size[2])
  {
    int grad_size;
    int i;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational divergence operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
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
        for (iRow = 0; iRow < nrows; iRow++) {
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
  void wls_var_div2(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int grad_size;
    int i;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational divergence operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
            b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
            b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_div3(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                    array<double, 2U> &vdops)
  {
    int grad_size;
    int i;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational divergence operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_func1(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int i;
    int iPoint;
    int iRow;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    u1 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial] - 1;
      if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[iPoint + b_wls->V.size
            (1) * j];
        }
      } else {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
            b_wls->V[iPoint + b_wls->V.size(1) * j];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 1.0 - s;
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
        for (iRow = 0; iRow < nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_func2(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 1U> &vdops)
  {
    int i;
    int iPoint;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    u1 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial];
      for (iPoint = 0; iPoint <= npoints; iPoint++) {
        u0 = iMonomial + b_wls->interp0;
        b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[iPoint + b_wls->V.size(1)
          * (j - 1)];
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 1.0 - s;
    }
  }

  static inline
  void wls_var_func3(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 1U> &vdops)
  {
    int i;
    int iPoint;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int u0;
    int u1;

    //  Compute variational WLS-fitting as weighted sum at quadrature points or at a single point.
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(((quad_pnts.size(0) + 3) / 4) << 2, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 0,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    u1 = b_wls->ncols - b_wls->interp0;
    for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
      int j;
      j = b_wls->jpvt[iMonomial] - 1;
      if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[iPoint + b_wls->V.size
            (1) * j];
        }
      } else {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
            b_wls->V[iPoint + b_wls->V.size(1) * j];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 1.0 - s;
    }
  }

  static inline
  void wls_var_grad1(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int grad_size;
    int i;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational gradient operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
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
          for (iRow = 0; iRow < nrows; iRow++) {
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
  void wls_var_grad2(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int grad_size;
    int i;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational gradient operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
            b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
            b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, grad_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                grad_size, b_wls->work);
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_grad3(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    int grad_size;
    int i;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char grad_data[3];

    //  Compute variational gradient operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = grad_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 1,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = grad_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (grad_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, grad_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_grad_div1(WlsObject *b_wls, const ::coder::array<double, 2U>
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
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size;
      int iPoint;
      int j;
      int nDiff;
      int nDims;
      int nOps;
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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, grad_div_size);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;

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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp <= nDiff; iOp++) {
        int offset;

        //  Skip padded zeros in the differential operator
        offset = (grad_div_data[iOp] - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        u1 = b_wls->ncols - b_wls->interp0;
        for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
          j = b_wls->jpvt[iMonomial] - 1;
          if (ws.size(0) == 0) {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              u0 = iMonomial + b_wls->interp0;
              b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
                b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                b_wls->V.size(1) * j];
            }
          } else {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              u0 = iMonomial + b_wls->interp0;
              b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
                b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint)
                - 1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
      vdops.set_size(nrows_vdops, grad_div_size);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = 0.0 - s;
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
  void wls_var_grad_div2(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int dim;
    int hess_size;
    int i;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int nOps;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = hess_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int b_dim{0}; b_dim <= nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
            b_wls->hs_inv.data[b_dim];
        }
      }
    } else {
      for (int b_dim{0}; b_dim <= nDims; b_dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = hess_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
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
    for (int iOp{0}; iOp < nOps; iOp++) {
      signed char b_i;

      //  Skip padded zeros in the differential operator
      b_i = hess_data[iOp];
      if (b_i > 0) {
        int offset;
        offset = (b_i - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        u1 = b_wls->ncols - b_wls->interp0;
        for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
          j = b_wls->jpvt[iMonomial];
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * (j - 1)];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                nOps, b_wls->work);
    vdops.set_size(nrows_vdops, hess_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
  void wls_var_grad_div3(WlsObject *b_wls, const ::coder::array<double, 2U>
    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::array<double, 2U>
    &vdops)
  {
    int dim;
    signed char grad_div_data[9];
    signed char hess_data[9];

    //  Variational grad-div operators as weighted sum at quadrature points
    dim = b_wls->us.size(1);

    //  compute and reorganize vdops
    if (ws.size(1) <= 1) {
      int grad_div_size;
      int i;
      int iPoint;
      int j;
      int nDiff;
      int nDims;
      int nOps;
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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size - 1;
      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = (quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[b_dim]) *
              b_wls->hs_inv.data[b_dim];
          }
        }
      } else {
        for (int b_dim{0}; b_dim <= nDims; b_dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[b_dim + b_wls->us.size(1) * iPoint] = quad_pnts[b_dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[b_dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
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
      for (int iOp{0}; iOp < nOps; iOp++) {
        signed char b_i;

        //  Skip padded zeros in the differential operator
        b_i = hess_data[iOp];
        if (b_i > 0) {
          int offset;
          offset = (b_i - 1) * stride;

          //  Sum up monomials weighted by weights for each component
          u1 = b_wls->ncols - b_wls->interp0;
          for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
            j = b_wls->jpvt[iMonomial] - 1;
            if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                  b_wls->V.size(1) * j];
              }
            } else {
              for (iPoint = 0; iPoint <= npoints; iPoint++) {
                u0 = iMonomial + b_wls->interp0;
                b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0
                  + b_wls->vdops.size(1) * iOp] + ws[ws.size(1) * iPoint] *
                  b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
              }
            }
          }
        }
      }

      //  Multiply by generalized inverse of Vandermonde matrix
      rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, nOps);
      rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank,
                  b_wls->vdops, nOps, b_wls->work);
      vdops.set_size(nrows_vdops, grad_div_size);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = static_cast<double>(hess_data[j] == 1) - s;
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
      int i;
      int iPoint;
      int iWeight;
      int j;
      int lenWs;
      int nDiff;
      int nDims;
      int npoints;
      int nrows;
      int nrows_vdops;
      int stride;
      int u0;
      int u1;

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
      npoints = quad_pnts.size(0) - 1;
      nDims = quad_pnts.size(1) - 1;
      nDiff = grad_div_size - 1;
      if (ws.size(0) == 0) {
        lenWs = 1;
      } else {
        lenWs = ws.size(1);
      }

      stride = ((quad_pnts.size(0) + 3) / 4) << 2;

      //  Scale the coordinates; use wls.us as buffer
      b_wls->us.set_size(stride, quad_pnts.size(1));
      if (b_wls->interp0 != 0) {
        //  Coordinate system is centered at first node in interp0 mode
        for (dim = 0; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
              quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
              b_wls->hs_inv.data[dim];
          }
        }
      } else {
        for (dim = 0; dim <= nDims; dim++) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
              quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
          }
        }
      }

      //  Compute the confluent Vandermonde matrix and right-hand side
      gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
                 b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
      u0 = b_wls->ncols;
      u1 = b_wls->nrows;
      if (u0 >= u1) {
        nrows_vdops = u0;
      } else {
        nrows_vdops = u1;
      }

      //  Force wls.vdops to be varsize and each operator to be stored contiguously
      m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
      u0 = nrows_vdops - b_wls->interp0;
      u1 = grad_div_size + b_wls->izero_;
      b_wls->vdops.set_size(u1, u0);
      u0 *= u1;
      for (u1 = 0; u1 < u0; u1++) {
        b_wls->vdops[u1] = 0.0;
      }

      //  Omit zeros in the diff operators
      iWeight = 1;

      //  Loop through the operators
      for (int iOp{0}; iOp <= nDiff; iOp++) {
        int offset;

        //  Skip padded zeros in the differential operator
        offset = (grad_div_data[iOp] - 1) * stride;

        //  Sum up monomials weighted by weights for each component
        u1 = b_wls->ncols - b_wls->interp0;
        for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
          j = b_wls->jpvt[iMonomial] - 1;
          if (ws.size(0) == 0) {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              u0 = iMonomial + b_wls->interp0;
              b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
                b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
                b_wls->V.size(1) * j];
            }
          } else {
            for (iPoint = 0; iPoint <= npoints; iPoint++) {
              u0 = iMonomial + b_wls->interp0;
              b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
                b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint)
                - 1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
      vdops.set_size(nrows_vdops, grad_div_size);

      //  Transpose the operator for row-major
      u1 = nrows_vdops - b_wls->interp0;
      for (i = 0; i < u1; i++) {
        for (j = 0; j <= nDiff; j++) {
          vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
            b_wls->vdops.size(1) * j];
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

      if (b_wls->interp0 != 0) {
        //  In interp0 mode, we set the first entry based on partition of unity
        for (j = 0; j <= nDiff; j++) {
          double s;
          s = 0.0;
          for (i = 2; i <= nrows; i++) {
            s += vdops[j + vdops.size(1) * (i - 1)];
          }

          vdops[j] = 0.0 - s;
        }
      }
    }

    //  compute output value
  }

  static inline
  void wls_var_hess1(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                     coder::array<double, 2U> &fs, ::coder::array<double, 2U>
                     &vdops, ::coder::array<double, 2U> &result)
  {
    int hess_size;
    int i;
    int iOp;
    int iPoint;
    int iRow;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];

    //  Compute variational Hessian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = hess_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = hess_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (iOp = 0; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, hess_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
      }
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int k{0}; k <= nDiff; k++) {
        for (iRow = 0; iRow < nrows; iRow++) {
          vdops[k + vdops.size(1) * iRow] = vdops[k + vdops.size(1) * iRow] *
            b_wls->rweights[iRow];
        }
      }
    }

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
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
          for (iRow = 0; iRow < nrows; iRow++) {
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
  void wls_var_hess2(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, ::coder::array<double, 2U> &vdops)
  {
    int hess_size;
    int i;
    int iPoint;
    int j;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];

    //  Compute variational Hessian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = hess_size - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = hess_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
            b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
            b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, hess_size);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                hess_size, b_wls->work);
    vdops.set_size(nrows_vdops, hess_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_hess3(WlsObject *b_wls, const ::coder::array<double, 2U>
                     &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                     array<double, 2U> &vdops)
  {
    int hess_size;
    int i;
    int iPoint;
    int iWeight;
    int j;
    int lenWs;
    int nDiff;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char hess_data[6];

    //  Compute variational Hessian operators as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    nDiff = hess_size - 1;
    if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
      lenWs = 1;
    } else {
      lenWs = ws.size(1);
    }

    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, 2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = hess_size + b_wls->izero_;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    iWeight = 1;

    //  Loop through the operators
    for (int iOp{0}; iOp <= nDiff; iOp++) {
      int offset;

      //  Skip padded zeros in the differential operator
      offset = (hess_data[iOp] - 1) * stride;

      //  Sum up monomials weighted by weights for each component
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0 + b_wls->vdops.size(1) * iOp] = b_wls->vdops[u0 +
              b_wls->vdops.size(1) * iOp] + ws[(iWeight + ws.size(1) * iPoint) -
              1] * b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
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
    vdops.set_size(nrows_vdops, hess_size);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      for (j = 0; j <= nDiff; j++) {
        vdops[j + vdops.size(1) * (i + b_wls->interp0)] = b_wls->vdops[i +
          b_wls->vdops.size(1) * j];
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

    if (b_wls->interp0 != 0) {
      //  In interp0 mode, we set the first entry based on partition of unity
      for (j = 0; j <= nDiff; j++) {
        double s;
        s = 0.0;
        for (i = 2; i <= nrows; i++) {
          s += vdops[j + vdops.size(1) * (i - 1)];
        }

        vdops[j] = 0.0 - s;
      }
    }
  }

  static inline
  void wls_var_lap1(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, const ::
                    coder::array<double, 2U> &fs, ::coder::array<double, 1U>
                    &vdops, ::coder::array<double, 2U> &result)
  {
    int i;
    int iPoint;
    int iRow;
    int lap_size_idx_1;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
              b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (iRow = 0; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
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
        for (iRow = 0; iRow < nrows; iRow++) {
          result[iFunc] = result[iFunc] + fs[iFunc + fs.size(1) * iRow] *
            vdops[iRow];
        }
      }
    }
  }

  static inline
  void wls_var_lap2(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, ::coder::array<double, 1U> &vdops, int
                    result_size[2])
  {
    int i;
    int iPoint;
    int lap_size_idx_1;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial];
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          u0 = iMonomial + b_wls->interp0;
          b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
            b_wls->V.size(1) * (j - 1)];
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
    }

    result_size[1] = 0;
    result_size[0] = 0;
  }

  static inline
  void wls_var_lap3(WlsObject *b_wls, const ::coder::array<double, 2U>
                    &quad_pnts, const ::coder::array<double, 2U> &ws, ::coder::
                    array<double, 1U> &vdops, int result_size[2])
  {
    int i;
    int iPoint;
    int lap_size_idx_1;
    int nDims;
    int npoints;
    int nrows;
    int nrows_vdops;
    int stride;
    int u0;
    int u1;
    signed char lap_data[3];

    //  Compute variational (or vector) Laplacian as weighted sum at quadrature points
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
    npoints = quad_pnts.size(0) - 1;
    nDims = quad_pnts.size(1) - 1;
    stride = ((quad_pnts.size(0) + 3) / 4) << 2;

    //  Scale the coordinates; use wls.us as buffer
    b_wls->us.set_size(stride, quad_pnts.size(1));
    if (b_wls->interp0 != 0) {
      //  Coordinate system is centered at first node in interp0 mode
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = (quad_pnts[dim +
            quad_pnts.size(1) * iPoint] - b_wls->origin.data[dim]) *
            b_wls->hs_inv.data[dim];
        }
      }
    } else {
      for (int dim{0}; dim <= nDims; dim++) {
        for (iPoint = 0; iPoint <= npoints; iPoint++) {
          b_wls->us[dim + b_wls->us.size(1) * iPoint] = quad_pnts[dim +
            quad_pnts.size(1) * iPoint] * b_wls->hs_inv.data[dim];
        }
      }
    }

    //  Compute the confluent Vandermonde matrix and right-hand side
    gen_vander(b_wls->us, quad_pnts.size(0), b_wls->degree, -2,
               b_wls->hs_inv.data, b_wls->hs_inv.size, b_wls->V);
    u0 = b_wls->ncols;
    u1 = b_wls->nrows;
    if (u0 >= u1) {
      nrows_vdops = u0;
    } else {
      nrows_vdops = u1;
    }

    //  Force wls.vdops to be varsize and each operator to be stored contiguously
    m2c_assert_msg(b_wls->izero_ == 0, "wls.izero_ must be equal to 0.");
    u0 = nrows_vdops - b_wls->interp0;
    u1 = b_wls->izero_ + 1;
    b_wls->vdops.set_size(u1, u0);
    u0 *= u1;
    for (u1 = 0; u1 < u0; u1++) {
      b_wls->vdops[u1] = 0.0;
    }

    //  Omit zeros in the diff operators
    for (int jDiff{0}; jDiff < lap_size_idx_1; jDiff++) {
      int offset;

      //  Loop through the operators
      offset = (lap_data[jDiff] - 1) * stride;

      //  Skip padded zeros in the differential operator
      u1 = b_wls->ncols - b_wls->interp0;
      for (int iMonomial{0}; iMonomial < u1; iMonomial++) {
        int j;
        j = b_wls->jpvt[iMonomial] - 1;
        if ((ws.size(0) == 0) || (ws.size(1) == 0)) {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + b_wls->V[(offset + iPoint) +
              b_wls->V.size(1) * j];
          }
        } else {
          for (iPoint = 0; iPoint <= npoints; iPoint++) {
            u0 = iMonomial + b_wls->interp0;
            b_wls->vdops[u0] = b_wls->vdops[u0] + ws[ws.size(1) * iPoint] *
              b_wls->V[(offset + iPoint) + b_wls->V.size(1) * j];
          }
        }
      }
    }

    //  Multiply by generalized inverse of Vandermonde matrix
    rrqr_rtsolve(b_wls->QR, b_wls->ncols, b_wls->rank, b_wls->vdops, 1);
    rrqr_qmulti(b_wls->QR, b_wls->nrows, b_wls->ncols, b_wls->rank, b_wls->vdops,
                1, b_wls->work);
    vdops.set_size(nrows_vdops);

    //  Transpose the operator for row-major
    u1 = nrows_vdops - b_wls->interp0;
    for (i = 0; i < u1; i++) {
      vdops[i + b_wls->interp0] = b_wls->vdops[i];
    }

    nrows = b_wls->nrows;
    if (b_wls->rweights.size(0) != 0) {
      for (int iRow{0}; iRow < nrows; iRow++) {
        vdops[iRow] = vdops[iRow] * b_wls->rweights[iRow];
      }
    }

    if (b_wls->interp0 != 0) {
      double s;

      //  In interp0 mode, we set the first entry based on partition of unity
      s = 0.0;
      for (i = 2; i <= nrows; i++) {
        s += vdops[i - 1];
      }

      vdops[0] = 0.0 - s;
    }

    result_size[1] = 0;
    result_size[0] = 0;
  }
}

// End of code generation (wls_internal.cpp)
