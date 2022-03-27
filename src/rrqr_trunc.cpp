// Copyright 2022 The NumGeom Group, Stony Brook University
// Main developers:
//     wlslib: Xiangmin Jiao, Qiao Chen, Jacob Jones
//     momp2cpp: Xiangmin Jiao, Qiao Chen
//
// rrqr_trunc.cpp
//
// Code generation for function 'rrqr_trunc'
//

// Include files
#include "rrqr_trunc.h"
#include "coder_array.h"
#include "m2c_lib.h"
#include <cstdio>
#include <stdexcept>

// Function Definitions
namespace wls {
//  rrqr_trunc  Truncate the output from QRCP based on DAG
void rrqr_trunc(const ::coder::array<unsigned char, 2U> &dag, int *n1, int rank,
                ::coder::array<int, 1U> &p, ::coder::array<int, 2U> &work,
                boolean_T *permuted)
{
  int b_i;
  int child;
  int n;
  int nChanged;
  n = dag.size(1) - 2;
  //  Last column of dag is for signature
  work.set_size(4, dag.size(1) - 1);
  //  compute inverse permutation and tag whether each column is truncated
  for (int i{0}; i <= n; i++) {
    work[4 * (p[i] - 1)] = i + 1;
    //  first column for inverse permute
    work[4 * (p[i] - 1) + 1] = i + 1 > rank;
    //  second column for tagging truncation
    work[4 * i + 3] = 0;
    //  third column for stack, and last column for marks
  }
  nChanged = -1;
  //  If a truncated monomial has an untruncated child, do not truncate it.
  b_i = rank + 1;
  for (int i{b_i}; i <= *n1; i++) {
    int c_tmp;
    c_tmp = p[i - 1] - 1;
    if (work[4 * c_tmp + 1] != 0) {
      int i1;
      i1 = dag.size(0);
      for (int j{0}; j < i1; j++) {
        child = c_tmp + dag[j + dag.size(0) * c_tmp];
        if ((dag[j + dag.size(0) * c_tmp] != 0) && (work[4 * child + 1] != 1)) {
          if (work[4 * child] <= rank) {
            work[4 * child + 1] = -1;
            if (work[4 * child + 3] == 0) {
              nChanged++;
              work[4 * nChanged + 2] = child + 1;
              work[4 * child + 3] = 1;
            }
          }
          work[4 * c_tmp + 1] = 0;
        }
      }
    }
  }
  *permuted = nChanged + 1 != 0;
  //  If a monomial is truncated but some of its children are not, then make
  while (nChanged + 1 != 0) {
    int c;
    boolean_T allChildrenTruncated;
    c = work[4 * nChanged + 2] - 1;
    nChanged--;
    //  Make sure all the untruncated children of a candidate monomial
    allChildrenTruncated = true;
    b_i = dag.size(0);
    for (int j{0}; j < b_i; j++) {
      if (dag[j + dag.size(0) * c] != 0) {
        child = c + dag[j + dag.size(0) * c];
        if (work[4 * child] <= rank) {
          allChildrenTruncated = false;
          if (work[4 * child + 1] != -1) {
            work[4 * child + 1] = -1;
            if (work[4 * child + 3] == 0) {
              nChanged++;
              work[4 * nChanged + 2] = child + 1;
              work[4 * child + 3] = 1;
            }
          }
        }
      }
    }
    if (!allChildrenTruncated) {
      work[4 * c + 1] = 0;
    }
  }
  if (*permuted) {
    int n2;
    //  permute the truncated columns to the end
    *n1 = 0;
    n2 = dag.size(1) - 1;
    for (int i{0}; i <= n; i++) {
      //  i is original ids
      if (work[4 * i + 1] == 0) {
        (*n1)++;
        p[*n1 - 1] = i + 1;
      } else {
        p[n2 - 1] = i + 1;
        n2--;
      }
    }
    m2cAssert(*n1 == n2, "Post-condition check failed");
  }
}

} // namespace wls

// End of code generation (rrqr_trunc.cpp)
