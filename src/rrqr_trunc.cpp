//
// rrqr_trunc.cpp
//
// Code generation for function 'rrqr_trunc'
//

// Include files
#include "rrqr_trunc.h"
#include "coder_array.h"

// Function Definitions
namespace wls {
static inline
boolean_T rrqr_trunc(const ::coder::array<unsigned char, 2U> &dag, int *n1,
                     int rank, ::coder::array<int, 1U> &p,
                     ::coder::array<int, 2U> &work)
{
  int b_i;
  int child;
  int i;
  int j;
  int n;
  int nChanged;
  boolean_T permuted;
  //  Truncate the output from QRCP based on DAG
  //
  //  [permuted, n1, p] = rrqr_tunc(dag, n1, rank, p)
  //  [permuted, n1, p, work] = rrqr_tunc(dag, n1, rank, p, work)
  //
  //  Parameters
  //  ----------
  //      dag:    Direct Acyclic Graph of the Vandermonde matrix
  //      n1:     The original length of the permutation vector
  //      rank:   Numerical rank of the matrix determined by rrqr
  //      p:      Column pivoting vector (1-based index) (length >= n)
  //      work:   Integer work space; you can start with zeros(1,0) and its
  //              allocation will grow automatically as needed
  //
  //  Output
  //  ------
  //      permuted: Whether the output is permuted
  //      n1:       Number of columns or refactorization
  //      p:        Updated column pivoting vector
  //      work:     Work space returned for reuse.
  //
  //  Note: Do NOT call this function manually. It is meant for code generation
  //        to be called by wls::rrqr_factor in C++ automatically.
  //
  //
  //  See also rrqr_factor, gen_vander_dag
  //  Force to be column major
  n = dag.size(1) - 1;
  //  Resize work space if needed
  if (work.size(1) < dag.size(1)) {
    work.set_size(4, dag.size(1));
  }
  //  compute inverse permutation and tag whether each column is truncated
  for (i = 0; i <= n; i++) {
    work[4 * (p[i] - 1)] = i + 1;
    //  first column for inverse permute
    work[4 * (p[i] - 1) + 1] = (i + 1 > rank);
    //  second column for tagging truncation
    work[4 * i + 3] = 0;
    //  third column for stack, and last column for marks
  }
  nChanged = -1;
  //  If a truncated monomial has an untruncated child, do not truncate it.
  b_i = rank + 1;
  for (i = b_i; i <= *n1; i++) {
    int c_tmp;
    c_tmp = p[i - 1] - 1;
    if (work[4 * c_tmp + 1] != 0) {
      int i1;
      i1 = dag.size(0);
      for (j = 0; j < i1; j++) {
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
  permuted = (nChanged + 1 != 0);
  //  If a monomial is truncated but some of its children are not, then make
  //  sure its untruncated children are all marked as candidate for truncation.
  while (nChanged + 1 != 0) {
    int c;
    boolean_T allChildrenTruncated;
    c = work[4 * nChanged + 2] - 1;
    nChanged--;
    //  Make sure all the untruncated children of a candidate monomial
    //  are considered for truncation.
    allChildrenTruncated = true;
    b_i = dag.size(0);
    for (j = 0; j < b_i; j++) {
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
  if (permuted) {
    int n2;
    //  permute the truncated columns to the end
    *n1 = 0;
    n2 = dag.size(1);
    for (i = 0; i <= n; i++) {
      //  i is original ids
      if (work[4 * i + 1] == 0) {
        (*n1)++;
        p[*n1 - 1] = i + 1;
      } else {
        p[n2 - 1] = i + 1;
        n2--;
      }
    }
  }
  return permuted;
}

} // namespace wls

// End of code generation (rrqr_trunc.cpp)
