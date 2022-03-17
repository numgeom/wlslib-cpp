// Copyright 2022 The NumGeom Group, Stony Brook University
//
// wls_internal_types.h
//
// Code generation for function 'wls_init'
//

#ifndef WLS_INTERNAL_TYPES_H
#define WLS_INTERNAL_TYPES_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "coder_bounded_array.h"

// Type Definitions
namespace wls {
struct WlsObject {
  int npoints;
  int degree;
  int order;
  int interp0;
  boolean_T use_dag;
  int stride;
  ::coder::array<double, 2U> us;
  ::coder::bounded_array<double, 3U, 2U> origin;
  ::coder::array<double, 1U> rweights;
  ::coder::bounded_array<double, 3U, 2U> hs_inv;
  ::coder::array<unsigned char, 2U> dag;
  ::coder::array<double, 2U> V;
  ::coder::array<double, 2U> QR;
  ::coder::array<double, 2U> vdops;
  int nrows;
  int ncols;
  int rank;
  ::coder::array<int, 1U> jpvt;
  ::coder::array<double, 1U> work;
  boolean_T rowmajor;
  ::coder::array<double, 2U> QRt;
};

struct WlsWeight {
  ::coder::array<char, 2U> name;
  ::coder::array<double, 1U> params_shared;
  ::coder::array<double, 2U> params_pointwise;
  ::coder::array<boolean_T, 1U> omit_rows;
};

} // namespace wls

#endif
// End of code generation (wls_internal_types.h)
