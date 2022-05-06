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
#ifndef CODER_ARRAY_SIZE_TYPE_DEFINED
#define CODER_ARRAY_SIZE_TYPE_DEFINED
namespace coder {
#ifdef M2C_SIZETYPE64
typedef int64_T SizeType;
#else
typedef int32_T SizeType;
#endif // M2C_SIZETYPE64
} // namespace coder
#endif // CODER_ARRAY_SIZE_TYPE_DEFINED

#include "coder_array.h"
#include "coder_bounded_array.h"

// Type Definitions
namespace wls {
struct WlsObject {
  int32_T npoints;
  int32_T degree;
  int32_T order;
  int32_T interp0;
  boolean_T use_dag;
  int32_T stride;
  ::coder::array<real_T, 2U> us;
  ::coder::bounded_array<real_T, 3U, 2U> origin;
  ::coder::array<real_T, 1U> rweights;
  ::coder::bounded_array<real_T, 3U, 2U> hs_inv;
  ::coder::array<uint8_T, 2U> dag;
  ::coder::array<real_T, 2U> V;
  ::coder::array<real_T, 2U> QR;
  ::coder::array<real_T, 2U> vdops;
  int32_T nrows;
  int32_T ncols;
  int32_T rank;
  boolean_T fullrank;
  ::coder::array<int32_T, 1U> jpvt;
  ::coder::array<real_T, 1U> work;
  boolean_T rowmajor;
  ::coder::array<real_T, 2U> QRt;
};

struct WlsWeight {
  ::coder::array<char_T, 2U> name;
  ::coder::array<real_T, 1U> params_shared;
  ::coder::array<real_T, 2U> params_pointwise;
  ::coder::array<boolean_T, 1U> omit_rows;
};

} // namespace wls

#endif
// End of code generation (wls_internal_types.h)
