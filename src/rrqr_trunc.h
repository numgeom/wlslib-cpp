// Copyright 2022 The NumGeom Group, Stony Brook University
//
// rrqr_trunc.h
//
// Code generation for function 'rrqr_trunc'
//

#ifndef RRQR_TRUNC_H
#define RRQR_TRUNC_H

// Include files
#include "rtwtypes.h"
#include "m2c_lib.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace wls {
static inline void rrqr_trunc(const ::coder::array<unsigned char, 2U> &dag, int *n1,
                       int rank, ::coder::array<int, 1U> &p,
                       ::coder::array<int, 2U> &work, boolean_T *permuted);

} // namespace wls

#include "rrqr_trunc.cpp"
#endif
// End of code generation (rrqr_trunc.h)
