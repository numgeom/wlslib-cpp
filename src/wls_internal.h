// Copyright 2022 The NumGeom Group, Stony Brook University
//
// wls_internal.h
//
// Code generation for function 'wls_internal'
//

#ifndef WLS_INTERNAL_H
#define WLS_INTERNAL_H

// Include files
#include "coder_array.h"
#include "m2c_lib.h"
#include "rtwtypes.h"
#include "wls_internal_types.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace wls {
struct WlsObject;

struct WlsWeight;

} // namespace wls

// Function Declarations
namespace wls {
static inline void
wls_func(WlsObject *b_wls, const ::coder::array<real_T, 2U> &pnts,
         const ::coder::array<real_T, 2U> &fs, coder::SizeType npoints,
         ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &result);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &pnts,
                            const ::coder::array<real_T, 2U> &fs,
                            ::coder::array<real_T, 2U> &vdops,
                            ::coder::array<real_T, 2U> &result);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &pnts,
                            ::coder::array<real_T, 2U> &vdops);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree,
                            coder::SizeType order, coder::SizeType interp0,
                            boolean_T use_dag, coder::SizeType npoints);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            coder::SizeType degree);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            coder::SizeType degree, coder::SizeType order);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            coder::SizeType degree, coder::SizeType order,
                            coder::SizeType interp0);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            coder::SizeType degree, coder::SizeType order,
                            coder::SizeType interp0, boolean_T use_dag);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree,
                            coder::SizeType order);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree,
                            coder::SizeType order, coder::SizeType interp0);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree,
                            coder::SizeType order, coder::SizeType interp0,
                            boolean_T use_dag);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            coder::SizeType degree, coder::SizeType order,
                            coder::SizeType interp0, boolean_T use_dag,
                            coder::SizeType npoints);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &quad_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 const ::coder::array<real_T, 2U> &varargin_2,
                                 coder::SizeType varargin_3,
                                 ::coder::array<real_T, 2U> &varargout_1,
                                 ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &quad_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 const ::coder::array<real_T, 2U> &varargin_2,
                                 ::coder::array<real_T, 2U> &varargout_1,
                                 ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &quad_pnts,
                                 ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &quad_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
             const ::coder::array<real_T, 2U> &ws,
             const ::coder::array<real_T, 2U> &fs, coder::SizeType varargin_1,
             ::coder::array<real_T, 2U> &vdops, const real_T result_data[],
             coder::SizeType result_size[2]);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                const ::coder::array<real_T, 2U> &fs,
                                ::coder::array<real_T, 2U> &vdops,
                                real_T result_data[],
                                coder::SizeType result_size[2]);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  const ::coder::array<real_T, 2U> &fs,
                  coder::SizeType varargin_1, ::coder::array<real_T, 2U> &vdops,
                  const real_T result_data[], coder::SizeType result_size[2]);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  const ::coder::array<real_T, 2U> &fs,
                  ::coder::array<real_T, 2U> &vdops, real_T result_data[],
                  coder::SizeType result_size[2]);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
                  ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               coder::SizeType varargin_3,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_grad_div(WlsObject *b_wls, const ::coder::array<real_T, 2U> &quad_pnts,
                 const ::coder::array<real_T, 2U> &ws,
                 const ::coder::array<real_T, 2U> &fs,
                 coder::SizeType varargin_1, ::coder::array<real_T, 2U> &vdops,
                 const real_T result_data[], coder::SizeType result_size[2]);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &quad_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    const ::coder::array<real_T, 2U> &fs,
                                    ::coder::array<real_T, 2U> &vdops,
                                    real_T result_data[],
                                    coder::SizeType result_size[2]);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &quad_pnts,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &quad_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &quad_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               coder::SizeType varargin_3,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &quad_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &varargout_1);

} // namespace wls

#include "wls_internal.cpp"
#endif
// End of code generation (wls_internal.h)
