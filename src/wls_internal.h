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
static inline void wls_eval_ops(const WlsObject *b_wls,
                                const ::coder::array<int8_T, 2U> &vdops,
                                const ::coder::array<real_T, 2U> &fs,
                                boolean_T isdiv, coder::SizeType nOps,
                                ::coder::array<real_T, 2U> &result);

static inline void
wls_func(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
         const ::coder::array<real_T, 2U> &fs, coder::SizeType nevpnts,
         ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &result);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &fs,
                            ::coder::array<real_T, 2U> &vdops,
                            ::coder::array<real_T, 2U> &result);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            ::coder::array<real_T, 2U> &vdops);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, coder::SizeType degree,
                            coder::SizeType order, coder::SizeType interp0,
                            boolean_T use_dag, coder::SizeType nstpnts);

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
                            coder::SizeType nstpnts);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight);

static inline void wls_solve_sys(WlsObject *b_wls,
                                 ::coder::array<real_T, 2U> &vdops);

static inline void wls_solve_sys(WlsObject *b_wls,
                                 const ::coder::array<int8_T, 2U> &diff_idx,
                                 ::coder::array<real_T, 2U> &vdops);

static inline void wls_solve_sys(WlsObject *b_wls,
                                 const ::coder::array<int8_T, 2U> &diff_idx,
                                 const ::coder::array<real_T, 2U> &ws,
                                 ::coder::array<real_T, 2U> &vdops);

static inline void wls_tab_mbasis(WlsObject *b_wls,
                                  const ::coder::array<real_T, 2U> &evpnts);

static inline void wls_tab_mbasis(WlsObject *b_wls,
                                  const ::coder::array<real_T, 2U> &evpnts,
                                  coder::SizeType order);

static inline void wls_tab_mbasis(WlsObject *b_wls,
                                  const ::coder::array<real_T, 2U> &evpnts,
                                  coder::SizeType order,
                                  coder::SizeType nevpnts);

static inline void wls_update_rhs(WlsObject *b_wls,
                                  const ::coder::array<int8_T, 2U> &diff_idx);

static inline void wls_update_rhs2(WlsObject *b_wls,
                                   const ::coder::array<int8_T, 2U> &diff_idx,
                                   const ::coder::array<real_T, 2U> &ws);

static inline void wls_update_rhs3(WlsObject *b_wls,
                                   const ::coder::array<int8_T, 2U> &diff_idx,
                                   const ::coder::array<real_T, 2U> &ws);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &eval_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 const ::coder::array<real_T, 2U> &varargin_2,
                                 coder::SizeType varargin_3,
                                 ::coder::array<real_T, 2U> &varargout_1,
                                 ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &eval_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 const ::coder::array<real_T, 2U> &varargin_2,
                                 ::coder::array<real_T, 2U> &varargout_1,
                                 ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &eval_pnts,
                                 ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &eval_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws_lap,
                                    const ::coder::array<real_T, 2U> &ws_grad,
                                    const ::coder::array<real_T, 2U> &fs,
                                    coder::SizeType nevpnts,
                                    ::coder::array<real_T, 2U> &vdops,
                                    ::coder::array<real_T, 2U> &result);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws_lap,
                                    const ::coder::array<real_T, 2U> &ws_grad,
                                    const ::coder::array<real_T, 2U> &fs,
                                    ::coder::array<real_T, 2U> &vdops,
                                    ::coder::array<real_T, 2U> &result);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws_lap,
                                    const ::coder::array<real_T, 2U> &ws_grad,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws_lap,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
             const ::coder::array<real_T, 2U> &ws,
             const ::coder::array<real_T, 2U> &fs, coder::SizeType varargin_1,
             ::coder::array<real_T, 2U> &vdops, const real_T result_data[],
             coder::SizeType result_size[2]);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                const ::coder::array<real_T, 2U> &fs,
                                ::coder::array<real_T, 2U> &vdops,
                                real_T result_data[],
                                coder::SizeType result_size[2]);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  const ::coder::array<real_T, 2U> &fs,
                  coder::SizeType varargin_1, ::coder::array<real_T, 2U> &vdops,
                  const real_T result_data[], coder::SizeType result_size[2]);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  const ::coder::array<real_T, 2U> &fs,
                  ::coder::array<real_T, 2U> &vdops, real_T result_data[],
                  coder::SizeType result_size[2]);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  const ::coder::array<real_T, 2U> &ws,
                  ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_div(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
            const ::coder::array<real_T, 2U> &varargin_1,
            const ::coder::array<real_T, 2U> &varargin_2,
            coder::SizeType varargin_3, ::coder::array<real_T, 2U> &varargout_1,
            real_T varargout_2_data[], coder::SizeType varargout_2_size[1]);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               ::coder::array<real_T, 2U> &varargout_1,
                               real_T varargout_2_data[],
                               coder::SizeType varargout_2_size[1]);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_grad(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_grad_div(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                 const ::coder::array<real_T, 2U> &ws,
                 const ::coder::array<real_T, 2U> &fs,
                 coder::SizeType varargin_1, ::coder::array<real_T, 2U> &vdops,
                 const real_T result_data[], coder::SizeType result_size[2]);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    const ::coder::array<real_T, 2U> &fs,
                                    ::coder::array<real_T, 2U> &vdops,
                                    real_T result_data[],
                                    coder::SizeType result_size[2]);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                coder::SizeType varargin_3,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_hess(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_kernel(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    coder::SizeType order, const ::coder::array<int8_T, 2U> &diff_idx,
    const ::coder::array<real_T, 2U> &ws, const ::coder::array<real_T, 2U> &fs,
    coder::SizeType nevpnts, ::coder::array<real_T, 2U> &vdops,
    ::coder::array<real_T, 2U> &result);

static inline void wls_var_kernel(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    coder::SizeType order, const ::coder::array<int8_T, 2U> &diff_idx,
    const ::coder::array<real_T, 2U> &ws, const ::coder::array<real_T, 2U> &fs,
    ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &result);

static inline void wls_var_kernel(WlsObject *b_wls,
                                  const ::coder::array<real_T, 2U> &eval_pnts,
                                  coder::SizeType order,
                                  const ::coder::array<int8_T, 2U> &diff_idx,
                                  ::coder::array<real_T, 2U> &vdops,
                                  ::coder::array<real_T, 2U> &result);

static inline void wls_var_kernel(WlsObject *b_wls,
                                  const ::coder::array<real_T, 2U> &eval_pnts,
                                  coder::SizeType order,
                                  const ::coder::array<int8_T, 2U> &diff_idx,
                                  const ::coder::array<real_T, 2U> &ws,
                                  ::coder::array<real_T, 2U> &vdops,
                                  ::coder::array<real_T, 2U> &result);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               coder::SizeType varargin_3,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               const ::coder::array<real_T, 2U> &varargin_2,
                               ::coder::array<real_T, 2U> &varargout_1,
                               ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_sym_grad(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    const ::coder::array<real_T, 2U> &fs,
                                    ::coder::array<real_T, 2U> &vdops,
                                    real_T result_data[],
                                    coder::SizeType result_size[2]);

static inline void wls_var_sym_grad(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &vdops,
                                    const real_T result_data[],
                                    coder::SizeType result_size[2]);

static inline void wls_var_sym_grad(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    const ::coder::array<real_T, 2U> &ws,
                                    ::coder::array<real_T, 2U> &vdops,
                                    const real_T result_data[],
                                    coder::SizeType result_size[2]);

} // namespace wls

#include "wls_internal.cpp"
#endif
// End of code generation (wls_internal.h)
