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
static inline void wls_curl(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::SizeType varargin_2,
                            ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_curl(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_curl(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_div(WlsObject *b_wls,
                           const ::coder::array<real_T, 2U> &eval_pnts,
                           const ::coder::array<real_T, 2U> &varargin_1,
                           ::coder::SizeType varargin_2,
                           ::coder::array<real_T, 2U> &varargout_1,
                           ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_div(WlsObject *b_wls,
                           const ::coder::array<real_T, 2U> &eval_pnts,
                           const ::coder::array<real_T, 2U> &varargin_1,
                           ::coder::array<real_T, 2U> &varargout_1,
                           ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_div(WlsObject *b_wls,
                           const ::coder::array<real_T, 2U> &eval_pnts,
                           ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::SizeType varargin_2,
                            ::coder::array<real_T, 2U> &varargout_1,
                            ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::array<real_T, 2U> &varargout_1,
                            ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_func(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_grad(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::SizeType varargin_2,
                            ::coder::array<real_T, 2U> &varargout_1,
                            ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_grad(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            const ::coder::array<real_T, 2U> &varargin_1,
                            ::coder::array<real_T, 2U> &varargout_1,
                            ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_grad(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &eval_pnts,
                            ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree,
                            ::coder::SizeType order, ::coder::SizeType interp0,
                            ::coder::SizeType unimono);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree, ::coder::SizeType order);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree, ::coder::SizeType order,
                            ::coder::SizeType interp0);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree, ::coder::SizeType order,
                            ::coder::SizeType interp0, boolean_T unimono);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree, ::coder::SizeType order,
                            ::coder::SizeType interp0, boolean_T unimono,
                            ::coder::SizeType nstpnts);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree,
                            ::coder::SizeType order);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree,
                            ::coder::SizeType order, ::coder::SizeType interp0);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree,
                            ::coder::SizeType order, ::coder::SizeType interp0,
                            boolean_T unimono);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const WlsWeight *weight, ::coder::SizeType degree,
                            ::coder::SizeType order, ::coder::SizeType interp0,
                            boolean_T unimono, ::coder::SizeType nstpnts);

static inline void wls_init(WlsObject *b_wls,
                            const ::coder::array<real_T, 2U> &us,
                            const ::coder::array<char_T, 2U> &weight,
                            ::coder::SizeType degree, ::coder::SizeType order,
                            ::coder::SizeType interp0,
                            ::coder::SizeType unimono);

static inline void wls_var_bilap(WlsObject *b_wls,
                                 const ::coder::array<real_T, 2U> &eval_pnts,
                                 const ::coder::array<real_T, 2U> &varargin_1,
                                 const ::coder::array<real_T, 2U> &varargin_2,
                                 ::coder::SizeType varargin_3,
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

static inline void wls_var_cdr(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 2U> &ws_rho,
                               const ::coder::array<real_T, 2U> &fs,
                               ::coder::SizeType nevpnts,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &result);

static inline void wls_var_cdr(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 2U> &ws_rho,
                               const ::coder::array<real_T, 2U> &fs,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &result);

static inline void wls_var_cdr(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 2U> &ws_rho,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_cdr4(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws_lap,
                                const ::coder::array<real_T, 2U> &ws_grad,
                                ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_cdr5(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws_lap,
                                ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_convdiff(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    ::coder::array<real_T, 2U> &ws_lap, ::coder::array<real_T, 2U> &ws_grad,
    const ::coder::array<real_T, 2U> &varargin_1, ::coder::SizeType varargin_2,
    ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_convdiff(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    ::coder::array<real_T, 2U> &ws_lap, ::coder::array<real_T, 2U> &ws_grad,
    const ::coder::array<real_T, 2U> &varargin_1,
    ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &ws_lap,
                                    ::coder::array<real_T, 2U> &ws_grad,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_convdiff(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &ws_lap,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::SizeType varargin_2,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                ::coder::array<real_T, 2U> &varargout_1,
                                ::coder::array<real_T, 2U> &varargout_2);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_curl(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &ws,
                                ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_curl_curl(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    ::coder::array<real_T, 1U> &ws,
    const ::coder::array<real_T, 2U> &varargin_1, ::coder::SizeType varargin_2,
    ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  ::coder::array<real_T, 1U> &ws,
                  const ::coder::array<real_T, 2U> &varargin_1,
                  ::coder::array<real_T, 2U> &vdops,
                  ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  ::coder::array<real_T, 2U> &vdops);

static inline void
wls_var_curl_curl(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                  ::coder::array<real_T, 1U> &ws,
                  ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::SizeType varargin_2,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_div(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 2U> &ws,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_func(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 2U> &varargin_1,
                                const ::coder::array<real_T, 2U> &varargin_2,
                                ::coder::SizeType varargin_3,
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
                                ::coder::SizeType varargin_3,
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

static inline void wls_var_grad_div(
    WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
    ::coder::array<real_T, 1U> &ws,
    const ::coder::array<real_T, 2U> &varargin_1, ::coder::SizeType varargin_2,
    ::coder::array<real_T, 2U> &vdops, ::coder::array<real_T, 2U> &varargout_1);

static inline void
wls_var_grad_div(WlsObject *b_wls, const ::coder::array<real_T, 2U> &eval_pnts,
                 ::coder::array<real_T, 1U> &ws,
                 const ::coder::array<real_T, 2U> &varargin_1,
                 ::coder::array<real_T, 2U> &vdops,
                 ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_grad_div(WlsObject *b_wls,
                                    const ::coder::array<real_T, 2U> &eval_pnts,
                                    ::coder::array<real_T, 1U> &ws,
                                    ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &ws,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::SizeType varargin_2,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &ws,
                               const ::coder::array<real_T, 2U> &varargin_1,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &varargout_1);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_lap(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               ::coder::array<real_T, 2U> &ws,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_uno(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 1U> &ws_graddiv,
                               const ::coder::array<real_T, 1U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 1U> &ws_rho,
                               const ::coder::array<real_T, 2U> &fs,
                               ::coder::SizeType nevpnts,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &result);

static inline void wls_var_uno(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 1U> &ws_graddiv,
                               const ::coder::array<real_T, 1U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 1U> &ws_rho,
                               const ::coder::array<real_T, 2U> &fs,
                               ::coder::array<real_T, 2U> &vdops,
                               ::coder::array<real_T, 2U> &result);

static inline void wls_var_uno(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 1U> &ws_graddiv,
                               const ::coder::array<real_T, 1U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               const ::coder::array<real_T, 1U> &ws_rho,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_uno(WlsObject *b_wls,
                               const ::coder::array<real_T, 2U> &eval_pnts,
                               const ::coder::array<real_T, 1U> &ws_graddiv,
                               const ::coder::array<real_T, 1U> &ws_lap,
                               const ::coder::array<real_T, 2U> &ws_grad,
                               ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_uno5(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 1U> &ws_graddiv,
                                const ::coder::array<real_T, 1U> &ws_lap,
                                ::coder::array<real_T, 2U> &vdops);

static inline void wls_var_uno6(WlsObject *b_wls,
                                const ::coder::array<real_T, 2U> &eval_pnts,
                                const ::coder::array<real_T, 1U> &ws_graddiv,
                                ::coder::array<real_T, 2U> &vdops);

} // namespace wls

#include "wls_internal.cpp"
#endif
// End of code generation (wls_internal.h)
