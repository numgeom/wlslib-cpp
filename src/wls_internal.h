//
// wls_internal.h
//
// Code generation for function 'wls_internal'
//

#ifndef WLS_INTERNAL_H
#define WLS_INTERNAL_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace wls {
struct WlsObject;

struct WlsWeight;

} // namespace wls

// Function Declarations
namespace wls {
static inline void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      const ::coder::array<double, 2U> &fs,
                      ::coder::array<double, 2U> &vdops);

static inline void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      const ::coder::array<double, 2U> &fs, int npoints,
                      ::coder::array<double, 2U> &vdops);

static inline void wls_func(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      ::coder::array<double, 2U> &vdops);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      boolean_T interp0, boolean_T use_dag, int npoints);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      boolean_T interp0);

static inline void wls_init(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      boolean_T interp0, boolean_T use_dag);

static inline void wls_var_bilap(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           const ::coder::array<double, 2U> &ws,
                           const ::coder::array<double, 2U> &fs,
                           ::coder::array<double, 1U> &vdops,
                           ::coder::array<double, 2U> &result);

static inline void wls_var_bilap(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           ::coder::array<double, 1U> &vdops);

static inline void wls_var_bilap(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           const ::coder::array<double, 2U> &ws,
                           ::coder::array<double, 1U> &vdops);

static inline void wls_var_curl(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs,
                          ::coder::array<double, 2U> &vdops,
                          double result_data[], int result_size[2]);

static inline void wls_var_curl(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl_curl(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               const ::coder::array<double, 2U> &ws,
                               const ::coder::array<double, 2U> &fs,
                               ::coder::array<double, 2U> &vdops,
                               double result_data[], int result_size[2]);

static inline void wls_var_curl_curl(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl_curl(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               const ::coder::array<double, 2U> &ws,
                               ::coder::array<double, 2U> &vdops);

static inline void wls_var_div(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &ws,
                         const ::coder::array<double, 2U> &fs,
                         ::coder::array<double, 2U> &vdops,
                         double result_data[], int result_size[2]);

static inline void wls_var_div(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         ::coder::array<double, 2U> &vdops);

static inline void wls_var_div(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &ws,
                         ::coder::array<double, 2U> &vdops);

static inline void wls_var_func(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs,
                          ::coder::array<double, 1U> &vdops,
                          ::coder::array<double, 2U> &result);

static inline void wls_var_func(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 1U> &vdops);

static inline void wls_var_func(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          ::coder::array<double, 1U> &vdops);

static inline void wls_var_grad(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs,
                          ::coder::array<double, 2U> &vdops,
                          ::coder::array<double, 2U> &result);

static inline void wls_var_grad(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_grad(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_grad_div(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              const ::coder::array<double, 2U> &ws,
                              const ::coder::array<double, 2U> &fs,
                              ::coder::array<double, 2U> &vdops,
                              double result_data[], int result_size[2]);

static inline void wls_var_grad_div(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              ::coder::array<double, 2U> &vdops);

static inline void wls_var_grad_div(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              const ::coder::array<double, 2U> &ws,
                              ::coder::array<double, 2U> &vdops);

static inline void wls_var_hess(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs,
                          ::coder::array<double, 2U> &vdops,
                          ::coder::array<double, 2U> &result);

static inline void wls_var_hess(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_hess(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_lap(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &ws,
                         const ::coder::array<double, 2U> &fs,
                         ::coder::array<double, 1U> &vdops,
                         ::coder::array<double, 2U> &result);

static inline void wls_var_lap(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         ::coder::array<double, 1U> &vdops, int result_size[2]);

static inline void wls_var_lap(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &ws,
                         ::coder::array<double, 1U> &vdops, int result_size[2]);

} // namespace wls

#include "wls_internal.cpp"
#endif
// End of code generation (wls_internal.h)
