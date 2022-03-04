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
static inline void wls_func1(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      const ::coder::array<double, 2U> &fs, int npoints,
                      ::coder::array<double, 2U> &vdops,
                      ::coder::array<double, 2U> &result);

static inline void wls_func2(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      const ::coder::array<double, 2U> &fs,
                      ::coder::array<double, 2U> &vdops,
                      ::coder::array<double, 2U> &result);

static inline void wls_func3(WlsObject *b_wls, const ::coder::array<double, 2U> &pnts,
                      ::coder::array<double, 2U> &vdops);

static inline void wls_init1(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      int interp0, boolean_T use_dag, int npoints);

static inline void wls_init2(WlsObject *b_wls, const ::coder::array<double, 2U> &us);

static inline void wls_init3(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight);

static inline void wls_init4(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree);

static inline void wls_init5(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order);

static inline void wls_init6(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      int interp0);

static inline void wls_init7(WlsObject *b_wls, const ::coder::array<double, 2U> &us,
                      const WlsWeight *weight, int degree, int order,
                      int interp0, boolean_T use_dag);

static inline void wls_var_bilap1(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           const ::coder::array<double, 2U> &varargin_1,
                           const ::coder::array<double, 2U> &varargin_2,
                           int varargin_3,
                           ::coder::array<double, 2U> &varargout_1,
                           ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_bilap2(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           const ::coder::array<double, 2U> &varargin_1,
                           const ::coder::array<double, 2U> &varargin_2,
                           ::coder::array<double, 2U> &varargout_1,
                           ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_bilap3(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_bilap4(WlsObject *b_wls,
                           const ::coder::array<double, 2U> &quad_pnts,
                           const ::coder::array<double, 2U> &varargin_1,
                           ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_curl1(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs, int varargin_1,
                          ::coder::array<double, 2U> &vdops,
                          const double result_data[], int result_size[2]);

static inline void wls_var_curl2(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          const ::coder::array<double, 2U> &fs,
                          ::coder::array<double, 2U> &vdops,
                          double result_data[], int result_size[2]);

static inline void wls_var_curl3(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl4(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &ws,
                          ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl_curl1(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               const ::coder::array<double, 2U> &ws,
                               const ::coder::array<double, 2U> &fs,
                               int varargin_1,
                               ::coder::array<double, 2U> &vdops,
                               const double result_data[], int result_size[2]);

static inline void wls_var_curl_curl2(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               const ::coder::array<double, 2U> &ws,
                               const ::coder::array<double, 2U> &fs,
                               ::coder::array<double, 2U> &vdops,
                               double result_data[], int result_size[2]);

static inline void wls_var_curl_curl3(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               ::coder::array<double, 2U> &vdops);

static inline void wls_var_curl_curl4(WlsObject *b_wls,
                               const ::coder::array<double, 2U> &quad_pnts,
                               const ::coder::array<double, 2U> &ws,
                               ::coder::array<double, 2U> &vdops);

static inline void wls_var_div1(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         const ::coder::array<double, 2U> &varargin_2,
                         int varargin_3,
                         ::coder::array<double, 2U> &varargout_1,
                         ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_div2(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         const ::coder::array<double, 2U> &varargin_2,
                         ::coder::array<double, 2U> &varargout_1,
                         ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_div3(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_div4(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_func1(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          int varargin_3,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_func2(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_func3(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_func4(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_grad1(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          int varargin_3,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_grad2(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_grad3(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_grad4(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_grad_div1(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              const ::coder::array<double, 2U> &ws,
                              const ::coder::array<double, 2U> &fs,
                              int varargin_1, ::coder::array<double, 2U> &vdops,
                              const double result_data[], int result_size[2]);

static inline void wls_var_grad_div2(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              const ::coder::array<double, 2U> &ws,
                              const ::coder::array<double, 2U> &fs,
                              ::coder::array<double, 2U> &vdops,
                              double result_data[], int result_size[2]);

static inline void wls_var_grad_div3(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              ::coder::array<double, 2U> &vdops);

static inline void wls_var_grad_div4(WlsObject *b_wls,
                              const ::coder::array<double, 2U> &quad_pnts,
                              const ::coder::array<double, 2U> &ws,
                              ::coder::array<double, 2U> &vdops);

static inline void wls_var_hess1(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          int varargin_3,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_hess2(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          const ::coder::array<double, 2U> &varargin_2,
                          ::coder::array<double, 2U> &varargout_1,
                          ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_hess3(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_hess4(WlsObject *b_wls,
                          const ::coder::array<double, 2U> &quad_pnts,
                          const ::coder::array<double, 2U> &varargin_1,
                          ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_lap1(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         const ::coder::array<double, 2U> &varargin_2,
                         int varargin_3,
                         ::coder::array<double, 2U> &varargout_1,
                         ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_lap2(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         const ::coder::array<double, 2U> &varargin_2,
                         ::coder::array<double, 2U> &varargout_1,
                         ::coder::array<double, 2U> &varargout_2);

static inline void wls_var_lap3(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         ::coder::array<double, 2U> &varargout_1);

static inline void wls_var_lap4(WlsObject *b_wls,
                         const ::coder::array<double, 2U> &quad_pnts,
                         const ::coder::array<double, 2U> &varargin_1,
                         ::coder::array<double, 2U> &varargout_1);

} // namespace wls

#include "wls_internal.cpp"
#endif
// End of code generation (wls_internal.h)
