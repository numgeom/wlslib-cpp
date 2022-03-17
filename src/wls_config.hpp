/*
    This file is part of wlslib project

    Copyright (C) 2020 NumGeom Group at Stony Brook University

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software Foundation,
    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*!
 * @file wls_config.hpp
 *
 * @brief Config LAPACK and alignment keywords
 */

#ifndef WLS_CONFIG_HPP_
#define WLS_CONFIG_HPP_

#define __WLS_COMP_NULL 100  // give a positive large value

#ifdef __GNUC__
#  define WLS_COMP_FLAG 1
#  if __GNUC__ < 7 && !defined(_WIN32)
#    warning "To get best performance, GCC v7.0+ is recommended."
#  endif
#else
#  ifdef __INTEL_COMPILER
#    define WLS_COMP_FLAG 2
#  else
#    define WLS_COMP_FLAG __WLS_COMP_NULL
#  endif
#endif

// handle restrict keyword for memory aliasing
#if WLS_COMP_FLAG == 1
#  define WLS_RESTRICT __restrict__
#elif WLS_COMP_FLAG == 2
#  define WLS_RESTRICT __restrict
#else
#  define WLS_RESTRICT
#endif

// check alignment
#ifndef WLS_ALIGN
#  define WLS_ALIGN 32
#elif WLS_ALIGN <= 0
#  undef WLS_ALIGN
#  define WLS_ALIGN 32
#endif

#if !(((WLS_ALIGN - 1) & WLS_ALIGN) == 0 && WLS_ALIGN > 7)
#  error "Alignment must be power of 2 and greater than 7"
#endif

// get alignment attribute
#if WLS_COMP_FLAG != __WLS_COMP_NULL
#  if defined(_WIN32) && !defined(__GNUC__)
#    define WLS_ALIGN_ATTR __attribute__((align(WLS_ALIGN)))
#  else
// assume *nix
#    define WLS_ALIGN_ATTR __attribute__((aligned(WLS_ALIGN)))
#  endif
#else
#  define WLS_ALIGN_ATTR
#endif

// determine runtime alignment
#if WLS_COMP_FLAG == 1
#  define WLS_ASSERT_ALIGNED(__type, __ptr) \
    assert(ptrdiff_t(__ptr) % WLS_ALIGN == 0); \
    __ptr = (__type *)__builtin_assume_aligned(__ptr, WLS_ALIGN)
#elif WLS_COMP_FLAG == 2
#  define WLS_ASSERT_ALIGNED(__type, __ptr) __assume_aligned(__ptr, WLS_ALIGN)
#else
#  define WLS_ASSERT_ALIGNED(__type, __ptr) (void)0
#endif

// determine integer type used in LAPACK/BLAS
#ifdef MATLAB_MEX_FILE
#  ifdef WLS_LAPACK_INT
#    undef WLS_LAPACK_INT
#  endif
#  include <stddef.h>
#  define WLS_LAPACK_INT ptrdiff_t

#ifdef __MINGW32__
#  define WLS_FC 1  // on Windows, use lowercase
#else
#  define WLS_FC 2  // default is lower case with one "_" appended
#endif

#endif
#ifndef WLS_LAPACK_INT
#  ifdef MKL_INT
#    define WLS_LAPACK_INT MKL_INT
#  else
#    define WLS_LAPACK_INT int
#  endif
#endif  // !WLS_LAPACK_INT

typedef WLS_LAPACK_INT lapack_int;  ///< integer type for LAPACK kernels

// determine Fortran name mangling
//  1: lower
//  2: lower_
//  3: lower__
//  4: UPPER
//  5: UPPER_
//  6: UPPER__
#ifndef WLS_FC
#  define WLS_FC 2  // default is lower case with one "_" appended
#endif
#if WLS_FC == 1
#  undef WLS_FC
#  define WLS_FC(__l, __U) __l
#elif WLS_FC == 2
#  undef WLS_FC
#  define WLS_FC(__l, __U) __l##_
#elif WLS_FC == 3
#  undef WLS_FC
#  define WLS_FC(__l, __U) __l##__
#elif WLS_FC == 4
#  undef WLS_FC
#  define WLS_FC(__l, __U) __U
#elif WLS_FC == 5
#  undef WLS_FC
#  define WLS_FC(__l, __U) __U##_
#elif WLS_FC == 6
#  undef WLS_FC
#  define WLS_FC(__l, __U) __U##__
#else
#  error "Unknown WLS_FC option, must be in (1,2,3,4,5,6)"
#endif

#endif  // WLS_CONFIG_HPP_
