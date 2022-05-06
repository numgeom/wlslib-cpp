/* Copyright 2022 The NumGeom Group, Stony Brook University. */

// Thread-safe wrappers for some I/O functions

#pragma once

#ifndef _m2c_lib_h
#define _m2c_lib_h

#ifdef MATLAB_MEX_FILE
#include "m2c_mex.h"

#else

#include <stdexcept>

#ifndef M2C_SKIP_THREAD_ID
#define M2C_PRINT_THREAD_ID 1
#else
#define M2C_PRINT_THREAD_ID 0
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

#ifndef m2cPrintf
#define m2cPrintf(...) std::fprintf(stdout, __VA_ARGS__)
#endif

#ifndef m2cPrintErrorNoflush
#define m2cPrintErrorNoflush(...) std::fprintf(stderr, __VA_ARGS__)
#endif

#ifndef m2cFlushError
#define m2cFlushError() std::fflush(stderr)
#endif

#ifndef m2cPrintError
#define m2cPrintError(...) { \
    m2cPrintErrorNoflush(__VA_ARGS__); \
    m2cFlushError(); \
}
#endif

#ifdef _OPENMP
#define m2cWarnMsgTxt(message) \
    if (!omp_in_parallel()) { \
        std::fprintf(stderr, "%s\n", message); \
        m2cFlushError(); \
    } else { \
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                m2cPrintErrorNoflush("Thread %d: ", omp_get_thread_num()); \
            m2cPrintErrorNoflush("Warning: "); \
            m2cPrintErrorNoflush("%s\n", message); \
            m2cFlushError(); \
        } \
    }
#else
#define m2cWarnMsgTxt(message) { \
    m2cPrintErrorNoflush("%s\n", message); \
    m2cFlushError(); \
}
#endif // _OPENMP


#ifdef _OPENMP
#define m2cWarnMsgIdAndTxt(warnId, ...) \
    if (!omp_in_parallel()) { \
        m2cPrintErrorNoflush("%s: ", warnId); \
        m2cPrintErrorNoflush(__VA_ARGS__); \
        m2cPrintErrorNoflush("\n"); \
        m2cFlushError(); \
    } else {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                m2cPrintErrorNoflush("Thread %d: ", omp_get_thread_num()); \
            m2cPrintErrorNoflush("%s: ", warnId); \
            m2cPrintErrorNoflush(__VA_ARGS__); \
            m2cPrintErrorNoflush("\n"); \
            m2cFlushError(); \
        } \
    }
#else
#define m2cWarnMsgIdAndTxt(warnId, ...) { \
    m2cPrintErrorNoflush("%s: ", warnId); \
    m2cPrintErrorNoflush(__VA_ARGS__); \
    m2cPrintErrorNoflush("\n"); \
    m2cFlushError(); \
}
#endif // _OPENMP

#ifdef _OPENMP
#define m2cErrMsgTxt(message) { \
    if (!omp_in_parallel()) { \
        m2cPrintErrorNoflush("%s\n", message); \
        m2cFlushError(); \
    } else {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                m2cPrintErrorNoflush("Thread %d: ", omp_get_thread_num()); \
            m2cPrintErrorNoflush("Error: "); \
            m2cPrintErrorNoflush("%s\n", message); \
            m2cFlushError(); \
        } \
    } \
    throw std::runtime_error("unnamedRuntimeError"); \
}
#else
#define m2cErrMsgTxt(message) { \
    m2cPrintErrorNoflush("%s\n", message); \
    m2cFlushError(); \
    throw std::runtime_error("unnamedRuntimeError"); \
}
#endif // _OPENMP


#ifdef _OPENMP
#define m2cErrMsgIdAndTxt(errId, ...) { \
    if (!omp_in_parallel()) {\
        m2cPrintErrorNoflush("%s: ", errId); \
        m2cPrintErrorNoflush(__VA_ARGS__); \
        m2cPrintErrorNoflush("\n"); \
        m2cFlushError(); \
    } else { \
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                m2cPrintErrorNoflush("Thread %d: ", omp_get_thread_num()); \
            m2cPrintErrorNoflush("%s: ", errId); \
            m2cPrintErrorNoflush(__VA_ARGS__); \
            m2cPrintErrorNoflush("\n"); \
            m2cFlushError(); \
        } \
    } \
    throw std::runtime_error(errId); \
}
#else
#define m2cErrMsgIdAndTxt(errId, ...) { \
    m2cPrintErrorNoflush("%s: ", errId); \
    m2cPrintErrorNoflush(__VA_ARGS__); \
    m2cPrintErrorNoflush("\n"); \
    m2cFlushError(); \
    throw std::runtime_error(errId); \
}
#endif // _OPENMP


#ifndef NDEBUG
#ifdef _OPENMP
#define m2cAssert(cond, message) { \
    if (!omp_in_parallel() && !(cond)) {\
        m2cPrintErrorNoflush("Assertion failed (%s) at line %d of file '%s'. %s\n", \
            #cond, __LINE__, __FILE__, message);\
        m2cFlushError(); \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message);\
    } else if (!(cond)) {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                m2cPrintErrorNoflush("Thread %d: ", omp_get_thread_num()); \
            m2cPrintErrorNoflush("Assertion failed (%s) at line %d of file '%s'. %s\n", \
                #cond, __LINE__, __FILE__, message); \
            m2cFlushError(); \
        } \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message);\
    } \
}
#else // _OPENMP
#define m2cAssert(cond, message) { \
    if (!(cond)) { \
        m2cPrintErrorNoflush("Assertion failed (%s) at line %d of file '%s'. %s\n", \
            #cond, __LINE__, __FILE__, message); \
        m2cFlushError(); \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message); \
    } \
}
#endif // _OPENMP
#else // NDEBUG
#define m2cAssert(cond, message) ((void)0)
#endif // NDEBUG

#endif // MATLAB_MEX_FILE
#endif // _m2c_lib_h
