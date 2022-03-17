/* Copyright 2022 The NumGeom Group, Stony Brook University. */

// Thread-safe wrappers for some I/O functions

#pragma once

#ifndef _m2c_lib_h
#define _m2c_lib_h

#ifdef MATLAB_MEX_FILE
#include "m2c_mex.h"
#else

#ifndef M2C_SKIP_THREAD_ID
#define M2C_PRINT_THREAD_ID 1
#else
#define M2C_PRINT_THREAD_ID 0
#endif

#ifdef _OPENMP
#include "omp.h"
#endif

#define m2cPrintf(...) \
    std::fprintf(stdout, __VA_ARGS__)

#define m2cPrintError(...) { \
        std::fprintf(stderr, __VA_ARGS__); \
        std::fflush(stderr); \
    }

#ifdef _OPENMP
#define m2cWarnMsgTxt(message) \
    if (!omp_in_parallel()) { \
        std::fprintf(stderr, "%s\n", message); \
        std::fflush(stderr); \
    } else { \
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                std::fprintf(stderr, "Thread %d: ", omp_get_thread_num()); \
            std::fprintf(stderr, "Warning: "); \
            std::fprintf(stderr, "%s\n", message); \
            std::fflush(stderr); \
        } \
    }
#else
#define m2cWarnMsgTxt(message) { \
    std::fprintf(stderr, "%s\n", message); \
    std::fflush(stderr); \
}
#endif // _OPENMP


#ifdef _OPENMP
#define m2cWarnMsgIdAndTxt(warnId, ...) \
    if (!omp_in_parallel()) { \
        std::fprintf(stderr, "%s: ", warnId); \
        std::fprintf(stderr, __VA_ARGS__); \
        std::fprintf(stderr, "\n"); \
        std::fflush(stderr); \
    } else {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                std::fprintf(stderr, "Thread %d: ", omp_get_thread_num()); \
            std::fprintf(stderr, "%s: ", warnId); \
            std::fprintf(stderr, __VA_ARGS__); \
            std::fprintf(stderr, "\n"); \
            std::fflush(stderr); \
        } \
    }
#else
#define m2cWarnMsgIdAndTxt(warnId, ...) { \
    std::fprintf(stderr, "%s: ", warnId); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fprintf(stderr, "\n"); \
    std::fflush(stderr); \
}
#endif // _OPENMP

#ifdef _OPENMP
#define m2cErrMsgTxt(message) { \
    if (!omp_in_parallel()) { \
        std::fprintf(stderr, "%s\n", message); \
        std::fflush(stderr); \
    } else {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                std::fprintf(stderr, "Thread %d: ", omp_get_thread_num()); \
            std::fprintf(stderr, "Error: "); \
            std::fprintf(stderr, "%s\n", message); \
            std::fflush(stderr); \
        } \
    } \
    throw std::runtime_error("unnamedRuntimeError"); \
}
#else
#define m2cErrMsgTxt(message) { \
    std::fprintf(stderr, "%s\n", message); \
    std::fflush(stderr); \
    throw std::runtime_error("unnamedRuntimeError"); \
}
#endif // _OPENMP


#ifdef _OPENMP
#define m2cErrMsgIdAndTxt(errId, ...) { \
    if (!omp_in_parallel()) {\
        std::fprintf(stderr, "%s: ", errId); \
        std::fprintf(stderr, __VA_ARGS__); \
        std::fprintf(stderr, "\n"); \
        std::fflush(stderr); \
    } else { \
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                std::fprintf(stderr, "Thread %d: ", omp_get_thread_num()); \
            std::fprintf(stderr, "%s: ", errId); \
            std::fprintf(stderr, __VA_ARGS__); \
            std::fprintf(stderr, "\n"); \
            std::fflush(stderr); \
        } \
    } \
    throw std::runtime_error(errId); \
}
#else
#define m2cErrMsgIdAndTxt(errId, ...) { \
    std::fprintf(stderr, "%s: ", errId); \
    std::fprintf(stderr, __VA_ARGS__); \
    std::fprintf(stderr, "\n"); \
    std::fflush(stderr); \
    throw std::runtime_error(errId); \
}
#endif // _OPENMP


#ifndef NDEBUG
#ifdef _OPENMP
#define m2cAssert(cond, message) { \
    if (!omp_in_parallel()) {\
        std::fprintf(stderr, "Assertion failed (%s) at line %d of file '%s'. %s\n", \
            #cond, __LINE__, __FILE__, message);\
        std::fflush(stderr); \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message);\
    } else if (!(cond)) {\
        _Pragma("omp critical") \
        { \
            if (M2C_PRINT_THREAD_ID) \
                std::fprintf(stderr, "Thread %d: ", omp_get_thread_num()); \
            std::fprintf(stderr, "Assertion failed (%s) at line %d of file '%s'. %s\n", \
                #cond, __LINE__, __FILE__, message); \
            std::fflush(stderr); \
        } \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message);\
    } \
}
#else // _OPENMP
#define m2cAssert(cond, message) { \
    if (!(cond)) { \
        std::fprintf(stderr, "Assertion failed (%s) at line %d of file '%s'. %s\n", \
            #cond, __LINE__, __FILE__, message); \
        std::fflush(stderr); \
        throw std::logic_error(std::string("Assertion failed (") + #cond + "). " + message); \
    }
}
#endif // _OPENMP
#else // NDEBUG
#define m2cAssert(cond, message) ((void)0)
#endif // NDEBUG

#endif // MATLAB_MEX_FILE
#endif // _m2c_lib_h
