#pragma once

#ifdef _MSC_VER
    #define FORCE_INLINE __forceinline
#else
    #define FORCE_INLINE inline __attribute__((always_inline))
#endif

#if defined(_MSC_VER)
    #define FORCE_UNROLL __pragma(loop(ivdep)) __pragma(unroll)
#elif defined(__clang__) || defined(__GNUC__)
    #define DO_PRAGMA(x) _Pragma(#x)
    #define FORCE_UNROLL DO_PRAGMA(unroll)
#else
    #define FORCE_UNROLL
#endif