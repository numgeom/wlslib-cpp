# Weighted Least Squares (WLS) Library in C++ #

This project implements the kernels for performing WLS computations by taking advantage of LAPACK optimiations. This library is template-baed header-only; just include `wls.hpp` and you are ready to go!

## Fortran name mangling and LAPACK integer customization ##

In `wls.hpp`, controling Fortran name mangling is done via a preprocessor macro `WLS_FC`, which takes values from one to six:

1. all lower cases without appending underscore,
2. all lower cases followed by a single underscore (default),
3. all lower cases followed by two underscores,
4. all upper cases without appending underscore,
5. all upper cases followed by a single underscore, and
6. all upper cases followed by two single underscores.

Example compiling command for appending a single understore and converting to all lower cases is

```cpp
g++ -c -DWLS_FC=2 proj.cpp
```

Also, depending on LAPACK configurations, the size of integer type may vary. For instance, MATLAB comes with a customized installation with integer type of `mwSignedIndex`, which is `ptrdiff_t` (i.e., 64-bit signed integer on 64-bit machines). `wls.hpp` takes an optional macro of `WLS_LAPACK_INT` that allows the user to customize the integer type, e.g.,

```cpp
g++ -c -DWLS_LAPACK_INT=mwSignedIndex proj.cpp
```

Notice that the default value of `WLS_LAPACK_INT` is `int.`

## Copyright ##

Copyright (C) 2020 NumGeom Group at Stony Brook University

`wls.hpp` is distributed under GNU Lesser General Public License v3 (LGPLv3). Refer to [LICENSE](./LICENSE) for more details.
