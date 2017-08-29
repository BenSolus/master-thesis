[![Travis.ci Shield](https://img.shields.io/travis/BenSolus/master-thesis/master.svg?style=plastic&label=Linux)](https://travis-ci.org/BenSolus/master-thesis)

# My master thesis

This is the git repository for the source code of my master thesis *"An
efficient solver for the Green cross approximation method on GPUs"*
based on [Approximation of integral operators by Green quadrature and nested
cross approximation](https://link.springer.com/article/10.1007/s00211-015-0757-y) and
building on the implementation of the [Green cross approximation for boundary element methods](https://arxiv.org/abs/1510.07244) for the
[H2Lib](http://www.h2lib.org).

# Dependencies

*   [Cairo](https://cairographics.org/) (1.14.10+)
*   [H2Lib](http://www.h2lib.org/) (3.0+)
*   [LAPACK](http://www.netlib.org/lapack/) (3.7.0+)
*   [OpenMP](http://www.openmp.org/) (gcc 7.1.1+)
*   [OpenCL](https://www.khronos.org/opencl/) (NVIDIA 1.2+)

Older versions or other vendors might also work, but untested.

# Compiling

Following bash script should retrieve the source and compile it in debug mode:

```shell
#!/usr/bin/bash

git clone https://github.com/BenSolus/master-thesis.git

cd master-thesis || exit 1

PROJECT_DIR=${PWD}

mkdir -p "${PROJECT_DIR}/builds/debug"

cd "${PROJECT_DIR}/builds/debug" || exit 1

cmake -DCMAKE_BUILD_TYPE=Debug -DH2INC:PATH=<H2Lib/Header> -DH2LINK:PATH=<H2Lib/Library/Object> -DDOC:BOOL=ON ../..

make -j2

```

where <H2Lib/Header> needs to be replaced with the path to the H2Lib header
files and <H2Lib/Library/Object> needs to be replaced with the path containing
the library object file (libh2.{a/so}). All executables reside in
${PROJECT_DIR}/builds/master-thesis-debug/Bin
