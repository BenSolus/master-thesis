CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(MasterThesis)

SET(CMAKE_VERBOSE_MAKEFILE OFF)

################################################################################
# Policy for relative paths                                                    #
################################################################################

CMAKE_POLICY(SET CMP0015 NEW)

SET(H2INC  ""    CACHE PATH "Path to the include directory of the H2Lib.")
SET(H2LINK ""    CACHE PATH "Path to the linking directory of the H2Lib.")
SET(DOC    OFF   CACHE BOOL "Whether to build the documentation")

################################################################################
# Setup for various paths                                                      #
################################################################################

SET(CMAKE_MODULE_PATH      "${PROJECT_SOURCE_DIR}/Moduls")
SET(LIBRARY_OUTPUT_PATH    "${PROJECT_BINARY_DIR}/Lib")
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/Bin")

################################################################################
# Find packages                                                                #
################################################################################

FIND_PACKAGE(LAPACK REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
FIND_PACKAGE(OpenCL REQUIRED)
FIND_PACKAGE(Cairo)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  FIND_PACKAGE(GTest REQUIRED)
ENDIF()

################################################################################
# Adding various include directories                                           #
################################################################################

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/Library)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/Tests)
INCLUDE_DIRECTORIES(${H2INC})

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
ENDIF()

################################################################################
# Adding various linking directories                                           #
################################################################################

LINK_DIRECTORIES(${H2LINK})

################################################################################
# Default compiler flags                                                       #
################################################################################

SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -std=c11 \
                                        -Wall \
                                        -Wextra \
                                        -pedantic \
                                        -Wno-overlength-strings \
                                        -pipe \
                                        -mavx2 \
                                        -DUSE_FLOAT=1 \
                                        -DUSE_OPENMP=1 \
                                        -DUSE_OPENCL=1 \
                                        ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 \
                                        -Wall \
                                        -Wextra \
                                        -pedantic \
                                        -Wno-overlength-strings \
                                        -pipe \
                                        -mavx2 \
                                        -DUSE_FLOAT=1 \
                                        -DUSE_OPENMP=1 \
                                        -DUSE_OPENCL=1 \
                                        ${OpenMP_CXX_FLAGS}")

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} \
                            ${OpenMP_EXE_LINKER_FLAGS}")

SET(CMAKE_C_FLAGS_RELEASE    "${CMAKE_C_FLAGS_RELEASE} -O0 -g")
SET(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -O0 -g")

################################################################################
# Debug compiler flags                                                         #
################################################################################

SET(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   -O0")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS_DEBUG    "${CMAKE_C_FLAGS_DEBUG}   -fprofile-arcs \
                                                       -ftest-coverage")
  SET(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs \
                                                       -ftest-coverage")
ENDIF()

################################################################################
# Adding various sub directories for compilation                               #
################################################################################

ADD_SUBDIRECTORY(Library)

IF(DOC MATCHES ON)
  ADD_SUBDIRECTORY(Doc)
ENDIF()

ADD_SUBDIRECTORY(Tests)
