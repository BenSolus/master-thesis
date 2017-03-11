CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(MasterThesis)

SET(H2INC  ""    CACHE PATH "Path to the include directory of the H2Lib.")
SET(H2LINK ""    CACHE PATH "Path to the linking directory of the H2Lib.")
SET(DOC    OFF   CACHE BOOL "Whether to build the documentation")

################################################################################
# Setup for various paths                                                      #
################################################################################

SET(CMAKE_MODULE_PATH      "${PROJECT_SOURCE_DIR}/Moduls")
SET(LIBRARY_OUTPUT_PATH    "${PROJECT_BINARY_DIR}/lib")
SET(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin")

################################################################################
# Find packages                                                                #
################################################################################

FIND_PACKAGE(OpenCL)

################################################################################
# Adding various include directories                                           #
################################################################################

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/Library)
INCLUDE_DIRECTORIES(${H2INC})

################################################################################
# Adding various linking directories                                           #
################################################################################

LINK_DIRECTORIES(${H2LINK})

################################################################################
# Default compiler flags                                                       #
################################################################################

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -Wall -Wextra -pedantic -pipe")

################################################################################
# Debug compiler flags                                                         #
################################################################################

SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} \
                         -O0 \
                         -fprofile-arcs \
                         -ftest-coverage")

################################################################################
# Adding various sub directories for compilation                               #
################################################################################

ADD_SUBDIRECTORY(Library)

IF(DOC MATCHES ON)
  ADD_SUBDIRECTORY(Doc)
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  ENABLE_TESTING()

  ADD_SUBDIRECTORY(Tests)
ENDIF()
