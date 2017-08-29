ADD_LIBRARY(MasterThesis
            gcopencl.c
            gcopencl.h
            greencross.c
            greencross.h
            laplace2d.c
            laplace2d.h
            laplace3d.c
            laplace3d.h
            ocl_system.c
            ocl_system.h
  oclworkpkgs.h oclworkpkgs.c)

################################################################################
# Convert OpenCL source code into a string.                                    #
################################################################################

IF(WIN32)
  # Write windows script to perform conversion.
else(WIN32)
  ADD_CUSTOM_TARGET(OPENCL_SRC
                    COMMAND ${PROJECT_SOURCE_DIR}/.cmake.sh
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
ENDIF(WIN32)

ADD_DEPENDENCIES(MasterThesis OPENCL_SRC)

IF(CMAKE_EXPORT_COMPILE_COMMANDS EQUAL 1)
  ADD_CUSTOM_COMMAND(TARGET MasterThesis
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy
                             ${PROJECT_BINARY_DIR}/compile_commands.json
                             ${PROJECT_SOURCE_DIR}/compile_commands.json)
ENDIF()
