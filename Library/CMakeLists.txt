ADD_LIBRARY(MasterThesis
            clamatrix.cl
            clavector.cl
            ie1d.c
            ie1d.h
            ocl_system.c
            ocl_system.h
)

IF(CMAKE_EXPORT_COMPILE_COMMANDS EQUAL 1)
  ADD_CUSTOM_COMMAND(TARGET MasterThesis
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy
                             ${PROJECT_BINARY_DIR}/compile_commands.json
                             ${PROJECT_SOURCE_DIR}/compile_commands.json)
ENDIF()
