ADD_LIBRARY(MasterThesis
            greencross.c
            greencross.h
            laplace2d.c
            laplace2d.h
            laplace3d.c
            laplace3d.h
            ocl_system.c
            ocl_system.h
)

TARGET_LINK_LIBRARIES(MasterThesis ${GSL_LIBRARIES})

IF(CMAKE_EXPORT_COMPILE_COMMANDS EQUAL 1)
  ADD_CUSTOM_COMMAND(TARGET MasterThesis
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy
                             ${PROJECT_BINARY_DIR}/compile_commands.json
                             ${PROJECT_SOURCE_DIR}/compile_commands.json)
ENDIF()
