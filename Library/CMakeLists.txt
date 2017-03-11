ADD_LIBRARY(MasterThesis
            objects.c
            objects.h
)

IF(CMAKE_EXPORT_COMPILE_COMMANDS EQUAL 1)
  ADD_CUSTOM_COMMAND(TARGET MasterThesis
                     POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy
                             ${PROJECT_BINARY_DIR}/compile_commands.json
                             ${PROJECT_SOURCE_DIR}/compile_commands.json)
ENDIF()
