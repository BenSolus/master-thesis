
################################################################################
# Find Boosts unit test framework                                              #
################################################################################

FIND_PACKAGE(Boost  COMPONENTS unit_test_framework program_options REQUIRED)

################################################################################
# Add test_*.cpp files to ctest                                                #
################################################################################

IF(CMAKE_BUILD_TYPE MATCHES Debug)
  ENABLE_TESTING()

  FILE(GLOB SOURCES RELATIVE "${PROJECT_SOURCE_DIR}/Tests"
                             "${PROJECT_SOURCE_DIR}/Tests/test_*.cpp")

  FOREACH(source ${SOURCES})
    STRING(REPLACE ".cpp" "" exe ${source})
    STRING(REPLACE "test_" "" name ${exe})

    ADD_EXECUTABLE(${exe} ${source})

    TARGET_LINK_LIBRARIES(${exe}
                          ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                          ${OpenCL_LIBRARIES}
                          MasterThesis
                          h2)
    IF(LAPACK_FOUND)
      TARGET_LINK_LIBRARIES(${exe} ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
    ENDIF()

    IF(CAIRO_FOUND)
      TARGET_LINK_LIBRARIES(${exe} ${CAIRO_LIBRARIES})
    ENDIF()

    ADD_TEST(${name} "${PROJECT_BINARY_DIR}/Bin/${exe}")
  ENDFOREACH(source ${SOURCES})

ENDIF()

ADD_EXECUTABLE(benchmark benchmark.cpp)

TARGET_LINK_LIBRARIES(benchmark ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
                                ${OpenCL_LIBRARIES}
                                MasterThesis
                                h2
                                ${CAIRO_LIBRARIES}
                                ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
                                ${OPENGL_LIBRARIES}
                                ${GLUT_LIBRARIES})