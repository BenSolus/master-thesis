
################################################################################
# Add test_*.cpp files to ctest                                                #
################################################################################

IF(CMAKE_BUILD_TYPE MATCHES Debug)

  ENABLE_TESTING()

  FILE(GLOB SOURCES RELATIVE "${PROJECT_SOURCE_DIR}/Tests"
                             "${PROJECT_SOURCE_DIR}/Tests/test_*.hpp")

  MESSAGE(STATUS "${SOURCES}")

  ADD_EXECUTABLE(tests main_test.cpp ${SOURCES})

  SET_TARGET_PROPERTIES(tests PROPERTIES LINKER_LANGUAGE CXX)

  TARGET_LINK_LIBRARIES(tests ${GTEST_LIBRARIES} ${OpenCL_LIBRARIES} MasterThesis h2)

  IF(LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(tests ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES})
  ENDIF()

  IF(CAIRO_FOUND)
    TARGET_LINK_LIBRARIES(tests ${CAIRO_LIBRARIES})
  ENDIF()

  ADD_TEST(tests "${PROJECT_BINARY_DIR}/Bin/${exe}")
ENDIF()



ADD_EXECUTABLE(benchmark benchmark.cpp)

TARGET_LINK_LIBRARIES(benchmark ${OpenCL_LIBRARIES}
                                MasterThesis
                                h2
                                ${CAIRO_LIBRARIES}
                                ${LAPACK_LINKER_FLAGS} ${LAPACK_LIBRARIES}
                                ${OPENGL_LIBRARIES}
                                ${GLUT_LIBRARIES})