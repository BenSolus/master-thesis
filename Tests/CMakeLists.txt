
################################################################################
# Find Boosts unit test framework                                              #
################################################################################

FIND_PACKAGE(Boost  COMPONENTS unit_test_framework REQUIRED)

################################################################################
# Add test_*.cpp files to ctest                                                #
################################################################################

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

  ADD_TEST(${name} "${PROJECT_BINARY_DIR}/Bin/${exe}")
ENDFOREACH(source ${SOURCES})
