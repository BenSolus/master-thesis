
################################################################################
# Suppress warnings                                                            #
################################################################################

CMAKE_POLICY(SET CMP0012 NEW)

################################################################################
# Find Doxygen                                                                 #
################################################################################

FIND_PACKAGE(Doxygen)

IF(${DOXYGEN_FOUND})

################################################################################
# Configure Doxygen                                                            #
################################################################################

  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
                 ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                 @ONLY)

################################################################################
# Add documentation build option                                               #
################################################################################

  ADD_CUSTOM_TARGET(doc
                    COMMAND ${DOXYGEN_EXECUTABLE}
                            ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                    COMMENT "Generating API documentation with Doxygen."
                    VERBATIM)
ENDIF()
