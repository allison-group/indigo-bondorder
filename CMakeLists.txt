# Should figure out minimum needed at some point. 3.5 is rather high
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)
CMAKE_POLICY(VERSION 3.5)

PROJECT(indigox VERSION 0.7.2.0)

# We require c++-14 compiler
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_EXTENSIONS OFF)
SET(CMAKE_CXX_VISIBILITY_PRESET hidden)

# Set defaults for available options
OPTION(BUILD_PYTHON "Build Python bindings" ON)
OPTION(USE_JAVA_QUICKBB "Build with Java interface to libtw so can use QuickBB tree-decomposition method" OFF)
OPTION(BUILD_EXAMPLES "Build a few example programs" ON)

# Default to release mode build. Does Debug actually do anything?
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
  MESSAGE(STATUS "Defaulting to Release mode build. Use -DCMAKE_BUILD_TYPE=Debug to build in debug mode, though it doesn't really do anything currently.")
ENDIF()

# Setup some compile flags, for all versions and specific versions
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-deprecated-declarations -Wno-deprecated")
IF("${CMAKE_SYSTEM_NAME}" STREQL "Linux")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF()
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -DDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/indigox)

# Set default install prefix to /usr/local/opt/indigox
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "/usr/local/opt/${PROJECT_NAME}" CACHE PATH "..." FORCE)
ELSE()
  SET(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME} CACHE PATH "..." FORCE)
ENDIF()
# Setup install directories
INCLUDE(GNUInstallDirs)
ADD_DEFINITIONS(-DDFT_DATA_DIR=${CMAKE_INSTALL_FULL_DATAROOTDIR})

# List all the files to compile into the indigox library
SET(INDIGOX_LIB_SRCS
    src/algorithm/electron_optimisation.cpp
    src/algorithm/formalbonds/astar.cpp
    src/algorithm/formalbonds/elimination_ordering.cpp
    src/algorithm/formalbonds/electron_optimisation_algorithm.cpp
    src/algorithm/formalbonds/fpt.cpp
    src/algorithm/formalbonds/local_optimisation.cpp
    src/classes/atom.cpp
    src/classes/bond.cpp
    src/classes/electron_graph.cpp
    src/classes/iters.cpp
    src/classes/molecular_graph.cpp
    src/classes/molecule.cpp
    src/classes/nicetreedecomp.cpp
    src/classes/periodictable.cpp
    src/classes/permutablegraph.cpp
    src/classes/treedecomp.cpp
    src/utils/common.cpp
    src/utils/filereader.cpp
    src/utils/helpers.cpp
    src/utils/java_interface.cpp
    src/utils/options.cpp)

# Build the library
ADD_LIBRARY(indigox STATIC ${INDIGOX_LIB_SRCS})

# We don't provide a boost version and it is required, so find it
# Sometime I need to figure out a minimum required version for it. :S
FIND_PACKAGE(Boost 1.66.0)
IF(NOT Boost_FOUND)
  MESSAGE(STATUS "Using provided Boost version.")
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external/boost/include)
  INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/external/boost/include/boost 
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
ELSE()
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

# Add interface to java if it's requested.
IF(USE_JAVA_QUICKBB)
  FIND_PACKAGE(JNI REQUIRED)
  ADD_DEFINITIONS(-DBUILD_JAVA)
  INCLUDE_DIRECTORIES(${JNI_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(indigox ${JNI_LIBRARIES})
ENDIF(USE_JAVA_QUICKBB)

# Build bindings if requested
IF(BUILD_PYTHON)
  ADD_SUBDIRECTORY(src/python)
ELSE(BUILD_PYTHON)
  MESSAGE(STATUS "Building of Python bindings disabled.")
ENDIF(BUILD_PYTHON)

IF(BUILD_EXAMPLES)
    ADD_EXECUTABLE(cxx_example examples/cxx_example.cpp)
    TARGET_INCLUDE_DIRECTORIES(cxx_example BEFORE PRIVATE ${PROJECT_SOURCE_DIR}/include)
    TARGET_LINK_LIBRARIES(cxx_example indigox)
    ADD_CUSTOM_TARGET(examples ALL)
    ADD_DEPENDENCIES(examples cxx_example)
    INSTALL(TARGETS cxx_example DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
    INSTALL(FILES examples/cxx_example.cpp examples/python_example.py
            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples)
ENDIF()

# Build documentation
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  SET(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
  SET(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  CONFIGURE_FILE(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  ADD_CUSTOM_TARGET(doc 
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with doxygen"
  VERBATIM )
  MESSAGE(STATUS "Doxygen located. Use 'make doc' to build documentation.")
ELSE(DOXYGEN_FOUND)
  MESSAGE(STATUS "Unable to locate doxygen. Documentation will NOT be built.")
ENDIF(DOXYGEN_FOUND)

FILE(COPY data DESTINATION .)
FILE(COPY external/libtw.jar DESTINATION external)
FILE(COPY examples/python_example.py DESTINATION examples)

# Install examples
# Install library
INSTALL(TARGETS indigox DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install header files. This should includes boost if the repo version is being used
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include/indigox
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install data files
INSTALL(FILES data/atomenergies_def2svpd.txt
        data/atomenergies_def2tzvppd.txt
        data/bondenergies_def2svpd.txt
        data/bondenergies_def2svpd_withbsse.txt
        data/bondenergies_def2tzvppd.txt
        data/bondenergies_def2tzvppd_withbsse.txt
        data/periodictable.txt
        external/libtw.jar
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
# Install documentation
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION ${CMAKE_INSTALL_DOCDIR} OPTIONAL)

# Uninstall target
IF(NOT TARGET uninstall)
  CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                  IMMEDIATE @ONLY)
  ADD_CUSTOM_TARGET(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake VERBATIM)
ENDIF()
