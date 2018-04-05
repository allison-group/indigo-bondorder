SET(INDIGOX_PYTHON_INTERFACE_SRCS
    interface.cpp
    pyatom.cpp
    pybond.cpp
    pyelement.cpp
    pymolecule.cpp
    pyoptions.cpp
    pyperiodictable.cpp)

# Setup pybind11 properly
SET(PYBIND11_PYTHON_VERSION 3.4)
ADD_SUBDIRECTORY(${PROJECT_SOURCE_DIR}/external/pybind11 ${PROJECT_BINARY_DIR}/pybind11)
# Find the site-packages location to install library to
EXECUTE_PROCESS( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
                 OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
IF(PYTHON_SITE_PACKAGES)
  MESSAGE(STATUS "Found Python site-packages location: ${PYTHON_SITE_PACKAGES}")
ELSE(PYTHON_SITE_PACKAGES)
  MESSAGE(WARNING "Unable to locate Python site-packages location. Python module will be installed in ${CMAKE_INSTALL_FULL_LIBDIR}")
  SET(PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_FULL_LIBDIR})
ENDIF(PYTHON_SITE_PACKAGES)

# Add the python module
PYBIND11_ADD_MODULE(pyindigox ${INDIGOX_PYTHON_INTERFACE_SRCS})
# Check compiler has required c++-14 features
# GCC versions < 5 pass the has -std=c++14 flag test, but don't have relaxed_constexpr as required
TARGET_COMPILE_FEATURES(pyindigox PRIVATE cxx_return_type_deduction)
TARGET_COMPILE_FEATURES(pyindigox PRIVATE cxx_relaxed_constexpr)
TARGET_LINK_LIBRARIES(pyindigox PRIVATE indigox)

# Set install rules
INSTALL(TARGETS pyindigox DESTINATION ${PYTHON_SITE_PACKAGES}/indigox)
INSTALL(FILES __init__.py DESTINATION ${PYTHON_SITE_PACKAGES}/indigox)
