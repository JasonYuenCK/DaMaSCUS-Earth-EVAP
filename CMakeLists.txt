cmake_minimum_required(VERSION 3.12)

project(DaMaSCUS-SUN VERSION "0.1.1")

# Default to Release on single-config generators so local builds use optimization.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Keep older vendored dependencies configurable under newer CMake.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5 CACHE STRING "Minimum policy version for vendored deps")

# Use C++-11
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

set(BIN_DIR       ${PROJECT_SOURCE_DIR}/bin)
set(GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(INCLUDE_DIR   ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR       ${PROJECT_SOURCE_DIR}/lib)
set(SRC_DIR       ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR     ${PROJECT_SOURCE_DIR}/tests)
set(EXTERNAL_DIR  ${PROJECT_SOURCE_DIR}/external)

#External projects
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
find_package(Threads REQUIRED)

# Prefer Homebrew Boost on local macOS builds, but do not force that path on
# clusters where Boost is supplied by the module environment.
if(NOT DEFINED BOOST_ROOT AND NOT DEFINED Boost_ROOT AND EXISTS "/opt/homebrew/opt/boost@1.75")
  set(BOOST_ROOT "/opt/homebrew/opt/boost@1.75" CACHE PATH "Preferred Boost root")
  set(Boost_NO_SYSTEM_PATHS ON CACHE BOOL "Prefer BOOST_ROOT over system Boost")
endif()
find_package(Boost 1.65 REQUIRED)

include(FetchContent)
# obscura
set(OBSCURA_DIR   ${EXTERNAL_DIR}/obscura)
FetchContent_Declare(obscura
  GIT_REPOSITORY  https://github.com/temken/obscura.git
  GIT_TAG         v1.0.1
  SOURCE_DIR     "${OBSCURA_DIR}"
  CMAKE_ARGS     "-DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF" )
FetchContent_GetProperties(obscura)
if(NOT obscura_POPULATED)
    FetchContent_Populate(obscura)
    add_subdirectory( ${obscura_SOURCE_DIR} )
endif()

# Googletest
set(GTEST_DIR   ${EXTERNAL_DIR}/googletest)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.10.0
  SOURCE_DIR     "${GTEST_DIR}/src"
  BINARY_DIR     "${GTEST_DIR}/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      "" )

# Refresh the build-tree version header on every build.
set(VERSION_HEADER_TEMPLATE ${INCLUDE_DIR}/version.hpp.in)
set(VERSION_HEADER ${GENERATED_DIR}/version.hpp)
add_custom_target(generate_version_header ALL
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DTEMPLATE=${VERSION_HEADER_TEMPLATE}
    -DOUTPUT=${VERSION_HEADER}
    -DPROJECT_NAME=${PROJECT_NAME}
    -DPROJECT_VERSION=${PROJECT_VERSION}
    -DPROJECT_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    -DPROJECT_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    -DPROJECT_VERSION_PATCH=${PROJECT_VERSION_PATCH}
    -P ${PROJECT_SOURCE_DIR}/cmake/GenerateVersionHeader.cmake
  BYPRODUCTS ${VERSION_HEADER}
  COMMENT "Refreshing generated version header"
  VERBATIM)

# Source and include directories
include_directories( ${INCLUDE_DIR} ${MPI_INCLUDE_PATH})
add_subdirectory( ${SRC_DIR} )

# Code Coverage Configuration
if(NOT TARGET coverage_config)
  add_library(coverage_config INTERFACE)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

option (BUILD_TESTING "Build the testing tree." ON)
# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if (BUILD_TESTING AND (PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()
  FetchContent_GetProperties(googletest)
  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(
      ${googletest_SOURCE_DIR} 
      ${googletest_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif()
  add_subdirectory(${TESTS_DIR})
endif()
