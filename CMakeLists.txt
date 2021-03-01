cmake_minimum_required(VERSION 3.10)

project(
    JustInterp
    VERSION 0.0.1
    LANGUAGES CXX
)

#================================================================
# Check if JustInterp is built as subproject or it's main project
#================================================================
set(JUST_INTERP_MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(JUST_INTERP_MAIN_PROJECT ON)
endif()

#================================================================
# Define options
#================================================================
if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

option(JUST_INTERP_BUILD_TESTS "Build the unit tests" ${JUST_INTERP_MAIN_PROJECT})

#================================================================
# Define configuration
#================================================================
set(JUST_INTERP_TARGET_NAME       ${PROJECT_NAME})
set(JUST_INTERP_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/include")

#================================================================
# Default build type
#================================================================
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("    -- Using default value of CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
else()
    message("    -- Obtained CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

#================================================================
# Add INTERFACE library
#================================================================
add_library(${JUST_INTERP_TARGET_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${JUST_INTERP_TARGET_NAME} ALIAS ${JUST_INTERP_TARGET_NAME})
target_compile_features(${JUST_INTERP_TARGET_NAME} INTERFACE cxx_std_17)
target_include_directories(
    ${JUST_INTERP_TARGET_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${JUST_INTERP_INCLUDE_BUILD_DIR}>
    $<INSTALL_INTERFACE:include>
)

#================================================================
# Buld unit tests
#================================================================
if (JUST_INTERP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()