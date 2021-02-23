cmake_minimum_required(VERSION 3.10)

add_library(JustInterp INTERFACE)

target_include_directories(JustInterp INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/include
                           $<INSTALL_INTERFACE:include>)