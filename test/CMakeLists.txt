cmake_minimum_required(VERSION 3.10)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(TestRunner ${SOURCES})

target_link_libraries(TestRunner JustInterp)
target_include_directories(TestRunner PUBLIC "thirdparty/doctest")

install(TARGETS TestRunner DESTINATION test)

include(thirdparty/doctest/scripts/doctest.cmake)
doctest_discover_tests(TestRunner)