cmake_minimum_required(VERSION 3.10)

project(JustInterp VERSION 0.0.1 LANGUAGES CXX)

enable_testing()

set(CMAKE_CXX_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
	message("Using default value of CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
else()
	message("Obtained CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

if (MSVC)
	message("Building on MSVC: ${MSVC_VERSION}")
	set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
else()
	message("Building on Unix")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -Wall -pedantic")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -Wall -pedantic")

	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/Installation/${CMAKE_BUILD_TYPE}")
endif()

add_subdirectory(src)
add_subdirectory(test)