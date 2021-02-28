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


set(JUST_INTERP_TARGET_NAME ${PROJECT_NAME})
set(JUST_INTERP_INCLUDE_BUILD_DIR "${PROJECT_SOURCE_DIR}/include")


add_library(${JUST_INTERP_TARGET_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${JUST_INTERP_TARGET_NAME} ALIAS ${JUST_INTERP_TARGET_NAME})
target_include_directories(
	${JUST_INTERP_TARGET_NAME}
	INTERFACE
    $<BUILD_INTERFACE:${JUST_INTERP_INCLUDE_BUILD_DIR}>
	$<INSTALL_INTERFACE:include>
)


add_subdirectory(test)