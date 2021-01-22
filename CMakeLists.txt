cmake_minimum_required(VERSION 3.1.0)

project(marching_cubes)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(GEOM_DIR "libgeom") # https://github.com/salvipeter/libgeom
include_directories("${GEOM_DIR}")

link_directories("${GEOM_DIR}/build")


add_library(marching_cubes STATIC
    mc.h
    mc.cc
    tables.cc
)

target_link_libraries(marching_cubes "geom")

add_executable(marching_cubes_test
    test.cc
)

target_link_libraries(marching_cubes_test marching_cubes)
