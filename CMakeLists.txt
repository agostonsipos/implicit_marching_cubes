cmake_minimum_required(VERSION 3.1.0)

project(marching_cubes)

if(CMAKE_VERSION VERSION_LESS "3.7.0")
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif()

if(UNIX)
        add_compile_options(-fPIC)
endif()

set(Transfinite_DIR "../transfinite/src") # https://github.com/salvipeter/transfinite
include_directories("${Transfinite_DIR}/geom")

if(WIN32)
        link_directories("${Transfinite_DIR}/geom/Release")
else()
        link_directories("${Transfinite_DIR}/../release/geom")
endif()

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


if(WIN32)
	set(CompilerFlags
			CMAKE_CXX_FLAGS
			CMAKE_CXX_FLAGS_DEBUG
			CMAKE_CXX_FLAGS_RELEASE
			CMAKE_C_FLAGS
			CMAKE_C_FLAGS_DEBUG
			CMAKE_C_FLAGS_RELEASE
			)
	foreach(CompilerFlag ${CompilerFlags})
	  string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
	endforeach()
endif()
