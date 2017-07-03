cmake_minimum_required(VERSION 2.8.4)
project(psede)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(fftw3)
find_package(LAPACK)

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
endif()

add_subdirectory(src)
add_subdirectory(tests)

