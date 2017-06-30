cmake_minimum_required(VERSION 2.8.4)
project(psede)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(fftw3)
find_package(lapack)

if (UNIX)
  find_library(MATH_LIBRARY m)
endif (UNIX)

add_subdirectory(src)
add_subdirectory(tests)

