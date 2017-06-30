
include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(lintest lintest.c)
target_link_libraries(lintest psede ${MATH_LIBRARY} ${fftw3_LIBRARIES} ${LAPACK_LIBRARIES})

