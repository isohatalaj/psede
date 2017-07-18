

include_directories(${PROJECT_SOURCE_DIR}/src)
link_directories(${PROJECT_BINARY_DIR}/src)

add_executable(lintest lintest.c)
target_link_libraries(lintest psede ${MATH_LIBRARY} ${fftw3_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(difftest difftest.c)
target_link_libraries(difftest psede ${MATH_LIBRARY} ${fftw3_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(integtest integtest.c)
target_link_libraries(integtest psede ${MATH_LIBRARY} ${fftw3_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(kfe2dtest kfe2dtest.c)
target_link_libraries(kfe2dtest psede ${MATH_LIBRARY} ${fftw3_LIBRARIES} ${LAPACK_LIBRARIES})

