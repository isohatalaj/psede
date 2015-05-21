
set(LIBPSEDE_SOURCE cps.c cps_ode.c)
set(LIBPSEDE_INCLUDE cps.h cps_ode.h)

add_library(psede ${LIBPSEDE_SOURCE})
install(TARGETS psede DESTINATION lib)
install(FILES ${LIBPSEDE_INCLUDE} DESTINATION include)
