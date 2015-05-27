
set(LIBPSEDE_SOURCE psede.c fct.c diff.c ode.c)
set(LIBPSEDE_INCLUDE psede.h psede_fct.h psede_diff.h psede_ode.h)

add_library(psede ${LIBPSEDE_SOURCE})
install(TARGETS psede DESTINATION lib)
install(FILES ${LIBPSEDE_INCLUDE} DESTINATION include)
