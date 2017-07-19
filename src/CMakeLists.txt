
set(LIBPSEDE_SOURCE psede.c fct.c diff.c ode.c util.c func.c multi.c)
set(LIBPSEDE_INCLUDE psede.h psede_fct.h psede_diff.h psede_ode.h psede_util.h psede_func.h psede_multi.h)

add_library(psede ${LIBPSEDE_SOURCE})
install(TARGETS psede DESTINATION lib)
install(FILES ${LIBPSEDE_INCLUDE} DESTINATION include)
