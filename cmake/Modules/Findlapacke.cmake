# Find the LAPACKE library

include(LibFindMacros)

# libfind_package(lapacke LAPACK)
libfind_pkg_check_modules(lapacke_PKGCONF)

find_path(lapacke_INCLUDE_DIR
  NAMES lapacke.h
  PATHS ${lapacke_PKGCONF_INCLUDE_DIRS}
)

find_library(lapacke_LIBRARY
  NAMES lapacke
  PATHS ${lapacke_PKGCONF_LIBRARY_DIRS}
)

set(lapacke_PROCESS_INCLUDES lapacke_INCLUDE_DIR)
set(lapacke_PROCESS_LIBS lapacke_LIBRARY)
libfind_process(lapacke)

