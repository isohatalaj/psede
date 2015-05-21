# Find the FFTW3 library

include(LibFindMacros)

# libfind_package(fftw3)
libfind_pkg_check_modules(fftw3_PKGCONF)

find_path(fftw3_INCLUDE_DIR
  NAMES fftw3.h
  PATHS ${fftw3_PKGCONF_INCLUDE_DIRS}
)

find_library(fftw3_LIBRARY
  NAMES fftw3
  PATHS ${fftw3_PKGCONF_LIBRARY_DIRS}
)

set(fftw3_PROCESS_INCLUDES fftw3_INCLUDE_DIR)
set(fftw3_PROCESS_LIBS fftw3_LIBRARY)
libfind_process(fftw3)

