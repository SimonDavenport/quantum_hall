##  This file generates the configure script for the entanglement spectrum program suite

# Process this file with autoconf to produce a configure script.

AC_INIT(Quantum Hall Entanglement and Monte Carlo programs, 1.0, simon.davenport2@gmail.com or ivan_ok2000@yahoo.com)
AC_CONFIG_SRCDIR(src/)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([subdir-objects])

# Option to use the mkl library
AC_ARG_ENABLE([intel_mkl],
  [AS_HELP_STRING([--enable-intel-mkl], [Set to yes to use intel mkl libray for cache efficient Lapack routines, if available])])

# Extra arguments to configure boost program options
AC_ARG_WITH([boost_include],
            [AC_HELP_STRING([--with-boost-include=DIR],
                            [boost program options include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([boost_lib],
            [AC_HELP_STRING([--with-boost-lib=DIR],
                            [boost program options lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure fftw 
AC_ARG_WITH([fftw_include],
            [AC_HELP_STRING([--with-fftw-include=DIR],
                            [fftw include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([fftw_lib],
            [AC_HELP_STRING([--with-fftw-lib=DIR],
                            [fftw lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure lapack
AC_ARG_WITH([lapack_lib],
            [AC_HELP_STRING([--with-lapack-lib=DIR],
                            [lapack lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure blas
AC_ARG_WITH([blas_lib],
            [AC_HELP_STRING([--with-blas-lib=DIR],
                            [blas lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure high precision libraries
AC_ARG_WITH([high_precision_include],
            [AC_HELP_STRING([--with-high-precision-include=DIR],
                            [mpc, mpfr and gmp include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([high_precision_lib],
            [AC_HELP_STRING([--with-high-precision-lib=DIR],
                            [mpc, mpfr and gmp lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Extra arguments to configure gsl
AC_ARG_WITH([gsl_include],
            [AC_HELP_STRING([--with-gsl-include=DIR],
                            [gsl include directory])],
            [CPPFLAGS="-I$withval $CPPFLAGS"])
            
AC_ARG_WITH([gsl_lib],
            [AC_HELP_STRING([--with-gsl-lib=DIR],
                            [gsl lib directory])],
            [LDFLAGS="-L$withval $LDFLAGS"])

# Option to use a speed optimized map
AC_ARG_ENABLE([speed_optimization],
  [AS_HELP_STRING([--enable-speed-optimization], [Set to yes to use a speed optimized matrix  construction method (requires google dense hash)])])
  
# Option to use a memory optimized map
AC_ARG_ENABLE([memory_optimization],
  [AS_HELP_STRING([--enable-memory-optimization], [Set to yes to use a memory optimized matrix  construction method (requires google sparse hash)])])

       
# Remove the default CXXFLAGS.
: ${CXXFLAGS=""}

# Add check for c++11 availability and set other compiler flags
CPPFLAGS="$CPPFLAGS -std=c++11 -O3  -Wall -Wfatal-errors"

# Check for a MPI c++ compiler, preferably one from intel
AC_PROG_CXX([mpicpc mpicxx mpicxxintel mpg++])

# Check for an intel fortran compiler, to be used for linking to intel mkl
# libraries only.

AS_IF([test "x$enable_intel_mkl" = "xyes"],[AC_PROG_FC([mpif90intel mpif90 mpifort])])

AS_IF([test "x$enable_intel_mkl" != "xyes"],[FC="$CXX"])

# Set the language to C++ in preparation for the further checks
AC_LANG(C++)

# Precompiler defs
AC_DEFINE(_DEBUG_,0)
AC_DEFINE(_ENABLE_MPI_,1)
AC_DEFINE(_ENABLE_HIGH_PRECISION_,1)
AC_DEFINE(_ENABLE_SPHERE_GEOMETRY_,1)
AC_DEFINE(_PYTHON_VERSION_,2)
AC_DEFINE(_ENABLE_DISC_GEOMETRY_,1)
AC_DEFINE(_ENABLE_TORUS_GEOMETRY_,1)
AC_DEFINE(_ENABLE_LATTICE_FQHE_,1)
AC_DEFINE(_WITH_GENERALISED_THETA_FUNCS_,1)
AC_DEFINE(_ENABLE_ANALYSIS_METHODS_,1)
AC_DEFINE(_CF_BENCHMARK_MODE_,0)

# Check for high precision library headers
AC_CHECK_HEADERS(gmp.h)
AC_CHECK_HEADERS(mpfr.h)
AC_CHECK_HEADERS(mpc.h)

# Check for fftw3 header
AC_CHECK_HEADERS(fftw3.h)

# Check for gsl_integration.h header
AC_CHECK_HEADERS(gsl/gsl_integration.h)

# Check for boost program options header
AC_CHECK_HEADERS(boost/program_options.hpp)

# Add fortran compiler libraries to the LIB variable
AS_IF([test "x$enable_intel_mkl" = "xyes"],[LIBS="$LIBS -mkl -nofor_main -cxxlib"])

AS_IF([test "x$enable_speed_optimization" = "xyes"],[AC_DEFINE(_SPEED_OPTIMIZED_MAP_,1) AC_DEFINE(_MEMORY_OPTIMIZED_MAP_,0) AC_CHECK_HEADERS(sparsehash/dense_hash_map)])

AS_IF([test "x$enable_memory_optimization" = "xyes"],[AC_DEFINE(_SPEED_OPTIMIZED_MAP_,0) AC_DEFINE(_MEMORY_OPTIMIZED_MAP_,1) AC_CHECK_HEADERS(sparsehash/sparse_hash_map)])

#   Add common library flags
LIBS="$LIBS -lblas -llapack -lm -lgsl -lfftw3 -lboost_program_options -lmpc -lmpfr -lgmp"

# Automatically generate Makefiles in each directory level
AC_OUTPUT(Makefile src/Makefile)
