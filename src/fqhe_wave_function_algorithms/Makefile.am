otherincludedir = $(includedir)/fqhe_wave_functions

nobase_otherinclude_HEADERS = *.hpp

lib_LIBRARIES = libfqhe_wave_functions.a libfqhe_wave_functions_MPI.a

libfqhe_wave_functions_a_SOURCES = \
fqhe_wave_function.cpp \
composite_fermion.cpp \
nass.cpp \
moore_read.cpp \
laughlin.cpp \
parafermion.cpp

libfqhe_wave_functions_MPI_a_SOURCES = \
fqhe_wave_function.cpp \
composite_fermion.cpp \
nass.cpp \
moore_read.cpp \
laughlin.cpp \
parafermion.cpp

libfqhe_wave_functions_a_CPPFLAGS = -D_ENABLE_MPI_=0
libfqhe_wave_functions_MPI_a_CPPFLAGS = -D_ENABLE_MPI_=1
