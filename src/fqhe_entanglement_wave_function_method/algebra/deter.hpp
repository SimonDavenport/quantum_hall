#ifndef INCLUDE_Deter
#define INCLUDE_Deter

#include <complex>
#include <iostream>

typedef std::complex<double> dcmplx;

extern "C"{
           int zgetrf_ ( int* , int* , dcmplx *, int* , int *, int * );
          }

dcmplx LogDeterminant(dcmplx** M, int dim) ;


#endif
