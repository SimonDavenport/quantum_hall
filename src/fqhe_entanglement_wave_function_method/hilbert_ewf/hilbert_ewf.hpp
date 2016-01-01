#ifndef INCLUDE_Hilbert_EWF
#define INCLUDE_Hilbert_EWF

#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <memory.h>
#include <complex>

typedef std::complex<double> dcmplx;

void EWF_Space(double, int, int &, int, int, std::vector <std::vector <double> > & , std::vector<std::vector <double> > &, int);
void partition(int, int, int, std::vector< std::vector <double> > &, int &);
void bin_complement(std::vector< std::vector <int> >, std::vector< std::vector <int> > &, int, int, int);
void sign_binary(std::vector< std::vector <int> >, int *, int, int, int);
void inv_binary(std::vector< std::vector <double> > &, std::vector< std::vector <int> >, int, int, int);
void binary(std::vector< std::vector <double> >, std::vector< std::vector <int> > &, int, int, int);


#endif
