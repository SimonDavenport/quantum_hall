#ifndef INCLUDE_Algebra
#define INCLUDE_Algebra

#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <memory.h>
#include <complex>
#include "../../utilities/mathematics/mt.hpp"
#include "../wave_functions/wave_functions.hpp"

typedef std::complex<double> dcmplx;
extern MersenneTwister mt;

//void Det(int, dcmplx **, dcmplx*);
void binomial (double **, double);
void rand_init(char * argv[1]);
double nfmod(double, double);
double factorial(unsigned int n);

double GAMMLN(double XX);
double BETACF(double A,double B,double X);
double BETAI (double A,double B,double X); 
double Gamma(double xx);
double Beta_I(double X, double a, double b);
void combinationUtil(int arr[], int data[], int start, int end, int index, int r, std::vector <int> & aux);
void Combinations(int n, int r, int** comb, int**comb_comp, int data[]);


#endif
