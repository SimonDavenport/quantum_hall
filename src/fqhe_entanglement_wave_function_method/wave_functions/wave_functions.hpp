#ifndef INCLUDE_Wave_function
#define INCLUDE_Wave_function

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include<math.h>
#include<memory.h>
#include<complex>



typedef std::complex<double> dcmplx;



class Wave_Functions
{
private:

/////////////// FQHE WAVEFUNCTIONS  ///////////////////////

void CF_EWF_Jain_Proj_3_LL (dcmplx &, double *, int, int, dcmplx *, dcmplx *, double **);
void CF_EWF_Jain_UnProj_3_LL (dcmplx &, double *, int, int, dcmplx *, dcmplx *, double **);//Jain book Proj//
void CF_EWF_Jain_Proj_2_LL (dcmplx &, double *, int, int, dcmplx *, dcmplx *, double **);//Jain book Proj//
void CF_EWF_Jain_UnProj_2_LL (dcmplx &, double *, int, int, dcmplx *, dcmplx *, double **);
void CF_EWF_Laugh (dcmplx &, double *, int, int, dcmplx *, dcmplx *);
void CF_EWF_Jain_negative_2_LL (dcmplx &, double *, int, int, dcmplx *, dcmplx *, double **);//Jain book Proj//


int *IPVT;
double *Lz_EWF_aux, **binom, *AT;
dcmplx **Matrix;

int ne_c;
int npart_c;
int LL_c;
int nu_c;
char flux;


/////////////// METROPOLIS WAVEFUNCTIONS /////////////////


public:

Wave_Functions(int, int, int, int, char);
~Wave_Functions();

void w_func (dcmplx &, double *, int , dcmplx *, dcmplx *, double **, int , int );   ///  evaluate w_func 
void Laugh_Test(dcmplx &, int, int, dcmplx *, dcmplx *, int, int, char);
void Laugh_Test_A(dcmplx &, int, int, dcmplx *, dcmplx *, int, int);
void Laugh_Test_B(dcmplx &, int, int, dcmplx *, dcmplx *, int, int);
void Laugh_Sphere(dcmplx &, int, int, dcmplx *, dcmplx *, int);

};

#endif
