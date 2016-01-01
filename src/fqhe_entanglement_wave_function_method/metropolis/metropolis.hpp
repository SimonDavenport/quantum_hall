#ifndef INCLUDE_Metropolis
#define INCLUDE_Metropolis

#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <memory.h>
#include <complex>
#include "../../utilities/mathematics/mt.hpp"	
#include "../wave_functions/wave_functions.hpp"
#include "../hilbert_ewf/hilbert_ewf.hpp"	
#include "../algebra/algebra.hpp"
#include <complex>
#include <vector>
#include "../wave_functions/wave_functions.hpp"

//#include "../Algebra/determinant.hpp"
#include "../algebra/deter.hpp"
//#include "../Algebra/determinant1.hpp"
#include "../gnuplot/gnuplot_i.hpp"

typedef std::complex<double> dcmplx;
extern MersenneTwister mt;
extern double  **binom;


class metropolis
     {
      private:
      char AorB;
      double *Lz_EWF;
      dcmplx *zu, *zv, *zu_move, *zv_move, *zz_move, **Matrix1;
      Wave_Functions *wf;  

      public:
      dcmplx *zz;
      int np, np_A, iter, LL, nu, nu_eff, accepted, rejected;
      double cut_sph, width;



      metropolis(dcmplx *zz1, int np1, int iter1, double width1, int np_A1, int nu1, int LL1, int nu_eff1, 
                 double cut_sph1);
      metropolis(dcmplx *zz1, int np1, int iter1, int width1);

      ~metropolis() ;



      void run_A();
      void run_B();
      void run_Sphere_Laughlin(int nu1);
      void run_Sphere();

      dcmplx IQHE_Prob_A(dcmplx *zu1, dcmplx *zv1);
      dcmplx IQHE_Prob_B(dcmplx *zu1, dcmplx *zv1);
      dcmplx Laugh_Sphere(dcmplx *zu1, dcmplx *zv1, int nu1);
      dcmplx Jain_negative_2_LL_A (dcmplx *zu1, dcmplx *zv1);
      dcmplx Jain_negative_2_LL_B (dcmplx *zu1, dcmplx *zv1);
      dcmplx EWF_Expansion_State(dcmplx *zu1, dcmplx *zv1, int npart_A);



      std::vector<double> correl(std::vector <double> Obs, int num_iter, double mean, int range);
      void plot(std::vector<double> Corr, int range);
      dcmplx Wfunction_Norm (double width_norm, int iter_norm, int init_norm, int metrop_norm, dcmplx &alpha);

      void Wfunction_Corr_IQHE_N(int num_metrop_corr, int range);
      void Wfunction_Corr_IQHE_NA(int num_metrop_corr, int range);
      void Wfunction_Corr_IQHE_NB(int num_metrop_corr, int range);
      void Wfunction_Corr_N(int num_metrop_corr, int range);





      void Wfunction_Norm_Thermal (double width_norm, int metrop_norm, int range);      
      void Coulomb_Thermal (double width_norm, int metrop_norm, int range);      
      dcmplx Wfunction_Coulomb_energy (int initial, int num_metrop_iter);

      void equilibration(int accep, int rejec);



      };








/*

void metropolis(int, dcmplx *, double *, double, int, int, int, char, double);
void metropolis_A(int, dcmplx *, double, int, int, int, int, int, double, int &, int &);
void metropolis_B(int, dcmplx *, double, int, int, int, int, int, double, int &, int &);
//void metropolis_AB(int, dcmplx *, double, int, double **, int, int, int, int, double);
void Laugh_Test_B(dcmplx &, int, int, dcmplx *, dcmplx *, int, int, double);
void Laugh_Test_A(dcmplx &, int, int, dcmplx *, dcmplx *, int, int, double);
void metropolis_Sphere(int, dcmplx *, double, int, int, int &, int &);
void Laugh_Sphere(dcmplx &, int, int, dcmplx *, dcmplx *, int);

void CF_EWF_Jain_negative_2_LL_A (dcmplx &, int , int , dcmplx *, dcmplx *);
void CF_EWF_Jain_negative_2_LL_B (dcmplx &, int , int , dcmplx *, dcmplx *);
*/

#endif
