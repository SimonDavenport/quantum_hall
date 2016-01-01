/********************************************
* Program to demonstrate the Gamma Function *
*                                           *
*               C++ version by J-P Moreau   *
*                   (www.jpmoreau.fr)       *
* ----------------------------------------- *
* Reference:                                * 
* "Numerical Recipes, By W.H. Press, B.P.   *
*  Flannery, S.A. Teukolsky and T. Vetter-  *
*  ling, Cambridge University Press, 1986"  *
*  [BIBLI 08].                              *
* ----------------------------------------- *
*                                           *
********************************************/
#include "algebra.hpp"

#define  half 0.5
#define  one  1.0
#define  fpf  5.5
#define  zero 0.0 

double x, y;
int i;  

double GAMMLN(double XX) {
// returns the value Ln(Gamma(XX) for XX>0
double STP,HALF,ONE,FPF,X,TMP,SER;
double COF[7];
int    J;
  COF[1]= 76.18009173; COF[2]=-86.50532033; COF[3]=24.01409822;
  COF[4]=-1.231739516; COF[5]=0.120858003E-2; COF[6]=-0.536382E-5;
  STP=2.50662827465;
  HALF=0.5; ONE=1.0; FPF=5.5;
  X=XX-ONE;
  TMP=X+FPF;
  TMP=(X+HALF)*log(TMP)-TMP;
  SER=ONE;
  for (J=1; J<=6; J++) {
    X=X+ONE;
    SER += COF[J]/X;
  }
  return (TMP+log(STP*SER));
}


/******************************************
*           FUNCTION  GAMMA(X)            *
* --------------------------------------- *
* Returns the value of Gamma(x) in double *
* precision as EXP(LN(GAMMA(X))) for X>0. *
******************************************/
double Gamma(double xx)  {
  double cof[7],stp,x,tmp,ser;
  int j;
  cof[1]=76.18009173;
  cof[2]=-86.50532033;
  cof[3]=24.01409822;
  cof[4]=-1.231739516;
  cof[5]=0.120858003e-2;
  cof[6]=-0.536382e-5;
  stp=2.50662827465;
  
  x=xx-one;
  tmp=x+fpf;
  tmp=(x+half)*log(tmp)-tmp;
  ser=one;
  for (j=1; j<7; j++) {
    x=x+one;
    ser=ser+cof[j]/x;
  }
  return (exp(tmp+log(stp*ser)));
}

// end of file gamma.cpp

