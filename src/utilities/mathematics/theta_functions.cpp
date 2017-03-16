////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		 This is the library file containing functions to evaluate
//!      Jacobi theta functions. The method is based simply on doing a
//!      power series expansion of the sum and then truncating when the
//!      errors fall to within some tolerance limit
//!      see http://mathworld.wolfram.com/JacobiThetaFunctions.html for 
//!      more details about theta functions 
//!
//!                    Copyright (C) Simon C Davenport
//!
//!		This program is free software: you can redistribute it and/or modify
//!		it under the terms of the GNU General Public License as published by
//!		the Free Software Foundation, either version 3 of the License,
//!		or (at your option) any later version.
//!
//!		This program is distributed in the hope that it will be useful, but
//!		WITHOUT ANY WARRANTY; without even the implied warranty of 
//!		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
//!		General Public License for more details.
//!
//!		You should have received a copy of the GNU General Public License
//!		along with this program. If not, see <http://www.gnu.org/licenses/>.
//! 
////////////////////////////////////////////////////////////////////////////////

#include "theta_functions.hpp"

//////////////////////////////////////////////////////////////////////////////
//!	\brief  This is the constructor for an object of type ThetaLookUp.
//!
//!  The object contains a lookup table of specified generalized theta function
//!  values in order to avoid repeated calculations.
//!  The table is for theta function values in a lattice of size Lx by Ly 
//!  For complex integer arguments Lx + I Ly. There is an option to include
//!  a complex double values offset to all values.
//////////////////////////////////////////////////////////////////////////////
utilities::ThetaLookUp::ThetaLookUp(
	const int Lx,					//!<	 X dimension of the lattice for which
									//!		 to generate A table of values Theta (Lx + I Ly) 
	const int Ly,					//!<	 Y dimension of the lattice for which
									//!		 to generate A table of values Theta (Lx + I Ly) 
	const double arg1,				//!<	 First argument of generalized theta function
	const double arg2,				//!<	 Second argument of generalized theta function
	const dcmplx quasiHoleOffset,	//!<	 A complex double offset to all the Theta
									//!		 function arguments in the table i.e we tabulate
									//!		 Theta_[arg1,arg2] (x + I y -quasiHoleOffset)
	const bool restrictSign,		//!<	 Restrict the sign of the arguments to only be 
									//!		 positive (then map other quadrants onto that)
	const int increaseDomain)		//!<	 Extend the range of values stored by this factor
									//!		 ALSO increase tau by this factor (for special
									//!		 application to torus quantum Hall wave functions)
{
    double tau = increaseDomain*(double)Ly/Lx;
    if(restrictSign)        //  only generate the positive quadrant (i>0 and j>0)
	{
	    thetaLookUpTable = new dcmplx[Lx*Ly*increaseDomain*increaseDomain];
	    for(int i=0; i<Lx*increaseDomain; ++i)
	    {
		    for(int j=0; j<Ly*increaseDomain; ++j)
		    {
			    thetaLookUpTable[i*Ly*increaseDomain+j] = thetaFunction::GeneralisedJacobi(
			        arg1, arg2, (dcmplx(i, j)-quasiHoleOffset)/(double)Lx, tau);
		    }
	    }
    }
    else    //  store all possible values (takes 4x the memory)
    {
	    thetaLookUpTable=new dcmplx[4*Lx*Ly*increaseDomain*increaseDomain];
	    for(int i=0; i<2*Lx*increaseDomain; ++i)
	    {
		    for(int j=0; j<2*Ly*increaseDomain; ++j)
		    {
			    thetaLookUpTable[i*(2*Ly*increaseDomain)+j] = thetaFunction::GeneralisedJacobi(
			        arg1, arg2, (dcmplx(i-Lx, j-Ly)-quasiHoleOffset)/(double)Lx, tau);
		    }
	    }
    }
}

//!
//!	Destructor for an object of type ThetaLookUp.
//!
utilities::ThetaLookUp::~ThetaLookUp()
{
	delete[] thetaLookUpTable;
}

//////////////////////////////////////////////////////////////////////////////
//!	\brief  The generalised Jacobi theta function is given by 
//!	\f$ \sum_{-\infty}^{\infty} exp(-\pi*t*(n+a)^2+2*I*\pi(n+a)(z+b)) \f$.
//!	
//!  The following a and b parameters in the generalized theta function 
//!  correspond to the standard definitions of Theta_1,...,Theta_4:
//!	
//!	\verbatim
//!			            a        b  
//!	    Theta_1 	   0.5		-0.5
//!	    Theta_2 	   0.5		0.0 
//!	    Theta_3 	   0.0		0.0 
//!	    Theta_4 	   0.0		0.5  
//!	
//!	\endverbatim
//!
//!  Note: the generalized theta functions are quite badly behaved.
//!  So it's important to check the convergence of the summations carefully
//!
//!	\return   complex LOG of theta function	
//////////////////////////////////////////////////////////////////////////////
dcmplx utilities::thetaFunction::GeneralisedJacobi(
	const double a,		//!<	First parameter in generalized theta function 
	const double b,		//!<	Second parameter in generalized theta function
	dcmplx z,			//!<	Complex argument used as written in the description
	const double t)		//!<	The 'nome' of the theta function, used as
						//!		written in the description. N.B. sometimes defined as exp(q) 
{
    int n;
    int nbrTerms;
    dcmplx thetaFunc,prevTheta,nextTerm;
    dcmplx xShiftFactor=0;
    dcmplx yShiftFactor=0;
    int extraTerms = 5+ceil(a)+ceil(b);			//	set the minimum number of terms
    //	Shift z to be as close as possible to zero and add in the periodic factors
    //	to compensate.
    while(real(z)>=1.0)
    {
    	z -= 1.0;
    	xShiftFactor += 2*PI*I*a;
    }
    while(real(z)<=-1.0)
	{
		z += 1.0;
		xShiftFactor -= 2*PI*I*a;
	}
    while(imag(z)>=t)
	{
		z -= I*t;
		yShiftFactor += (PI*t-2*PI*I*(z+b));
	}
    while(imag(z)<=-t)
	{
		z += I*t;
		yShiftFactor -= (PI*t-2*PI*I*(z+b));
	}
    //  calculate first term in the sum (n=0 term)
    thetaFunc = exp(PI*a*(-a*t+2.0*I*(z+b)));
    prevTheta = thetaFunc;
    //  calculate the next pair of terms (n=+1 and n=-1)
    nextTerm = exp(PI*(a+1.0)*(-(a+1.0)*t+2.0*I*(z+b)))+exp(PI*(a-1.0)*(-(a-1.0)*t+2.0*I*(z+b)));
    thetaFunc += nextTerm;
    n = 2;
    nbrTerms = 2;
    //	Keep going until until adding terms does not change the result
    while(abs(nextTerm/prevTheta)>thetaTol || extraTerms>0)
    {
        prevTheta += nextTerm;
        //  add next pair of terms to the sum
        nextTerm = exp(PI*(a+n)*(-(a+n)*t+2.0*I*(z+b)));
        nextTerm += exp(PI*(a-n)*(-(a-n)*t+2.0*I*(z+b)));
        thetaFunc += nextTerm;
        ++n;
        ++nbrTerms;
        --extraTerms;
    }
	thetaFunc = log(thetaFunc)+xShiftFactor+yShiftFactor;
	//	Use high precision if the low precision range overflows
    if(std::isinf(real(thetaFunc)) || std::isnan(real(thetaFunc)))
    {
		//	if the double range is overflowed, then we need to use high precision variables
	    //  Declare local variables
	    int prec = 256;
	    mpc_t runTot;
	    mpc_t prev,next1,next2;
	    mpc_t compare;
	    mpfr_t absCompare;
		mpc_init2(runTot, prec);
		mpc_init2(prev, prec);
		mpc_init2(next1, prec);
		mpc_init2(next2, prec);
		mpc_init2(compare, prec);
		mpfr_init(absCompare);
		//  calculate first term in the sum (n=0 term)
		nextTerm = PI*a*(-a*t+2.0*I*(z+b));
		mpc_set_d_d(runTot, real(nextTerm), imag(nextTerm), MPC_RNDNN);
		mpc_exp(runTot, runTot, MPC_RNDNN);
		mpc_set(prev, runTot, MPC_RNDNN);
		//  calculate the next pair of terms (n=+1 and n=-1)
		nextTerm = PI*(a+1.0)*(-(a+1.0)*t+2.0*I*(z+b));
		mpc_set_d_d(next1, real(nextTerm), imag(nextTerm), MPC_RNDNN);
		mpc_exp(next1, next1, MPC_RNDNN);
		mpc_add(runTot, runTot, next1, MPC_RNDNN);
		nextTerm = PI*(a-1.0)*(-(a-1.0)*t+2.0*I*(z+b));
		mpc_set_d_d(next2, real(nextTerm), imag(nextTerm), MPC_RNDNN);
		mpc_exp(next2, next2, MPC_RNDNN);
		mpc_add(runTot, runTot, next2, MPC_RNDNN);
		n = 2;
		nbrTerms = 2;
		mpc_sub(compare, runTot, prev, MPC_RNDNN);
		mpc_div(compare, compare, runTot, MPC_RNDNN);
		mpc_abs(absCompare, compare, MPFR_RNDN);
		while(mpfr_get_d(absCompare, MPFR_RNDN)>thetaTol)
		{
			mpc_add(prev, prev, next1, MPC_RNDNN);
			mpc_add(prev, prev, next2, MPC_RNDNN);
			 //  add next pair of terms to the sum
			nextTerm = PI*(a+n)*(-(a+n)*t+2.0*I*(z+b));
			mpc_set_d_d(next1, real(nextTerm), imag(nextTerm), MPC_RNDNN);
			mpc_exp(next1, next1, MPC_RNDNN);
			mpc_add(runTot, runTot, next1, MPC_RNDNN);
			nextTerm = PI*(a-n)*(-(a-n)*t+2.0*I*(z+b));
			mpc_set_d_d(next2, real(nextTerm), imag(nextTerm), MPC_RNDNN);
			mpc_exp(next2, next2, MPC_RNDNN);
			mpc_add(runTot, runTot, next2, MPC_RNDNN);
			mpc_sub(compare, runTot, prev, MPC_RNDNN);
			mpc_div(compare, compare, runTot, MPC_RNDNN);
			mpc_abs(absCompare, compare, MPFR_RNDN);
			++n;
			++nbrTerms;
		}
		mpc_log(runTot, runTot, MPC_RNDNN);
		mpc_real(absCompare, runTot, MPFR_RNDN);
		thetaFunc = mpfr_get_d(absCompare, MPFR_RNDN);
		mpc_imag(absCompare, runTot, MPFR_RNDN);
		thetaFunc += dcmplx(0, mpfr_get_d(absCompare, MPFR_RNDN));
		thetaFunc += xShiftFactor+yShiftFactor;
		mpc_clear(runTot);
		mpc_clear(prev);
		mpc_clear(next1);
		mpc_clear(next2);
		mpc_clear(compare);
		mpfr_clear(absCompare);
    }
    return thetaFunc;
}
