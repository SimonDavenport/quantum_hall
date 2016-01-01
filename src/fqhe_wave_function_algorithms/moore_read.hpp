////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 04/12/2013
//!
//!  \file 
//!		This is the header file for the Moore-Read class
//!
//!                    Copyright (C) 2013 Simon C Davenport
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

#ifndef _MOORE_READ_HPP_INCLUDED_
#define _MOORE_READ_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <utilities/mathematics/dense_linear_algebra.hpp>

#if _ENABLE_TORUS_GEOMETRY_
	#include <utilities/mathematics/theta_functions.hpp>
#endif

//////////////////////////////////////////////////////////////////////////////
//!	\brief	A namespace containing all functions and classes used in 
//!	fractional quantum Hall effect calculations. 
//////////////////////////////////////////////////////////////////////////////
namespace FQHE
{

///////     STATIC CONST DEFINITIONS      //////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_TORUS_GEOMETRY_	

//	Define a small shift factor to deal with theta function zeros
static const double smallShift = 0.001;

//	Define COM zero(s)
static const double mooreReadComZero=sqrt(2.0);

#ifndef _TORUS_SAMETOL_DEFINED_
#define _TORUS_SAMETOL_DEFINED_

static const double sameTol = pow(10.0,-14.0);

#endif
//	This constant determines the minimum difference to define continuous 
//	quasi-hole positions to be coincident

#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!	\brief	This class defines the basics data associated with, and functions 
//!	to evaluate, Moore--Read or Pfaffian type FQHE wave functions.
//!	
//!	For more information see e.g. Nuclear Physics B360 (1991) 362-396 and do 
//!	a forward citation search!
////////////////////////////////////////////////////////////////////////////////

class MooreRead : public WaveFunction
{
	private:
	
	dcmplx *m_pfaff;		//!<	An array to store the Pfaffian matrix
							//!
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_
	
	//  The following variables are used for theta function pre-calculations 
    //	for quasi holes on a lattice

    utilities::ThetaLookUp **m_thetaFuncs;	//!< 	To store a table of theta function values
										//!	  	on the lattice with no offset
    utilities::ThetaLookUp **m_thetaFuncsQh;//!<	Store a list of tables of theta 
										//!		function values offset by each 
										//!  	unique quasi hole position
    dcmplx *m_uniqueQhPositions;		//!<	A list of unique quasi hole coordinates
										//!
    int m_nbrQhPositions;		        //!<	The number of unique quasi hole coordinates
										//!		in the list
    utilities::ThetaLookUp **m_thetaFuncsW1W2;//!<	Store a list of tables of theta function
										//!		values offset by values of (W1-W2)/2
    dcmplx *m_uniqueW1W2;				//!<	A list of unique values of (W1-W2)/2	
										//!
    int m_nbrW1W2;				        //!<	The number of unique (W1-W2)/2 values 
										//!		in the list
    utilities::ThetaLookUp **m_thetaFuncsW1W2W3W4A;//!<	Store a list of tables of theta function
										//!		values offset by values of (W1+W2-W3-W4)/2
    dcmplx *m_uniqueW1W2W3W4A;			//!<	A list of unique values of (W1+W2-W3-W4)/2
										//!
    int m_nbrW1W2W3W4A;		            //!<	The number of unique (W1+W2-W3-W4)/2 values 
										//!		in the list
    utilities::ThetaLookUp **m_thetaFuncsW1W2W3W4B;//!<	Store a list of tables of theta function
										//!		values offset by values of (W1-W2+W3-W4)/2
    dcmplx *m_uniqueW1W2W3W4B; 			//!<	A list of unique values of (W1-W2+W3-W4)/2
										//!
    int m_nbrW1W2W3W4B;		            //!<	The number of unique (W1-W2+W3-W4)/2 values 	
										//!		in the list
    utilities::ThetaLookUp **m_thetaFuncsCom;//!<	Store a list of tables of theta function
										//!		values offset by the centre of mass coordinate
										//!	 	(minus the COM zero depending on whether we use
										//!		the COM zeros formalism or not)
    dcmplx *m_uniqueCom;				//!<	A list of unique COM (centre-of-mass) coordinates
										//!
    int m_nbrCom;				        //!<	The number of unique COM coordinates
										//!		in the list
				
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_
	
    double m_timeJastrow;	//!<	Total time taken to evaluate Jastrow factors
							//!
    double m_timePfaffian;	//!<	Total time taken to evaluate Pfaffian terms
							//!
    double m_timeGauss;		//!<	Total time taken to evaluate Gaussian factor
							//!
    double m_timeQuasihole;	//!<	Total time taken to evaluate quasi hole part
							//!
    double m_timeCM;		//!<	Total time taken to evaluate COM terms
							//!
    double m_timer;			//!<	Timer variable
							//!
	#endif
	//////////////////////////////////////////////////////////////////////////////
	
	//	Function declarations
	
	public:
	
	MooreRead(WaveFunctionData*);	//	constructor
	~MooreRead();					//	destructor
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_SPHERE_GEOMETRY_
	
	//	Evaluate the Moore-Read wave function in the sphere geometry
	dcmplx EvaluateWfSphere(const int,dcmplx*,dcmplx*) const;	
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_DISC_GEOMETRY_
	
	//	Evaluate the Moore-Read wave function in the disc geometry		
	dcmplx EvaluateWfDisc(const int,dcmplx*) const;		
		
	#endif
	//////////////////////////////////////////////////////////////////////////////////	
		
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_	
		
	void PrintPfaffian(const int dim) const;	//	Print out the Pfaffian matrix elements

	#endif
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_

	//	Used for pre-calculations of theta function values
	void LatticeDataPrecalculation(const int,const int,dcmplx*);
	
	//	Evaluate the Moore-Read wave function in the torus geometry
	dcmplx EvaluateWfTorus(const int,std::complex<int>*,const int) const;			
					
    dcmplx EvaluateQuasiholeWfTorus(const int,std::complex<int>*,const int,const int,dcmplx*) const;
    //	Evaluate Moore-Read wave function with quasi-particles in sphere geometry

	#endif
	//////////////////////////////////////////////////////////////////////////////////

};

}

#endif
