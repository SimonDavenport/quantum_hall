////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This is the header file for the Moore-Read class
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

#ifndef _MOORE_READ_HPP_INCLUDED_
#define _MOORE_READ_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <utilities/mathematics/dense_linear_algebra.hpp>
#if _ENABLE_TORUS_GEOMETRY_
	#include <utilities/mathematics/theta_functions.hpp>
#endif
namespace FQHE
{
    ///////     STATIC CONST DEFINITIONS      //////////////////////////////////////
    #if _ENABLE_TORUS_GEOMETRY_	
    static const double smallShift = 0.001;
    static const double mooreReadComZero = sqrt(2.0);
    #ifndef _TORUS_SAMETOL_DEFINED_
    #define _TORUS_SAMETOL_DEFINED_
    static const double sameTol = pow(10.0,-14.0);
    #endif
    #endif
    
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
	    #if _ENABLE_TORUS_GEOMETRY_
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

	    public:
	    MooreRead(WaveFunctionData*);
	    ~MooreRead();
	    #if _ENABLE_SPHERE_GEOMETRY_
	    dcmplx EvaluateWfSphere(const int, dcmplx*, dcmplx*) const;	
	    #endif
	    #if _ENABLE_DISC_GEOMETRY_		
	    dcmplx EvaluateWfDisc(const int, dcmplx*) const;		
	    #endif
	    #if _TEST_MODE_	
	    void PrintPfaffian(const int dim) const;
	    #endif
	    #if _ENABLE_TORUS_GEOMETRY_
	    void LatticeDataPrecalculation(const int, const int, dcmplx*);
	    dcmplx EvaluateWfTorus(const int, std::complex<int>*, const int) const;				
        dcmplx EvaluateQuasiholeWfTorus(const int, std::complex<int>*, const int, const int, dcmplx*) const;
	    #endif
    };
}   //  End namespace FQHE
#endif
