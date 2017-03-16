////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This is the header file for the Laughlin wave function class
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

#ifndef _LAUGHLIN_HPP_INCLUDED_
#define _LAUGHLIN_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <iostream>
#include <time.h>
#if _ENABLE_TORUS_GEOMETRY_	
#include <utilities/mathematics/theta_functions.hpp>
#endif

namespace FQHE
{
    ///////     STATIC CONST DEFINITIONS      //////////////////////////////////////
    #if _ENABLE_TORUS_GEOMETRY_	
    static const double laughlinComZero = sqrt(2.0);
    #ifndef _TORUS_SAMETOL_DEFINED_
    #define _TORUS_SAMETOL_DEFINED_
    static const double sameTol = pow(10.0, -14.0);
    #endif
    #endif
 
    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief	This class defines the basics data associated with, and functions 
    //!	to evaluate, Laughlin type FQHE wave functions.
    //!	
    //!	For more information see e.g. PRL 50 1395--1398 (1983).
    ////////////////////////////////////////////////////////////////////////////////
    class Laughlin : public WaveFunction
    {
        #if _ENABLE_TORUS_GEOMETRY_
        private:
        //  The following variables are used for theta function pre-calculations
        //	for quasiholes on a lattice
        utilities::ThetaLookUp *m_thetaFuncs;//!< 	To store a table of theta function values
									    //!	  	on the lattice with no offset
        utilities::ThetaLookUp **m_thetaFuncsQh;//!<	Store a list of tables of theta 
									    //!		function values offset by each 
									    //!  	unique quasi hole position
        dcmplx *m_uniqueQhPositions;	//!<	A list of unique quasi hole coordinates
									    //!
        int m_nbrQhPositions;	        //!<	The number of unique quasi hole coordinates
									    //!		in the list
        utilities::ThetaLookUp **m_thetaFuncsCom;//!<	Store a list of tables of theta function
									    //!		values offset by the centre of mass coordinate
									    //!	 	(minus the COM zero depending on whether we use
									    //!		the COM zeros formalism or not)
        dcmplx *m_uniqueCom;			//!<	A list of unique COM (centre-of-mass) coordinates
									    //!
        int m_nbrCom;			        //!<	The number of unique COM coordinates
									    //!		in the list
	    #if _BENCHMARK_MODE_
        double m_timeJastrow;	        //!<	Total time taken to evaluate Jastrow factors
							            //!
        double m_timePfaffian;	        //!<	Total time taken to evaluate Pfaffian terms
							            //!
        double m_timeGauss;		        //!<	Total time taken to evaluate Gaussian factor
							            //!
        double m_timeQuasihole;	        //!<	Total time taken to evaluate quasi hole part
							            //!
        double m_timeCM;		        //!<	Total time taken to evaluate COM terms
							            //!
        double m_timer;			        //!<	Timer variable
							            //!
	    #endif
	    #endif
        public:

	    Laughlin(WaveFunctionData*);
	    ~Laughlin();
	    #if _ENABLE_SPHERE_GEOMETRY_
	    dcmplx EvaluateWfSphere(const int, dcmplx*, dcmplx*) const;
	    dcmplx EvaluateQuasiholeWfSphere(const int, dcmplx*, dcmplx*, const int, dcmplx*, dcmplx*) const;     
	    #endif
	    #if _ENABLE_DISC_GEOMETRY_		    
	    dcmplx EvaluateWfDisc(const int, dcmplx*) const;
	    #endif
	    #if _ENABLE_TORUS_GEOMETRY_
	    void LatticeDataPrecalculation(const int, const int, dcmplx*);
	    dcmplx EvaluateWfTorus(const int, std::complex<int>*, const int) const;
        dcmplx EvaluateQuasiholeWfTorus(const int, std::complex<int>*, const int, const int, dcmplx*) const;
        #endif
    };
}   //  End namespace FQHE
#endif
