////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		This is the header file for the thetaFunctions.cpp library
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

#ifndef _THETA_FUNCTIONS_HPP_INCLUDED_
#define _THETA_FUNCTIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////       
#include "../general/dcmplx_type_def.hpp"
#include "../general/pi_const_def.hpp"
#include "../general/i_const_def.hpp"
#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#if _DEBUG_
#include "../general/debug.hpp"
#endif

///////		STATIC CONSTANT DECLARATIONS		    ////////////////////////////

//  Tolerance level required to stop calculation and return a result;
static const double thetaTol = pow(10.0, -15.0);

namespace utilities
{
    //!
    //!	A function Namespace for theta function routines.
    //!
    namespace thetaFunction
    {
        //	Generalised Jacobi theta function
        dcmplx   GeneralisedJacobi(const double, const double, dcmplx, const double);
    }

    //!
    //!	To contain theta look-up table data.
    //!
    class ThetaLookUp
    {
	    public:
	    dcmplx *thetaLookUpTable;	//!<	Look-up table for theta function values on
								    //!		a lattice of size Lx by Ly
	    ThetaLookUp(const int Lx, const int Ly, const double arg1, const double arg2,
				    const dcmplx quasiHoleOffset, const bool restrictSign, const int increaseDomain);
	    ~ThetaLookUp();
    };
}   //  End namespace utilities
#endif
