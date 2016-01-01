////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 04/12/2013
//!
//!  \file 
//!		This is the header file for the NASS class
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

#ifndef _NASS_HPP_INCLUDED_
#define _NASS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <sstream>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////
//!	\brief	A namespace containing all functions and classes used in 
//!	fractional quantum Hall effect calculations. 
//////////////////////////////////////////////////////////////////////////////
namespace FQHE
{

////////////////////////////////////////////////////////////////////////////////
//!	\brief	This class defines the basics data associated with, and functions 
//!	to evaluate, non-abelian spin singlet (NASS) type FQHE wave functions.
//!	
//!	See e.g. Nucl. Phys. B 607, 549--576 for more details about the NASS state
////////////////////////////////////////////////////////////////////////////////

class NonAbelianSpinSinglet : public WaveFunction
{
	private:
	//	Declare variables used in wave function evaluation algorithm

	int m_nbrPerms;		//!<	Number of unique co-ordinate permutations which
						//!		the NASS wf is not invariant under
	int *m_permsList;	//!<	An array containing a list of the indices of 
						//!		permuted co-ordinates
	int m_nOverTwo;		//!<	Dimensions of each row of array of permutations
						//!		i.e Half the number of particles
	int *m_sizeListA;	//!<	A list of the number of the number of 
						//!		transpositions required to get from one permutation
						//!		in the list to the following one (group A)
	int *m_sizeListB;	//!<	A list of the number of the number of 
						//!		transpositions required to get from one permutation
						//!		in the list to the following one (group B)
	int	*m_indexA1;		//!<	Array to store the indices of those transpositions (group A1)
						//!
	int	*m_indexA2;		//!<	Array to store the indices of those transpositions (group A2)
						//!
	int	*m_indexB1;		//!<	Array to store the indices of those transpositions (group B1)
						//!
	int	*m_indexB2;		//!<	Array to store the indices of those transpositions (group B2)
						//!
	dcmplx *m_zwDiff;	//!<	An array of terms z[i]-w[j]	
						//!
	dcmplx *m_cmplxCtr1;//!<	An array for complex number accumulation for group 1
						//!
	dcmplx *m_cmplxCtr2;//!<	An array for complex number accumulation for group 2	
						//!

	public:
	
	NonAbelianSpinSinglet(WaveFunctionData*);	//	constructor	
	~NonAbelianSpinSinglet();					//	destructor
	
	//////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_DISC_GEOMETRY_
	
	//	Evaluate NASSwavefunction in the disc geometry
	
	dcmplx EvaluateWfDisc(const int,dcmplx*) const;	
	
	#endif
	//////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_SPHERE_GEOMETRY_
	
	//	Evaluate NASS wave function in the sphere geometry
	
	dcmplx EvaluateWfSphere(const int,dcmplx*,dcmplx*) const;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////
	
};

}

#endif
