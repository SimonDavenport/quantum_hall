////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 03/11/2014
//!
//!  \file 
//!   	This file contains functions for generating the occupation basis of 
//!     Fractional Quantum Hall Quasihole wave functions for the construction
//!     of the associated real space entanglement spectrum.
//!
//!                    Copyright (C) 2014 Simon C Davenport
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

#ifndef _OCCUPATION_BASIS_HPP_INCLUDED_
#define _OCCUPATION_BASIS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include <math.h>   //  For fabs
#include <vector>   //  For std::vector
#include "../utilities/general/cout_tools.hpp"
#include "../utilities/mathematics/integer_partitions.hpp"
#include "../utilities/mathematics/combinatorics.hpp"
#include <bitset>   //  For bitwise representation of Landau level occupations

#if _DEBUG_
#include "../utilities/debug.hpp"
#endif

///////     STATIC CONSTANTS     ///////////////////////////////////////////////

//  Define the maximum number of particles that can be represented
static const int _MAX_NUMBER_ = 64;

///////     FUNCTION DECLARATIONS     //////////////////////////////////////////

namespace occupationBasis
{

int FindNearestValue(const int b,const int a1,const int a2);

double FindLzA(std::vector<int>& occupations,int nbrLowest);

void GenerateOccupationData(const int nbrParticles,const int nbrA,const int lzA, 
	const int nbrColours,std::vector<std::vector<int> >& deltaNa0List,
    std::vector<std::vector<int> >& nbrEachLL,std::vector<int>& deltaLza,
    std::vector<std::vector<int> >& maxExcitationSize);

void GenerateLevelOccupations(std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,
    std::vector<int>& spectrumLabels,const int nbr,const int nbrA,const int nbrLevels,const int sector);

void AssignEnergies(std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,double* energyLevels,
    double* spectrum,const double normalization,const int nbrStates,const int nbrLevels);
    
void AssignSquaredEnergies(std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,double* energyLevels,
    double* spectrum,const double normalization,const int nbrStates,const int nbrLevels);
    
};  //  End occuaptionBasis namespace
    
#endif
