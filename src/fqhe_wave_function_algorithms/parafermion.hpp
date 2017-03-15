////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This is the header file for the Parafermion wave function class
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

#ifndef _PARAFERMION_HPP_INCLUDED_
#define _PARAFERMION_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <utilities/mathematics/binomial_table.hpp>
#include <utilities/mathematics/binary_number_tools.hpp>
#include <utilities/mathematics/dense_linear_algebra.hpp>
#include <vector> 
#include <sstream>
#if _DEBUG_
#include <utilities/general/debug.hpp>
#include <bitset>
#endif

namespace FQHE
{
    ////////////////////////////////////////////////////////////////////////////////
    //!	\brief	This class defines the basics data associated with, and functions 
    //!	to evaluate, parafermion type FQHE wave functions.
    //!	
    //!	For more information see e.g. 
    ////////////////////////////////////////////////////////////////////////////////
    class Parafermion : public WaveFunction
    {
        private:
        unsigned int m_nCluster;
        unsigned int m_nLastCluster;
        dcmplx* m_pfaff;
        std::vector<std::vector<unsigned char> > m_clusterPermutations;
        void BuildSubset(const unsigned int nRemain, std::vector<unsigned int>* nCluster, 
                         std::vector<unsigned char>* remainList, std::vector<unsigned char>* currPerm, 
                         std::vector<std::vector<unsigned char> >* clusterPermutations,
                         const unsigned int level) const;
        public:
	    Parafermion(WaveFunctionData*);	
	    ~Parafermion();
	    #if _ENABLE_SPHERE_GEOMETRY_
	    dcmplx EvaluateWfSphere(const int, dcmplx*, dcmplx*) const;
	    #endif
	    #if _ENABLE_DISC_GEOMETRY_	    
	    dcmplx EvaluateWfDisc(const int, dcmplx*) const;
	    #endif
    };
}   //  End namespace FQHE
#endif
