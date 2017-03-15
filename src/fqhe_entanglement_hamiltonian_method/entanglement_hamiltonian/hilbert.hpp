////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!	 	Header file for the HilbertSpace class. Used to generate and diagonalize
//!	 	a matrix representation of an entanglement Hamiltonian written in terms
//!	 	of U(1) current operators.
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

#ifndef _HILBERT_HPP_INCLUDED_
#define _HILBERT_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "term.hpp"
#include "list_of_terms.hpp"
#include "../../occupation_basis/occupation_basis.hpp"
#include "../../utilities/wrappers/lapack_wrapper.hpp"
#include "../../utilities/algorithms/quick_sort.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

//!
//!	The HilbertSpace class contains data about the Hilbert space built from 
//!	U(1) current operators, it provides functions to populate a matrix representation
//!	of the Hamiltonian and to diagonalize it, writing the output data to files.
//!
class HilbertSpace
{
	private:
	std::vector<Term> m_hilbertSpace;     //!<   A vector of Hilbert space kets written
	                                      //!    in terms of U(1) current operators and
	                                      //!    stored as a vector of Term objects 		
	double* m_matrixRepresentation;       //!<   Double array to store the matrix
	                                      //!    representation of the Hamiltonian
	int m_dimension;                      //!<   Dimension of the Hilbert space
	                                      //!
	bool m_useBlocks;                     //!<   Option to use block optimizations of not
	                                      //!
	double** m_blockRepresentation;	      //!<	 A list of arrays to store the blocks
										  //!	 comprising the Hamiltonian matrix
	std::vector<int> m_blockDims;	      //!<	 A list of block dimensions for the 
										  //!    case where the matrix contains a block structure
	std::vector<int> m_blockLabels;       //!<   A list of labels to distinguish different 
	                                      //!    branch types 
	double* m_eigenvalues;                //!<   Array containing eigenvalues
	                                      //!
	int* m_eigenvalueLabels;              //!<   A list of eigenvalue labels
	                                      //!
	bool m_hilbertSpaceBuilt;             //!<   Flag to specify that the Hilbert 
	                                      //!<   space is built
	bool m_matrixGenerated;               //!<   Flag to specify that the matrix is generated
	                                      //!
	bool m_matrixDiagonalized;            //!<   Flag to specify that the matrix is diagonalized
    std::string GenFileName(boost::program_options::variables_map* optionList,
	                        const int lzA) const; 
	public:
	HilbertSpace();
	~HilbertSpace();
	void BuildHilbertSpace(boost::program_options::variables_map* optionList,
	                       const int lzA);
	void GenerateMatrixElements(ListOfTerms hamiltonian);
	void Diagonalize();
    void EigenvaluesToFile(boost::program_options::variables_map* optionList,
                           const int lzA) const;
};
#endif
