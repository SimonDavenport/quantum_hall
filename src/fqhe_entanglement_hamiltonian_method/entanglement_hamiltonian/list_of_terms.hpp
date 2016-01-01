////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 10/04/2014
//!
//!  \file 
//!		The ListOfTerms class contains a list of
//!     Term objects and functions to perform commutations.
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

#ifndef _LIST_OF_TERMS_HPP_INCLUDED_
#define _LIST_OF_TERMS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "term.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include <vector>
#include <iostream>
#include <iomanip>

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief The ListOfTerms class contains functions to perform recursive commutation
//!	operations on matrix elements of Hamiltonians constructed form U(1)
//!	current operators and Majora operators. 
//! 
//!	The Hamiltonian is represented by a vector of Term objects. The bra and
//! ket forming the matrix element is defined in the constructor. The class
//! assumes the following relations are true:  
//! 
//!      - [J_n,J_m]= delta_{n+m,0} n  for all U(1) currents J,K,... 
//!      - [J_n,K_m]=0                 for all U(1) currents J,K,... 
//!		 - J_{n>0} |0> = 0  and <0| J_{n<0}   
//!      - J_0 |0> = N   
//!      - [Y_n,Y_m] = delta_{n+m,0}   for all Majorana operators (cahrged currents)
//!      = [Y_n,J_m] = 0
//!
////////////////////////////////////////////////////////////////////////////////

class ListOfTerms
{
	private:

	std::vector<Term> m_operator;     //!<   The Hamiltonian is represented by a vector of Term objects 

	bool FullyCommutedTest() const;
	void FullRightCommute();
	void ClearRightZeroTerms();
	void ClearZeroTerms();
	void EvaluateNumberOperators();
	ListOfTerms operator*=(const Term& rhs);
	void SortAllVectors();
	friend double InnerProduct(const Term&  bra,const ListOfTerms&  ket);
	
	public:

	//ListOfTerms(Term& basisBra,std::vector<Term>& op,Term& basisKet);
	//ListOfTerms(Term& basisBra,Term& basisKet);
	
	//	Default (empty) constructor
	ListOfTerms();
	
	//	Copy constructor
	
	ListOfTerms(const ListOfTerms& other);
	ListOfTerms(const std::vector<Term>& op);
	
	//  Destructor
	~ListOfTerms();

	void CombineSameTerms();
	std::vector<ListOfTerms> OperateOnKet(std::vector<Term>& basisKet);
	std::vector<Term> GetOperator();

	void PrintOperator() const;
	void PrintMatrixElement(const bool asKet) const;

	ListOfTerms operator+=(const Term& rhs);
	ListOfTerms operator+=(const ListOfTerms& rhs);
};

#endif
