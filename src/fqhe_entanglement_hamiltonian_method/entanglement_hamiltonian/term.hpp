////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		Header file for the Term class. Used to represent terms in an
//!		entanglement Hamiltonian written in terms of Virasoro algebra U(1)
//!		current operators.
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

#ifndef _TERM_HPP_INCLUDED_
#define _TERM_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/mathematics/combinatorics.hpp"
#include "../../utilities/algorithms/quick_sort.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include <vector>
#include <algorithm>    //  for std::reverse
#include <iostream>
#include <iomanip>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif
//	Some U(1) and charged current labels are defined to distinguish 
//	between different integer partition colours
#define _U1_LABEL_1_ 'J'    //!<    Define the first U(1) label
#define _U1_LABEL_2_ 'K'    //!<    Define a second U(1) label
#define _U1_LABEL_3_ 'I'    //!<    Define a third U(1) label
#define _MAJORANA_CURRENT_ 'Y'//!<  Also define a majorana current operator

////////////////////////////////////////////////////////////////////////////////
#if _BENCHMARK_MODE_==1	
//	Count of function calls
extern int g_callsCtor;
extern int g_callsCommuteRight;
extern int g_callsSameTest;
extern int g_callsZeroTest;
extern int g_callsCommutator;
#endif
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!	\brief The Term class contains data to represent a single term in the
//!	entanglement Hamiltonian, or a bra or ket vector in its Hilbert space.
//!			
//!	All terms are built from U(1) current operators in a Virasoro algebra.
//!
//! If Term is constructed with m_deltaNa0 assigned, then it represents a
//! bra or ket vector with reference to a particular particle cut. Alternatively,
//! it can be constructed without m_deltaNa0 assigned, in which case it 
//! represents an operator.
//!
//! Combinations of such objects will be combined in the correct way such
//! that e.g. operators acting on bras/ke4.08248e-01 K_-1 K_-1 K_-1 |-1 1 >
//!											
//!	Example:
//!		The operator 4 J_1 K_1 K_-2 J_-3 with N_A=4 would be represented as:
//!		- m_coefficient = 4.0
//!		- m_opArray = {1,1,-2,-3}
//!		- m_opType = {J,K,K,J}
//!		- m_deltaNa0 = {4}
////////////////////////////////////////////////////////////////////////////////
class Term
{
	private:
	double 	m_coefficient;	                //!< Numerical coefficient for the stored term
											//!
	std::vector <short int> m_opArray;		//!< A list of integers representing the U(1) indices	
											//!  of the U(1) operators making the term	
	std::vector <char> m_opType;		    //!< A list of labels of U(1) operators to
											//!  enable more than one type of U(1) operator
											//!  to be present
	std::vector <int> m_occupations;	    //!< Number of each constituent colour. This is 
											//!	 used to impose orthogonality between states 
											//!	 with different m_occupations.
	std::vector <int> m_deltaNa0;	    	//!< A list of integers to represent the shifted 
											//!  number of particles in the Na cut, labelled by
											//!	 their different colours. Given by the minimum of
											//!	 (N_A - N_A0) where N_A0 is any cut associated
											//!	 with the lowest pseudo-energy
	bool m_isOperator;  	                //!< Specify whether we're representing an operator 
											//!  or a bra/ket vector with our Term
	bool SameTest(const Term& rhs) const;
	bool ZeroTest() const;
	bool RightZeroTest() const;
	void CommuteRight(std::vector<Term>& newTerm);
	void CommuteLeft(std::vector<Term>& newTerm);
	void EvaluateNumberOperator();
	bool AllNegativeIndex() const;
	void Conjugate();
	bool IsAscendingForm();
	void ExactResult();
	void SortVector();
    double Commutator(const char label1, const int index1, const char label2, const int index2) const;
    void CountOperators(std::vector<int>& labelCounterPositive, std::vector<int>& labelCounterNegative,
                        std::vector<int>& shiftIndex, int* nbrColours, int* maxIndex) const;
	friend double InnerProduct(const Term& bra, const Term& ket);		
	friend class ListOfTerms;
	friend class HilbertSpace;														
	public:
	Term(const double coefficient, const std::vector<int>& occupations, const std::vector<int>& deltaNa0);	
	Term(const double coefficient);
	Term(const double coefficient, const std::vector<short int>& opArray, const std::vector<char>& opType);
	Term(const double coefficient, const std::vector<short int>& opArray, const std::vector<char>& opType,
		 const std::vector<int>& occupations, const std::vector<int>& deltaNa0);
	Term(const Term& other);
	~Term();
	Term operator*=(const Term& rhs);
	void PrintTerm() const;
	void PrintDeltaNa0() const;
	void SetCoefficient(double value);
	double GetCoefficient() const;
	int GetNbrU1() const;
};
#endif
