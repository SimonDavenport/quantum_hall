////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!		Library file for the Term class. Used to represent terms in an
//!		entanglement Hamiltonian written in terms of  U(1)
//!		current operators and Majorana operators. 
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

#include "term.hpp"

////////////////////////////////////////////////////////////////////////////////
//!	\brief Calculate the inner product between Term objects
//!	
//! \return Numerical value of inner product
////////////////////////////////////////////////////////////////////////////////
double InnerProduct(
    const Term&  bra,     //!<    Bra vector of inner product
    const Term&  ket)     //!<    Ket vector of inner product
{
    //  if the operators are the same (not including their coefficient)
    //  then we can exactly evaluate their inner product
    //  Make a copy of each vector, sort both to the same form
    //  and then check if they are the same or not (since the
    //  basis is orthogonal, if they differ then we get 0)
    Term tempBra(bra);
    Term tempKet(ket);
    tempBra.SortVector();
    tempKet.SortVector();
    if(tempBra.SameTest(tempKet))
    {
        //  Combine the bra and ket together
        tempBra.Conjugate();
        tempBra *= tempKet;
        //  Evaluate the operator exactly
        tempBra.ExactResult();
        return tempBra.GetCoefficient();
    }
    else
    {
        return 0.0;
    }
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Empty constructor for an object of type Term. Represents a bra/ket
//!	vector in the Hilbert space of an entanglement Hamiltonian built from U(1)
//!	current operators.
////////////////////////////////////////////////////////////////////////////////
Term::Term(
    const double coefficient,
    const std::vector<int>& occupations,
    const std::vector<int>& deltaNa)
 :   m_coefficient(coefficient),
	 m_opArray(0),
	 m_opType(0),
	 m_occupations(occupations),
	 m_deltaNa0(deltaNa),
	 m_isOperator(false)                 //  The object behaves like a bra or ket, 
										 //  rather than an operator
{
    ////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsCtor+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//!
//!	Constructor for the Term object, representing an operator term
//!
Term::Term(const double coefficient,const std::vector<short int>& opArray,const std::vector<char>& opType)
 :   m_coefficient(coefficient),
	 m_opArray(opArray),
	 m_opType(opType),
	 m_occupations(0),
	 m_deltaNa0(0),
	 m_isOperator(true)                 //  The object behaves like an operator
{
    ////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1
	g_callsCtor+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Full constructor for an object of type Term. Represents a bra/ket
//!	vector in the Hilbert space of an entanglement Hamiltonian built from U(1)
//!	current operators.
//!
//! The constructor sets all internal variables.
//!	An error is produced if the opArray and opType vectors are not the same size.
////////////////////////////////////////////////////////////////////////////////
Term::Term(
	const double coefficient,			    //!< Numerical coefficient for the stored term
	const std::vector<short int>& opArray,	//!< List of integers representing the U(1) indices						
	const std::vector<char>& opType,		//!< List of labels of U(1) operators
	const std::vector<int>& occupations,    //!< Number of particles of each colour
	const std::vector<int>& deltaNa0)	    //!< List of particle number of each colour
										    //!	 minus the particle number for the lowest
										    //!	 pseudo-energy cut
    :   m_coefficient(coefficient),
        m_opArray(opArray),
        m_opType(opType),
        m_occupations(occupations),
        m_deltaNa0(deltaNa0),
        m_isOperator(false)                 //  The object behaves like a bra or ket, 
                                            //  rather than an operator
{
	if(m_opArray.size() != m_opType.size())
	{
		std::cerr<<"TERM CONSTRUCTION ERROR: index list and label list do not match. "<<std::endl;
		getchar();
	}
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1
	g_callsCtor+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//!
//!	Copy constructor
//!
Term::Term(
    const Term& other)		//!< Term object to be copied
	
	:   m_coefficient(other.m_coefficient),
	    m_opArray(other.m_opArray),
	    m_opType(other.m_opType),
		m_occupations(other.m_occupations),
		m_deltaNa0(other.m_deltaNa0),
	    m_isOperator(other.m_isOperator) 
{
    ////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsCtor+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//!
//!	Constant term constructor - to represent a global shift operator
//!
Term::Term(const double coefficient)
 :   m_coefficient(coefficient),
	 m_opArray(0),
	 m_opType(0),
	 m_occupations(0),
	 m_deltaNa0(0),
	 m_isOperator(true)                 //  The object behaves like an operator
{

    ////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsCtor+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//!
//!	Destructor for an object of type Term. Currently empty.
//!	
Term::~Term()
{}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Calculates the commutator between the J_n and K_m U(1) current operators
//!	
//!	Here we are implementing the following commutation relations
//!		- [J_n,J_m] = n  delta_{m+n},0
//!		- [K_n,K_m] = n  delta_{m+n},0	... etc.
//!		- [J_n,K_m] = 0	... etc.
//!
//! Also evaluates commutation relations for charged operators Y_k as
//!
//!     - {Y_n,Y_m} = delta_{m+n},0
//!     - [Y_n,J_m] = 0 etc.
//!
//!	\return	Value of the commutator
////////////////////////////////////////////////////////////////////////////////
double Term::Commutator(
	const char label1,  //!< The colour label of the 1st U(1) operator (e.g. J, K etc.)
	const int  index1,  //!< The U(1) index of the 1st U(1) operator (e.g. J_index1)
	const char label2,  //!< The colour label of the 2nd U(1) operator (e.g. J, K etc.)
	const int  index2)  //!< The U(1) index of the 2nd U(1) operator (e.g. K_index2)
	const
{
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsCommutator+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
	if(label1==label2)
	{
		if(index1+index2==0)
		{
		    if(_MAJORANA_CURRENT_ == label1)
		    {
		        return 1.0;
		    }
		    else
		    {
			    return (double)index1;
			}
		}
		else
		{
			return 0.0;
		}
	}
	else
	{
		return 0.0;
	}
}

//!
//!	Prints a representation of the stored term to the command line.
//!
void Term::PrintTerm() const
{
	utilities::cout.SecondaryOutput() << std::scientific << std::setprecision(5) << m_coefficient<<" ";
	int length = m_opArray.size();
	for(int i=0; i<length; ++i)
	{
		utilities::cout.SecondaryOutput() << m_opType[i] << "_" << m_opArray[i] << " ";
	}
}

//!
//! Prints a representation of the shifted particle cut associated with 
//!	the current state
//!
void Term::PrintDeltaNa0() const
{
	int length = m_deltaNa0.size();
	for(int i=0; i<length; ++i)
	{
		utilities::cout.SecondaryOutput() << m_deltaNa0[i] << " ";
	}
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Returns the current value of the term's coefficient.
//!
//!	\return	The current value of the term's coefficient
////////////////////////////////////////////////////////////////////////////////
double Term::GetCoefficient() const
{
	return m_coefficient;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief	Returns the current number of U(1) operators in the term.
//!
//!	\return	The current number of U(1) operators in the term
////////////////////////////////////////////////////////////////////////////////
int Term::GetNbrU1() const
{
	return m_opArray.size();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief	Determines if the term contains only J_m<0
//!
//!	\return	true if there are only J_m<0 in the term
////////////////////////////////////////////////////////////////////////////////
bool Term::AllNegativeIndex() const
{
	int signCounter = 0;
	const int size = m_opArray.size();
	for(int i=0; i<size; ++i)
	{
		if(m_opArray[i]<0)	
		{
		    ++signCounter;
		}
	}
	return signCounter == size;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief	Resets the current value of the term's coefficient to a new value.
//!
//!	This might be used when adding up terms of the same type.
////////////////////////////////////////////////////////////////////////////////
void Term::SetCoefficient(
	double value)	//!<	New value of the coefficient
{
	m_coefficient = value;
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Overloads the * operator such that objects of type Term 'multiply'
//!	together in the usual way i.e. the coefficients multiply and the operator
//!	list gets concatenated.
//!
//!  We want "multiplication" to act such that for "this" * rhs, if "this" is an operator
//!  and rhs is a ket then "this" picks up the m_occupations value from the rhs object and
//!  then the result "lhs" is a ket vector.
//!  If "this" is a bra and rhs is an operator then rhs should pick up the m_occupations
//!   value of "this" and the result "lhs" is a bra vector.
//!  If "this" and rhs are bra and ket then we shall define the result as orthogonal
//!  unless their m_occupations values are the same.
//!
//!	\return A new object of type Term which is a combination of the multiple objects. 	
////////////////////////////////////////////////////////////////////////////////
Term Term::operator*=(
    const Term& rhs)	//!<	Right-hand side of the expression as this = this * rhs
{
	m_coefficient *= rhs.m_coefficient;
	m_opArray.insert(m_opArray.end(), rhs.m_opArray.begin(), rhs.m_opArray.end());
	m_opType.insert(m_opType.end(), rhs.m_opType.begin(), rhs.m_opType.end());
    if(m_isOperator && !rhs.m_isOperator)
    {
        //  If "this" is an operator and rhs is a ket then "this" picks up the 
        //  m_occupations value from the rhs object and then the result "lhs" is 
        //  a ket vector
        m_occupations = rhs.m_occupations;
        m_deltaNa0    = rhs.m_deltaNa0;
        m_isOperator = false;     //  It's now a ket
    }
    else if(!m_isOperator && rhs.m_isOperator)
    {
        //  If "this" is a bra and rhs is an operator then rhs should pick up the 
        //  m_occupations value of "this" and the result "lhs" is a bra vector.    
	}
	else if(!m_isOperator && !rhs.m_isOperator)
	{
	    //  If "this" and rhs are bra and ket then we shall define the result as 
	    //  orthogonal unless their m_occupations values are the same.
	    if(m_occupations != rhs.m_occupations)
	    {
	        m_coefficient=0.0;
	    }
	}
	else
	{
	    //  Complain if something confusing happened
	    if(m_occupations != rhs.m_occupations)
	    {
	        std::cerr << "ERROR with Term::operator* : Illegal operator combination!" << std::endl;
	        exit(EXIT_FAILURE);
	    }
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Takes the right-most operator with positive index and commutes it to the
//!	right. This function takes into account that we are evaluating a matrix element
//!	with the vacuum state.  
//!
//! This is done regardless of the operator label (e.g. J or K).
//!	If the right-most operator is already on the far right, then do nothing 
//!	(it will be removed by a call to the 'ZeroTest' function)
//!
//!	When the commutation takes place, it will potentially produce two terms.
//!	The original uncommuted term will be overwritten with operators J_n and
//!	J_m reversed. The additional term from the commutator [J_n,J_m] will be
//!	appended to the list 'newTerm', as long as it is non-zero.
//!	
//!	Different commutation relations can be called depending on the U(1)
//!	current label e.g. [J_n,k_m] can be different to [J_n,J_m] etc.
////////////////////////////////////////////////////////////////////////////////
void Term::CommuteRight(
	std::vector<Term>& newTerm)	//!<	A list where new terms can be added if they are non-zero	
{
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsCommuteRight+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
	int nbr = m_opArray.size();
	// If the array is empty (just a constant) then return
	if(nbr==0)
	{
		return;
	}
	//	If the term is of "ascending form" and it contains an even number of operators
	//	then return the analytical result
	if(((nbr&1)==0) && (IsAscendingForm()))
	{
		this->ExactResult();
		return;
	}
	//////    Determine the right-most operator with positive index    /////////////
	int toCommute=nbr-1;
	while(toCommute>0)
	{
		if(m_opArray[toCommute]>0)
		{
			break;
		}
		toCommute--;
	}
	//////    Do nothing if the right-most operator is already on the far right 
	if(toCommute==nbr-1)
	{
		return;
	}
	else	//	Commute it to the right using J_n J_m = J_m J_n + [J_n,J_m]
	{
		//	Determine how many places we can commute the operator through "for free"
		//	exploiting the commutation relation being 0 when n+m !=0;
		//	From the position toCommute, search right until an operator is found which will give a non-zero
		//	commutation. Otherwise, move all the way to the right.
		int commuteWith=toCommute+1;
		while(m_opArray[commuteWith]+m_opArray[toCommute]!=0)
		{
			++commuteWith;
			if(commuteWith>=nbr)
			{
				commuteWith = nbr-1;
				break;
			}
		}
		// 	Generate a possible term given by the [J_n,J_m] part
		double commutator = Commutator(m_opType[toCommute], m_opArray[toCommute], m_opType[commuteWith], m_opArray[commuteWith]);
		if(commutator!=0)
		{
			Term returnTerm(*this);		//	make a copy of the original term
			//	remove the two U(1) current operators that we replaced by a commutator
			returnTerm.m_opArray.erase(returnTerm.m_opArray.begin()+toCommute);
			returnTerm.m_opType.erase(returnTerm.m_opType.begin()+toCommute);
			returnTerm.m_opArray.erase(returnTerm.m_opArray.begin()+commuteWith-1);
			returnTerm.m_opType.erase(returnTerm.m_opType.begin()+commuteWith-1);
			//	multiply the m_coefficient by the value of the commutator
			returnTerm.m_coefficient *= commutator;
			newTerm.push_back(returnTerm);
		}
		//	Overwrite the original term swapping J_n and J_m:
        //  TODO implement - sign in the case of fermion modes
		int tempInt = m_opArray[toCommute];
		m_opArray[toCommute] = m_opArray[commuteWith];
		m_opArray[commuteWith] = tempInt;
		char tempChar = m_opType[toCommute];
		m_opType[toCommute] = m_opType[commuteWith];
		m_opType[commuteWith] = tempChar;
	}
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Test if two terms are the 'same', meaning that they contain exactly
//!	the	 same list of U(1) current operators, regardless of the coefficient.
//!	(coefficients of such identical terms will be added together)
//!
//!	Works as a.SameTest(b) returns true if a and b are the same according
//!	to the above definition.
//!
//!	\return	true/false depending on comparison
////////////////////////////////////////////////////////////////////////////////
bool Term::SameTest(
	const Term& rhs)	//!<	Term object to be compared with
    const 
{
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsSameTest+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
	if(m_occupations == rhs.m_occupations)
	{
	    //  Sort the arrays in preparation for comparison
		if(m_opArray==rhs.m_opArray)
		{
			if(m_opType==rhs.m_opType)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to determine whether the current term will evaluate to zero
//!	when taking a matrix element of the term with the vacuum.
//!
//!	There are 5 possible cases:	
//!		- If there are an odd number of operators
//!		- If the first term is such that we have <0| J_n<0	...
//!		- If the last term is such that we have ... J_n>0 |0>
//!		- If the m_coefficient is zero (due to e.g. added terms cancelling)
//!		- If there are no pairs of J_n J_-n (otherwise a single one would commute throught and annihilate the vacuum).  
//!
//!	\return If the term is zero then return true, otherwise false
////////////////////////////////////////////////////////////////////////////////
bool Term::ZeroTest() const
{
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsZeroTest+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
	const int nbrOps = m_opArray.size();
	//	the matrix element of the term will be zero if it contains an odd number of operators
	if((nbrOps & 1)==1)
	{
		return true;
	}
	//	The matrix element will be zero if it has a vanishing coefficient
	const double zeroTol = 0.0000000000000001;
	if(fabs(m_coefficient)<zeroTol)
	{
		return true;
	}
	//	The matrix element will be zero if J_n>0 |0> or <0|J_n<0 occurs in it
	if(nbrOps>0 && (m_opArray.front()<0 || m_opArray.back()>0) )
	{
		return true;
	}
	//	Check for pairs of +n and -n
	if(nbrOps>0)
	{
		std::vector<int> labelCounterPositive;
        std::vector<int> labelCounterNegative;
        std::vector<int> shiftIndex;
        int nbrColours;
        int maxIndex;
	    this->CountOperators(labelCounterPositive, labelCounterNegative, shiftIndex, &nbrColours, &maxIndex);
		//	Check to see if the number of J_n<0 and J_n>0 match up for each n
		for(int i=0; i<nbrColours*maxIndex; ++i)
		{
			if(labelCounterPositive[i] != labelCounterNegative[i])
			{
				return true;
			}
		}
	}
	//	otherwise it is non-zero
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to determine whether the current term will evaluate to zero
//!	when evaluating a ket vector
//!
//!	There are 3 possible cases:	
//!		- If the last term is such that we have ... J_n>0 |0>
//!		- If the m_coefficient is zero (due to e.g. added terms cancelling)
//!		- Not every operator J_n>0 can be paired with a J<n to its right
//!
//!	\return If the term is zero then return true, otherwise false
////////////////////////////////////////////////////////////////////////////////
bool Term::RightZeroTest() const
{
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	g_callsZeroTest+=1;
	#endif
	////////////////////////////////////////////////////////////////////////////////
	const int nbrOps=m_opArray.size();
	//	The matrix element will be zero if it has a vanishing coefficient
	const double zeroTol = 0.0000000000000001;
	if(fabs(m_coefficient)<zeroTol)
	{
		return true;
	}
	//	The matrix element will be zero if J_n>0 |0>  occurs in it
	if(nbrOps>0 && m_opArray.back()>0 )
	{
		return true;
	}
	//	Check that the number of J_n>0 is not greater than the number of J_n<0 
	//	for a given n and colour
	if(nbrOps>0)
	{
	    std::vector<int> labelCounterPositive;
        std::vector<int> labelCounterNegative;
        std::vector<int> shiftIndex;
        int nbrColours;
        int maxIndex;
	    this->CountOperators(labelCounterPositive, labelCounterNegative, shiftIndex, &nbrColours, &maxIndex);
		//	Check to see if the number of J_n<0 is equal to or exceeds the number
		//	of J_n>0 for the same n and colour 
		for(int i=0; i<nbrColours*maxIndex; ++i)
		{
			if(labelCounterPositive[i]>labelCounterNegative[i])
			{
				return true;
			}
		}
	}
	//	Additionally check that every J_n>0 there is a J_n<0 to its right,
	//	of the same colour
	if(nbrOps>0)
	{
		for(int i=0; i<nbrOps; ++i)
		{
			short int value = m_opArray[i];
			char type = m_opType[i];
			if(value>0)
			{
				//	Check if there is a corresponding n<0 term to its left
				for(int j=i+1; j<nbrOps; ++j)
				{
					if(m_opType[j] == type && m_opArray[j]+value==0)
					{
						//	it's non-zero
						return false;
					}
				}
				return false;
			}
		}
	}
	//	otherwise it is non-zero
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief Take the Hermitian conjugate of the Hamiltonian. 				 
//!
//! This swaps the operator ordering and reverses the sign of all the operators.
////////////////////////////////////////////////////////////////////////////////
void Term::Conjugate()
{
	std::reverse(m_opArray.begin(), m_opArray.end());
	std::reverse(m_opType.begin(), m_opType.end());
	int length = m_opArray.size();
	for(int i=0; i<length; ++i)
	{
		m_opArray[i] =- m_opArray[i];
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Evaluates the number operators e.g. J_0 for each term. 
//!
//! Since these operators commute with everything, they can simply
//!	be replaced with the shifted number of particles in the cut of a given colour, 
//!	divided by the square root of the filling factor. NOTE: if we only consider
//!	CF-EWF states, then the effective filling factor is nu=1 (a single LL of a 
//!	given colour).
//!
//!	After doing this, we can remove the number operator(s) from the term.
//!
//! Note the number operator contains a factor of 1/sqrt(filling factor), however
//! that can be factored out, so it must be included in the original definition
//! of the number operator term's coefficient. 							
////////////////////////////////////////////////////////////////////////////////
void Term::EvaluateNumberOperator()
{
	if(m_deltaNa0.size()==0)
	{
		std::cerr << "\n\tERROR in EvaluateNumberOperator() : m_deltaNa0 not assigned!" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::vector<short int> arr;
	std::vector<char> type;
	int nbrOps = m_opType.size();
	if(nbrOps>0)
	{
		for(int i=0; i<nbrOps; ++i)
		{
			if(m_opArray[i] == 0)
			{
			    if(m_opType[i] == _U1_LABEL_1_)
			    {
			        m_coefficient *= m_deltaNa0[0];
			    }
			    else if(m_opType[i] == _U1_LABEL_2_)
			    {
			        m_coefficient *= m_deltaNa0[1];
			    }
			    else if(m_opType[i] == _U1_LABEL_3_)
			    {
			        m_coefficient *= m_deltaNa0[2];
			    }
			    //	NOTE: extend this if statement if more colours are required
			}
			else
			{
				arr.push_back(m_opArray[i]);
				type.push_back(m_opType[i]);
			}
		}
	}
	m_opType  = type;
	m_opArray = arr;
	return;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Test whether the term is in "ascending form" meaning that all 
//! J_n<0 operators appear to the right of all J_n>0 operators
//!
//! If the term is in ascending form then a matrix element can be evaluated
//! using a simple analytical expression.
//!
//! \return true if in "ascending form", false otherwise
////////////////////////////////////////////////////////////////////////////////
bool Term::IsAscendingForm()
{
    int nbrOps = m_opType.size();
    if(nbrOps>0)
	{
	    int nbrColours = 4;
        std::vector<char> possibleColours(nbrColours);
        possibleColours[0] = _U1_LABEL_1_;
        possibleColours[1] = _U1_LABEL_2_;
        possibleColours[2] = _U1_LABEL_3_;
        possibleColours[3] = _MAJORANA_CURRENT_;
	    for(int c=0; c<nbrColours; ++c)
	    {
	        int signChanges = 0;
		    for(int i=0; i<nbrOps; ++i)
		    {
		        if(m_opType[i] == possibleColours[c])
		        {
		            if(m_opArray[i]>0 && signChanges == 0)
		            {
                        continue;
		            }
		            else if(m_opArray[i]<0 && signChanges == 1)
		            {
		                continue;
		            }
		            else if(m_opArray[i] == 0)
		            {
		                continue;
		            }
		            else if(signChanges>1)
		            {
                        return false;
		            }
		            else
		            {
		                signChanges+=1;
		            }
		        }
		    }
		    if(signChanges>1)
            {
                ++signChanges;
                return false;
            }
        }
    }
    else
    {
        //  Just a constant term on its own is in ascending form
        return true;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Evaluates the exact result for the case when the term is in "ascending
//!	form" (see the function IsAscendingForm for the definition)
//!
//!	The coefficient of this term is then simply given by
//!	n1! (n2!) 2^n2 (n3!) 3^n3 ...
//!	where n1 is the number of pairs of J_1 J_-1 etc. of each colour
//!	(ignoring J_0 terms here)
//!
//! In addition, for the charged modes we get the factor
//!
//! n1! n2! n3! ...
//!
//! TODO implement correct sign in the case of majorana modes
////////////////////////////////////////////////////////////////////////////////
void Term::ExactResult()
{
	//	We need to keep track of the counting of each type of
	//	label (colour, and index combined)
	//	To do that we need to know how many different types of
	//	these terms occur
	std::vector<int> labelCounterPositive;
    std::vector<int> labelCounterNegative;
    std::vector<int> shiftIndex;
    int nbrColours;
    int maxIndex;
	this->CountOperators(labelCounterPositive, labelCounterNegative, shiftIndex,
	                     &nbrColours, &maxIndex);
	//  If the matrix element/ket does not contain an equal number of
	//  J_n>0 and J_n<0 of the same type, then exact evaluation is
	//  invalid
	for(unsigned int i=0; i<labelCounterPositive.size(); ++i)
	{
	    if(labelCounterPositive[i] != labelCounterNegative[i])
	    {
	        return;
	    }
	}
	//	Now that we have the count, simply return the combined result
	long int result = 1;
	bool usedLabel1 = false;
    bool usedLabel2 = false;
    bool usedLabel3 = false;
    bool usedLabel4 = false;
    for(unsigned int i=0; i<m_opArray.size(); ++i)
    {
	    if(m_opType[i]==_U1_LABEL_1_)		usedLabel1 = true;
	    if(m_opType[i]==_U1_LABEL_2_)		usedLabel2 = true;
	    if(m_opType[i]==_U1_LABEL_3_)		usedLabel3 = true;
	    if(m_opType[i]==_MAJORANA_CURRENT_) usedLabel4 = true;
    }
	for(int i=0; i<maxIndex; ++i)
	{
	    if(usedLabel1)
		{
	        result *= utilities::Factorial<double>(labelCounterPositive[i])*
	                                               pow(i+1, labelCounterPositive[i]);
	    }
		if(usedLabel2)
		{
			result *= utilities::Factorial<double>(labelCounterPositive[shiftIndex[0]+i])*
			                                       pow(i+1, labelCounterPositive[shiftIndex[0]+i]);
		}
		if(usedLabel3)
		{
			result *= utilities::Factorial<double>(labelCounterPositive[shiftIndex[1]+i])*
			                                       pow(i+1,labelCounterPositive[shiftIndex[1]+i]);
		}
		if(usedLabel4)
		{
			result *= utilities::Factorial<double>(labelCounterPositive[shiftIndex[2]+i]);
		}
	}
	m_opArray.clear();
	m_opType.clear();
	m_coefficient *= result;
	return;
}

//!
//! Put the operator content of the term in ascending order e.g.
//! J_-1 J_-2 J_-2 J_-3
//!
void Term::SortVector()
{
    utilities::QuickSort<short int, char, _ASCENDING_ORDER_>(&m_opArray[0], &m_opType[0], m_opType.size());
    std::reverse(m_opArray.begin(), m_opArray.end());
    std::reverse(m_opType.begin(), m_opType.end());
    //  If there is more than one colour, then also sort by colour
    if(m_occupations.size()>1)
    {
        utilities::QuickSort<char, short int, _ASCENDING_ORDER_>(&m_opType[0], &m_opArray[0], m_opType.size());
    }
    return;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Make arrays containing the count of each operator J_k for a given
//! colour and k index. Also return the total number of colours used
//! and the maximum k index appearing (these numbers define the array dimensions
//! of labelCounterPositive and labelCounterNegative).
////////////////////////////////////////////////////////////////////////////////
void Term::CountOperators(
    std::vector<int>& labelCounterPositive, //!<    On completion contains the counts
                                            //!     of all J_k with k positive
    std::vector<int>& labelCounterNegative, //!<    On completion contains the counts
                                            //!     of all J_k with k negative
    std::vector<int>& shiftIndex,           //!<    List of starting positions
                                            //!     of different colour operator countrs
                                            //!     in the labelCounter outputs
    int* nbrColours,                        //!<    On completion returns the number
                                            //!     of different colours used
    int* maxIndex)                          //!<    On completion returns the maximum
                                            //!     k index appearing
    const
{	
	const unsigned int nbrOps = m_opArray.size();
    int usedLabel1 = 0;
    int usedLabel2 = 0;
    int usedLabel3 = 0;
    int usedLabel4 = 0;
    *maxIndex      = 0;
    shiftIndex.resize(3);
    for(unsigned int i=0; i<nbrOps; ++i)
    {
	    if(abs(m_opArray[i])>*maxIndex)	    *maxIndex=abs(m_opArray[i]);
	    if(m_opType[i]==_U1_LABEL_1_)		usedLabel1 = 1;
	    if(m_opType[i]==_U1_LABEL_2_)		usedLabel2 = 1;
	    if(m_opType[i]==_U1_LABEL_3_)		usedLabel3 = 1;
	    if(m_opType[i]==_MAJORANA_CURRENT_) usedLabel4 = 1;
    }
    *nbrColours = usedLabel1 + usedLabel2 + usedLabel3 + usedLabel4;
    labelCounterNegative.resize(*nbrColours*(*maxIndex));
    labelCounterPositive.resize(*nbrColours*(*maxIndex));
    for(int i=0; i<*nbrColours*(*maxIndex); ++i)
    {
	    labelCounterNegative[i]=0;
	    labelCounterPositive[i]=0;
    }
    shiftIndex[0] = usedLabel1*(*maxIndex);
    shiftIndex[1] = usedLabel1*(*maxIndex)+usedLabel2*(*maxIndex);
    shiftIndex[2] = usedLabel1*(*maxIndex)+usedLabel2*(*maxIndex)+usedLabel3*(*maxIndex);
    for(unsigned int i=0; i<nbrOps; ++i)
    {
	    if(m_opType[i] == _U1_LABEL_1_)
	    {
		    if(m_opArray[i]<0)
		    {
			    ++labelCounterNegative[abs(m_opArray[i])-1];
		    }
		    else if(m_opArray[i]>0)
		    {
			    ++labelCounterPositive[abs(m_opArray[i])-1];
		    }
	    }
	    if(m_opType[i] == _U1_LABEL_2_)
	    {
		    if(m_opArray[i]<0)
		    {
			    ++labelCounterNegative[shiftIndex[0]+abs(m_opArray[i])-1];
		    }
		    else if(m_opArray[i]>0)
		    {
			    ++labelCounterPositive[shiftIndex[0]+abs(m_opArray[i])-1];
		    }
	    }
	    if(m_opType[i] == _U1_LABEL_3_)
	    {
		    if(m_opArray[i]<0)
		    {
			    ++labelCounterNegative[shiftIndex[1]+abs(m_opArray[i])-1];
		    }
		    else if(m_opArray[i]>0)
		    {
			    ++labelCounterPositive[shiftIndex[1]+abs(m_opArray[i])-1];
		    }
	    }
	    if(m_opType[i] == _MAJORANA_CURRENT_)
	    {
		    if(m_opArray[i]<0)
		    {
			    ++labelCounterNegative[shiftIndex[2]+abs(m_opArray[i])-1];
		    }
		    else if(m_opArray[i]>0)
		    {
			    ++labelCounterPositive[shiftIndex[2]+abs(m_opArray[i])-1];
		    }
	    }
    }
    return;
}
