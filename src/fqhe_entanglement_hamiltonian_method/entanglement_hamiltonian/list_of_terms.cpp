////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 11/04/2014
//!
//!  \file 
//!		The ListOfTerms class contains a list of Term objects
//!      and functions to perform commutations.
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

#include "list_of_terms.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Calculate the inner product between a Term object and a Virasoro object
//!	
//! \return Numerical value of inner product
//!
////////////////////////////////////////////////////////////////////////////////

double InnerProduct(
    const Term&  bra,            //!<    Bra vector of inner product
    const ListOfTerms&  ket)     //!<    Ket vector of inner product
{
    //  if the operators are the same (not including their coefficient)
    //  then we can exactly evaluate their inner product
    
    const int size = ket.m_operator.size();
    
    double innerProduct = 0.0;
    const double zeroTol = 0.0000000000000001;
    
    if(!fabs(ket.m_operator[0].GetCoefficient())<zeroTol)
    {
        //  Add terms to the inner product only if the ket is non-zero
    
        for(int i=0;i<size;i++)
        {
            innerProduct += InnerProduct(bra,(ket.m_operator[i]));        
        }
    }
    
    return innerProduct;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Default constructor for an object of type ListOfTerms. Generates
//! an empty class
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms::ListOfTerms()
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Constructor for an object of type ListOfTerms from a vector
//! of Term objects
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms::ListOfTerms(
	const std::vector<Term>& hamiltonian)     //!<    vector of Term objects representing the Hamiltonian
	: m_operator(hamiltonian)
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Copy constructor
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms::ListOfTerms(
	const ListOfTerms& other)      //!<    The instance of the class to copy
	:   m_operator(other.m_operator)    //      Simply generate a copy of the operator
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	Destructor for an object of type ListOfTerms. Currently empty
////////////////////////////////////////////////////////////////////////////////

ListOfTerms::~ListOfTerms()
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	Define an operation to combine two virasoro operators into one
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms ListOfTerms::operator+=(
    const ListOfTerms& rhs)    //!<    Operator to add to the one in the current object
{
    m_operator.insert(m_operator.end(),rhs.m_operator.begin(),rhs.m_operator.end());

    return *this;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	Define an operation to combine a virasoro operator and one additional 
//! term
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms ListOfTerms::operator+=(
    const Term& rhs)    //!<    Operator to add to the one in the current object
{
    m_operator.push_back(rhs);

    return *this;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief	Define an operation to combine a U(1) current operator and one additional 
//! term by multiplying each term the operator by the additional term
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms ListOfTerms::operator*=(
    const Term& rhs)
{
    for(unsigned int i=0;i<m_operator.size();i++)
    {
        m_operator[i] *= rhs;
    }
    
    return *this;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Prints out to the screen a representation of the currently stored
//!	matrix element.
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::PrintMatrixElement(const bool asKet) const
{
	utilities::cout.SecondaryOutput()<<std::endl;

	const int nbrTerms=m_operator.size();

	for(int i=0;i<nbrTerms;i++)
	{
        if(!asKet) 
        {
            utilities::cout.DebuggingInfo()<<"<";
            m_operator[i].PrintDeltaNa0();
            utilities::cout.DebuggingInfo()<<"| ";
        }
        
		m_operator[i].PrintTerm();

		utilities::cout.DebuggingInfo()<<"|";
		
		m_operator[i].PrintDeltaNa0();
		
		utilities::cout.DebuggingInfo()<<">";

		if(i<nbrTerms-1)
		{
			utilities::cout.DebuggingInfo()<<" + ";
		}
	}
	
	utilities::cout.DebuggingInfo()<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Prints out to the screen a representation of the currently stored
//!	operator
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::PrintOperator() const
{
	const int nbrTerms=m_operator.size();

	for(int i=0;i<nbrTerms;i++)
	{
		m_operator[i].PrintTerm();

		if(i<nbrTerms-1)
		{
			utilities::cout.SecondaryOutput()<<" + ";
		}
	}
	
	utilities::cout.SecondaryOutput()<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to multiply the stored Hamiltonian by a ket vector,
//!	returning a modified ket vector
//!
//!  \return    H |ket>
//!
////////////////////////////////////////////////////////////////////////////////

std::vector<ListOfTerms> ListOfTerms::OperateOnKet(std::vector<Term>& basisKet)
{
    //utilities::cout.DebuggingInfo()<<"OPERATE ON KET"<<std::endl;
    
	//	make basisKet.size() copies of the original Hamiltonian
	
	std::vector<ListOfTerms> storeValue(basisKet.size(),m_operator);

	for(unsigned int i=0;i<basisKet.size();i++)
	{
	    storeValue[i] *=  basisKet[i];  //	Multiply the Hamiltonian by each basis ket

		//storeValue[i].PrintMatrixElement(true);
		
		storeValue[i].FullRightCommute();		//	Evaluate |newKet> = H |ket>

		//storeValue[i].PrintMatrixElement(true);
	}
	
    //utilities::cout.DebuggingInfo()<<"DONE"<<std::endl;
    
	return storeValue;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to initiate the recursive commutation routine for ket
//!	 vectors.
//!
//!  This will completely evaluate the numerical value of the matrix element
//!  loaded into the object. Some short cuts are used for cases when the matrix
//!  element is known to be zero based on the analytic properties of the 
//!  commutation relations e.g J_1 J_2 J_-1 |0> = 0 because J_2 commutes with 
//!  all the other terms. In fact, most matrix elements tend to be zero.
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::FullRightCommute()
{
	//utilities::cout.DebuggingInfo()<<"\n\tFULL RIGHT COMMUTE:"<<std::endl;

	//	Evaluate e.g. J_0, K_0,... (number operator) terms first, since they commute with everything
	
	//utilities::cout.DebuggingInfo()<<"INITIAL THING:"<<std::endl;
	//this->PrintMatrixElement(true);
	
	this->EvaluateNumberOperators();

	//utilities::cout.DebuggingInfo()<<"AFTER NUMBER OPERATORS EVALUATED:"<<std::endl;
	//this->PrintMatrixElement(true);
	
	//	First remove terms that evaluate to zero i.e. with J_n>0 |0> or <0| J_n<0
	//	or terms where there are no pairs of n and -n (if that's not the case then
	//	a J_n could be pulled through to annihilate the vacuum).
	
	this->ClearRightZeroTerms();
	
	//utilities::cout.DebuggingInfo()<<"AFTER CLEARED ZERO TERMS:"<<std::endl;
	//this->PrintMatrixElement(true);

	bool doneTest = this->FullyCommutedTest();

	while(!doneTest)
	{
		//	For each term in the operator:

		std::vector<Term> newTerms;

		//	Take the right-most operator with positive label and commute right
		//	if it is already on the far right, then evaluate the term to be zero

		//	Alternatively, convert the result to the form
		//	n1! (n2!) 2^n2 (n3!) 3^n3 ...
		//	where n1 is the number of pairs of J_1 J_-1 etc.

		const int nbrTerms = m_operator.size();

        //utilities::cout.DebuggingInfo()<<"BEFORE COMMUTE RIGHT"<<std::endl;
		//this->PrintMatrixElement(true);
		
		for(int i=0;i<nbrTerms;i++)
		{
			m_operator[i].CommuteRight(newTerms);
		}
		
		//utilities::cout.DebuggingInfo()<<"AFTER COMMUTE RIGHT"<<std::endl;
		//this->PrintMatrixElement(true);

		//	Clear terms that are now zero
		
		this->ClearRightZeroTerms();
		
		//utilities::cout.DebuggingInfo()<<"AFTER CLEARED ZERO TERMS:"<<std::endl;
	    //this->PrintMatrixElement(true);
		
		//	update the operator

		m_operator.insert(m_operator.end(),newTerms.begin(),newTerms.end());

		//newTerms.erase(newTerms.begin(),newTerms.end());

        //utilities::cout.DebuggingInfo()<<"AFTER NEW TERMS ADDED"<<std::endl;
        //this->PrintMatrixElement(true);

		//	Combine any terms that are the same

		this->CombineSameTerms();
		
		//utilities::cout.DebuggingInfo()<<"RESULT SO FAR"<<std::endl;
		//this->PrintMatrixElement(true);getchar();

		//	Test if there is only a constant term remaining, otherwise keep on commuting...

		doneTest = this->FullyCommutedTest();
		
	}

	this->ClearRightZeroTerms();

    this->SortAllVectors();
    
    this->CombineSameTerms();

	//utilities::cout.DebuggingInfo()<<"FINAL THING:"<<std::endl;
	//this->PrintMatrixElement(true);
	
	//	return the final value of the matrix element

	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to determine whether each term will evaluate to zero or
//!	not when taking a ket of the term with the vacuum.
//!
//!	To do this, for each term we call ZeroTest from the Term class 
//! (see term.cpp). 
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::ClearRightZeroTerms()
{
	//utilities::cout.DebuggingInfo()<<"Clear Right Zero Terms"<<std::endl;
	//this->PrintMatrixElement(false);

	std::vector<Term> newTerms;

	if(m_operator.size()>0)
	{
		for(unsigned int i=0;i<m_operator.size();i++)
		{
			//	For each term check for 
			//	1. J_n |0> with n>0
			//	3. 	0 coefficient
			//	4. Pairs of J_n J_-n

			if(!m_operator[i].RightZeroTest())
			{
				newTerms.push_back(m_operator[i]);
			}
		}
	}

	//	make sure that we have at least one constant term (even if it is zero)

	if(newTerms.size()==0)
	{
		newTerms.push_back(Term(0.0,m_operator[0].m_occupations,m_operator[0].m_deltaNa0));
	}

	m_operator = newTerms;
	
	//this->PrintMatrixElement(false);getchar();

	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to combine terms that have the same operator content.
//!	We just generate a new term with the coefficients of the old like terms added.
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::CombineSameTerms()
{
	//utilities::cout.DebuggingInfo()<<"COMBINE TERMS"<<std::endl;
	//this->PrintOperator();

	std::vector<Term> newTerms;
	std::vector<int> blackList;

	const int nbrTerms=m_operator.size();

	for(int i=0;i<nbrTerms;i++)
	{
		//m_operator[i].PrintTerm();

		//utilities::cout.DebuggingInfo()<<"blackList contains:";

		int size=blackList.size();

		int j=0;
		bool proceed=true;

		while(j<size)
		{
			//utilities::cout.DebuggingInfo()<<blackList[j]<<" ";

			if(blackList[j]==i)
			{
				proceed=false;
				break;
			}

			j++;
		}

		//utilities::cout.DebuggingInfo()<<std::endl;

		if(proceed)
		{
			double coefficient=m_operator[i].GetCoefficient();

			for(int j=i+1;j<nbrTerms;j++)
			{
				if(m_operator[i].SameTest(m_operator[j]))
				{
					//utilities::cout.DebuggingInfo()<<"ADD TO BLACKLIST ";

					blackList.push_back(j);

					//utilities::cout.DebuggingInfo()<<blackList[0]<<std::endl;

					//	combine coefficients
					coefficient+=m_operator[j].GetCoefficient();
				}
			}

			newTerms.push_back(m_operator[i]);

			(newTerms.back()).SetCoefficient(coefficient);

			//(newTerms.back()).PrintTerm();
		}
	}

	m_operator = newTerms;

	//utilities::cout.DebuggingInfo()<<"nbr new terms: "<<newTerms.size()<<std::endl;
	//this->PrintOperator();

	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Returns the currently stored operator. 
//!
//! \return Currently stored operator as a vector of Term objects
//!
////////////////////////////////////////////////////////////////////////////////

std::vector<Term> ListOfTerms::GetOperator()
{
	return m_operator;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Evaluates the number operators e.g. J_0 for each term. 
//! 
//! Since these operators commute with everything, they can simply be
//!	replaced with the occupation numbers.	
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::EvaluateNumberOperators()
{
	const int nbrTerms=m_operator.size();

	if(nbrTerms>0)
	{
		for(int i=0;i<nbrTerms;i++)
		{
			m_operator[i].EvaluateNumberOperator();
		}
	}

	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Put the operator content of each term in ascending order	e.g.
//! J_-1 J_-2 J_-2 J_-3
//!
////////////////////////////////////////////////////////////////////////////////

void ListOfTerms::SortAllVectors()
{
    for(unsigned int i=0;i<m_operator.size();i++)
    {
        m_operator[i].SortVector();
    }
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Test if the operator is fully right commuted i.e. it contains no
//!	J_n>0 terms or is just a constant
//!
////////////////////////////////////////////////////////////////////////////////

bool ListOfTerms::FullyCommutedTest() const
{
    //utilities::cout.DebuggingInfo()<<"FULLY COMMUTED TEST"<<std::endl;

	const unsigned int nbrTerms = m_operator.size();

	//	It's fully commuted if the only remaining term is a constant

	if(nbrTerms==1)
	{
		if((m_operator[0].GetNbrU1())==0)
		{

			return true;
            //utilities::cout.DebuggingInfo()<<"FULLY COMMUTED"<<std::endl;
		}
	}
	
	//	It's fully commuted if every remaining term contains no J_n>0
	if(nbrTerms>0)
	{
		for(unsigned int i=0;i<nbrTerms;i++)
		{
			if(!m_operator[i].AllNegativeIndex())
			{
				return false;
			}
		}

		return true;
        //utilities::cout.DebuggingInfo()<<"FULLY COMMUTED"<<std::endl;
	}

	//	Otherwise it's not fully commuted
	
	return false;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
