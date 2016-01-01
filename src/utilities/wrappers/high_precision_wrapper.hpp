////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/12/2013
//!
//!  \file
//!		This file contains a list of operator overloads that can be called
//!		when using different types of real or complex high precision variable
//!		It is based on the template wrapper class template 
//!		<typename U,int P> class HpWrap
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

#ifndef _HIGH_PREC_WRAPPER_HPP_INCLUDED_
#define _HIGH_PREC_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "../general/dcmplx_type_def.hpp"  // type dcmplx as complex double type

//	Include high precision c-libraries
#include <gmp.h>        //  See http://gmplib.org/
#include <stdint.h>     //  Needed by mpc.h
#include <mpfr.h>       //  See http://www.mpfr.org/
#include <mpc.h>        //  See http://www.multiprecision.org/

//////      STATIC CONSTANT AND TYPEDEF DECLARATIONS       ////////////////////

#ifndef _MPF_COMPLX_TYPE_DEFINED_
#define _MPF_COMPLX_TYPE_DEFINED_
//!
//!	Define a special complex mpf wrapper type
//!	NOTE this only acts as a place holder to describe complex mpf variables
//!
typedef std::complex<mpf_t*> mpfCmplx;
#endif

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//


////////////////////////////////////////////////////////////////////////////////
//! \brief A namespace to contain any functions and utilities that I have 
//! written for use with any c++ program.
//!
////////////////////////////////////////////////////////////////////////////////

namespace utilities
{

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief A class template to implement the Kahan summation algorithm for 
    //! objects of the high precision wrapper type
    //!
    //////////////////////////////////////////////////////////////////////////////// 

    template<class T>
    struct HighPrecKahanAccumulation
    {
        T m_sum;                //!<    Current accumulated sum 
        T m_correction;         //!<    Correction term
        T m_tempSum;
        T m_tempCorrection;
        
        void Init()
        {
            m_sum.Set(0.0);
            m_correction.Set(0.0);
        }
        
        T& GetTotal()
        {
            return m_sum;
        }
    };

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief Defines a global  template += operator overload for the Kahan 
    //! accumulator class. Call as simply a += b where a is an existing 
    //! HighPrecKahanAccumulation object and b is the numerical value to be added to 
    //! the sum (should be of the same type as the class type)
    //!
    //! \return address to the updated HighPrecKahanAccumulation object
    //!
    //////////////////////////////////////////////////////////////////////////////// 
    
    template<typename T>
    HighPrecKahanAccumulation<T>& operator+=(HighPrecKahanAccumulation<T>& a,const T& b)
    {
        a.m_tempCorrection = b;
        a.m_tempSum = a.m_sum;
        
        a.m_tempCorrection -= a.m_correction;
        a.m_tempSum += a.m_tempCorrection;

        a.m_correction = a.m_tempSum;
        a.m_correction -= a.m_sum;
        a.m_correction -= a.m_tempCorrection;

        a.m_sum = a.m_tempSum;
        
        //std::cout<<"m_sum  = "<<a.m_sum.Get()<<std::endl;
        //std::cout<<"m_correction  = "<<a.m_correction.Get()<<std::endl;
        
        return a;
    }
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

    //////////////////////////////////////////////////////////////////////////////// 
    //! \brief Defines a global  template *= operator overload for the Kahan 
    //! accumulator class. Call as simply a *= b where a is an existing 
    //! HighPrecKahanAccumulation object and b is the numerical value to be multiplied
    //! by (should be of the same type as the class type)
    //!
    //! \return address to the updated HighPrecKahanAccumulation object
    //!
    //////////////////////////////////////////////////////////////////////////////// 
    
    template<typename T>
    HighPrecKahanAccumulation<T>& operator*=(HighPrecKahanAccumulation<T>& a,const T& b)
    {
        a.m_sum *= b;
        a.m_correction *= b;
        
        //std::cout<<"m_sum  = "<<a.m_sum.Get()<<std::endl;
        //std::cout<<"m_correction  = "<<a.m_correction.Get()<<std::endl;

        return a;
    }

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	A generic template class HpWrap<complex type,int>
//!	
//!	This class is used to declare function overloads for +=, -=, *=, /= and = 
//!	(most of which are not declared for e.g. complex<double> variables)	
//!	it also defines multiplication of a complex number by an integer
//!
//!	The following specializations are defined
//!	
//!	HpWrap<complex<double>,0> for complex doubles
//!
//!	HpWrap<complex<__float128>,0> complex __float128 (from quadmath.h)	
//!
//!	HpWrap<mpc_t,prec> complex arbitary precision (mpc_t)
//!								 precision is set with the integer 'prec'
//!
//!	cCmplexWrap<mpfr_t,prec> real arbitary precision (mpfr_t)
//!								 precision is set with the integer 'prec'
//!
//!		cCmplexWrap<__float128>  real __float128  (from quadmath.h)	
//!	
//!		Notes
//!	
//!		Each call to the class template will instruct the compiler to generate
//!		A new class. So HpWrap<mpc_t,128> and HpWrap<mpc_t,256>
//!		will be complied as separate classes (that's how a template works!)
//!	
////////////////////////////////////////////////////////////////////////////////

template <typename U,int P>
class HpWrap
{
	private:

	//!	\brief	Internal data type variable
	//!
	U m_x;

	public:

	//!
	//!	\brief	Overload for the += operator
	//!
	HpWrap& operator+= (const HpWrap& rhs)
	{
		m_x=m_x+rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Overload for the -= operator
	//!
	HpWrap& operator-= (const HpWrap& rhs)
	{
		m_x=m_x-rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Overload for the *= operator
	//!
	HpWrap& operator*= (const HpWrap& rhs)
	{
		m_x=m_x*rhs.m_x;
		return *this;
	}

	//!
	//!	\brief	Overload for the /= operator
	//!
	HpWrap& operator/= (const HpWrap& rhs)
	{
		m_x=m_x/rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Overload for the = operator
	//!
	HpWrap& operator= (const HpWrap& rhs)
	{
		m_x=rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Overload the < operator
	//!
	bool operator< (const HpWrap& rhs) const
	{
		return (m_x.real()<rhs.m_x.real());
	}
	
	//!
	//!	\brief	Multiplication by some other type
	//!
	HpWrap& MulBy (int lhs,const HpWrap& rhs)
	{
		m_x=U(lhs*rhs.m_x.real(),lhs*rhs.m_x.imag());
		return *this;
	}
	
	//!
	//!	\brief	Define a function to multiply by long int
	//!
	HpWrap& MulBy (long int lhs,const HpWrap& rhs)
	{
		m_x=U(lhs*rhs.m_x.real(),lhs*rhs.m_x.imag());

		return *this;
	}
	
	//!
	//!	\brief	Overload function MulBy for double* type
	//!
	HpWrap& MulBy (double* lhs,const HpWrap& rhs)
	{
		m_x=*lhs*rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Overload function MulBy for double type
	//!
	HpWrap& MulBy (double lhs,const HpWrap& rhs)
	{
		m_x=lhs*rhs.m_x;
		return *this;
	}
	
	//!
	//!	\brief	Divide by an integer
	//!
	HpWrap& DivBy (int lhs,const HpWrap& rhs)
	{
		m_x=U(rhs.m_x.real()/lhs,rhs.m_x.imag()/lhs);
		return *this;
	}
	
	//!
	//!	\brief Define a log function
	//!
	HpWrap& LogOf (const HpWrap& rhs)
	{
		m_x=log(rhs.m_x);
		return *this;
	}

	//!
	//!	\brief Define an abs function
	//!
	HpWrap& AbsOf (const HpWrap& rhs)
	{
		m_x=U(abs(rhs.m_x),0.0);
		return *this;
	}

	//!
	//!	\brief Set the internal variable to some external value
	//!	redeclared in specializations

	void Set(dcmplx in){m_x=(U)in;}
	
	//!
	//!	\brief Set the internal variable to some external value
	//!	redeclared in specializations
	void Set(double in){m_x=(U)in;}

	//!
	//!	\brief	Return the value of the internal variable as a dcmplx type.
	//!	Redeclared in specializations
	//!

	dcmplx Get() const {return (dcmplx)m_x;}

	//!
	//!	\brief	Return the value of the internal variable
	//!

	U* GetReference() const {return &m_x;}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Specialization of the cpWrap template class to double variables
//!	
////////////////////////////////////////////////////////////////////////////////

template <int P> class HpWrap<double,P>
{
	private:

	//	internal double variable
	double m_x;

	public:

	HpWrap()
    {
        m_x = 0;
    }

	~HpWrap(){}

	//!
	//!	\brief	Overload for the *= operator
	//!
	HpWrap& operator*= (const HpWrap& rhs)
	{
		m_x=m_x*rhs.m_x;
		return *this;
	}

	//!
	//!	\brief	Overload for the /= operator
	//!
	HpWrap& operator/= (const HpWrap& rhs)
	{
		m_x=m_x/rhs.m_x;
		return *this;
	}

	//!
	//!	\brief	Overload function MulBy for double type
	//!
	HpWrap& MulBy (double lhs,const HpWrap& rhs)
	{
		m_x=lhs*rhs.m_x;
		return *this;
	}

	//!	
	//!	\brief Set the internal variable to some external value
	//!	redeclared in specializations
	void Set(dcmplx in){m_x=in.real();}
	
	//!	
	//!	\brief Set the internal variable to some external value
	//!	redeclared in specializations
	void Set(double in){m_x=in;}

	//!
	//!	\brief	Return the value of the internal variable
	//!
	double Get() const {return m_x;}

	//!
	//!	\brief	Return a pointer to the internal variable
	//!
	
	double* GetReference(){return &m_x;}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Specialization of the cpWrap template class to mpc_t variables
//!	
////////////////////////////////////////////////////////////////////////////////

template <int P> class HpWrap<mpc_t,P>
{
	private:

	//	internal mpc_t variable
	mpc_t m_x;

	public:

	//	constructor initializes with precision P
	HpWrap(){mpc_init2(m_x,P);}

	//	destructor clears variable to prevent all memory leaks
	~HpWrap(){mpc_clear(m_x);}

	//	overload for the += operator
	HpWrap& operator+= (const HpWrap& rhs)
	{
		mpc_add(m_x,m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	overload for the -= operator
	HpWrap& operator-= (const HpWrap& rhs)
	{
		mpc_sub(m_x,m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	overload for the *= operator
	HpWrap& operator*= (const HpWrap& rhs)
	{
		mpc_mul(m_x,m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	overload for the /= operator
	HpWrap& operator/= (const HpWrap& rhs)
	{
		mpc_div(m_x,m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	overload for the = operator
	HpWrap& operator= (const HpWrap& rhs)
	{
		mpc_set(m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	defined an abs function
	HpWrap& AbsOf (const HpWrap& rhs)
	{
		mpfr_t tmp;
		mpfr_init2(tmp,P);
		mpc_abs(tmp,rhs.m_x,MPFR_RNDN);

		mpc_set_fr(m_x,tmp,MPC_RNDNN);

		mpfr_clear(tmp);
		return *this;
	}

	//	defined a log function
	HpWrap& LogOf (const HpWrap& rhs)
	{
		mpc_log(m_x,rhs.m_x,MPC_RNDNN);
		return *this;
	}

	//	overload the < operator
	bool operator< (const HpWrap& rhs) const
	{
		return (mpc_cmp(m_x,rhs.m_x)<0);
	}

	//	overload function MulBy for integer type
	HpWrap& MulBy (int lhs,const HpWrap& rhs)
	{
		mpc_mul_si(m_x,rhs.m_x,lhs,MPC_RNDNN);
		return *this;
	}

	//	overload function MulBy for integer type
	HpWrap& MulBy (long int lhs,const HpWrap& rhs)
	{
		mpc_mul_ui(m_x,rhs.m_x,lhs,MPC_RNDNN);
		return *this;
	}

	//	overload function MulBy for mpfr_t* type
	HpWrap& MulBy (mpfr_t* lhs,const HpWrap& rhs)
	{
		mpc_mul_fr(m_x,rhs.m_x,*lhs,MPC_RNDNN);
		return *this;
	}

	//	divide by an integer
	HpWrap& DivBy (int lhs,const HpWrap& rhs)
	{
		mpc_div_ui(m_x,rhs.m_x,lhs,MPC_RNDNN);
		return *this;
	}

	//	set the internal mpc_t variable equal to a complex double
	void Set(dcmplx in)
	{
		mpc_set_d_d(m_x,real(in),imag(in),MPC_RNDNN);
	}

	//	return a complex double variable
	dcmplx Get()  const
	{
		mpfr_t rePart,imPart;
		dcmplx returnVal;

		mpfr_init2(rePart,P);
		mpfr_init2(imPart,P);

		mpc_real(rePart,m_x,MPFR_RNDN);
		mpc_imag(imPart,m_x,MPFR_RNDN);

		returnVal=dcmplx(mpfr_get_ld(rePart,MPFR_RNDN),mpfr_get_ld(imPart,MPFR_RNDN));

		mpfr_clear(rePart);
		mpfr_clear(imPart);

		return returnVal;
	}

	//	return the value of the internal variable

	mpc_t* GetReference() const {return &m_x;}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Specialization of the cpWrap template class to mpfr_t variables
//!	
////////////////////////////////////////////////////////////////////////////////

template <int P> class HpWrap<mpfr_t,P>
{
	private:

	//	internal mpfr_t variable
	mpfr_t m_x;

	public:

	//	constructor initializes with precision P
	HpWrap(){mpfr_init2(m_x,P);}

	//	destructor clears variable to prevent all memory leaks
	~HpWrap(){mpfr_clear(m_x);}

	//	overload for the += operator
	HpWrap& operator+= (const HpWrap& rhs)
	{
		mpfr_add(m_x,m_x,rhs.m_x,MPFR_RNDN);
		return *this;
	}

	//	overload for the -= operator
	HpWrap& operator-= (const HpWrap& rhs)
	{
		mpfr_sub(m_x,m_x,rhs.m_x,MPFR_RNDN);
		return *this;
	}

	//	overload for the *= operator
	HpWrap& operator*= (const HpWrap& rhs)
	{
		mpfr_mul(m_x,m_x,rhs.m_x,MPFR_RNDN);
		return *this;
	}

	//	overload for the /= operator
	HpWrap& operator/= (const HpWrap& rhs)
	{
		mpfr_div(m_x,m_x,rhs.m_x,MPFR_RNDN);
		return *this;
	}

	//	overload for the = operator
	HpWrap& operator= (const HpWrap& rhs)
	{
		mpfr_set(m_x,rhs.m_x,MPFR_RNDN);
		return *this;
	}

	//	multiplication by double type
	HpWrap& MulBy (double lhs,const HpWrap& rhs)
	{
		mpfr_mul_d(m_x,rhs.m_x,lhs,MPFR_RNDN);
		return *this;
	}

	//	set the internal mpfr_t variable equal to real part of a complex double
	void Set(double in)
	{
		mpfr_set_d(m_x,in,MPFR_RNDN);
	}

	//	return a  double variable containing the real part of mpfr_t m_x
	double Get() const
	{
		return mpfr_get_ld(m_x,MPFR_RNDN);
	}

	//	return the value of the internal variable

	mpfr_t* GetReference(){return &m_x;}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Specialization of the cpWrap template class to mpf_t variables
//!	
////////////////////////////////////////////////////////////////////////////////

template <int P> class HpWrap<mpf_t,P>
{
	private:

	//	internal mpf_t variable
	mpf_t m_x;

	public:

	//	constructor initializes with precision P
	HpWrap(){mpf_init2(m_x,P);}

	//	destructor clears variable to prevent all memory leaks
	~HpWrap(){mpf_clear(m_x);}

	//	overload for the += operator
	HpWrap& operator+= (const HpWrap& rhs)
	{
		mpf_add(m_x,m_x,rhs.m_x);
		return *this;
	}

	//	overload for the -= operator
	HpWrap& operator-= (const HpWrap& rhs)
	{
		mpf_sub(m_x,m_x,rhs.m_x);
		return *this;
	}

	//	overload for the *= operator
	HpWrap& operator*= (const HpWrap& rhs)
	{
		mpf_mul(m_x,m_x,rhs.m_x);
		return *this;
	}

	//	overload for the /= operator
	HpWrap& operator/= (const HpWrap& rhs)
	{
		mpf_div(m_x,m_x,rhs.m_x);
		return *this;
	}

	//	overload for the = operator
	HpWrap& operator= (const HpWrap& rhs)
	{

		mpf_set(m_x,rhs.m_x);
		return *this;
	}

	//	multiplication by double type
	HpWrap& MulBy (double lhs,const HpWrap& rhs)
	{
		mpf_t tmp;
		mpf_init2(tmp,P);
		mpf_set_d(tmp,lhs);
		mpf_mul(m_x,rhs.m_x,tmp);
		mpf_clear(tmp);
		return *this;
	}

	//	overload function MulBy for integer type
	HpWrap& MulBy (long int lhs,const HpWrap& rhs)
	{
		mpf_mul_ui(m_x,rhs.m_x,lhs);

		return *this;
	}

	//	set the internal mpfr_t variable equal to real part of a complex double
	void Set(double in)
	{
		mpf_set_d(m_x,in);
	}
	
	//	return a  double variable containing the real part of mpfr_t m_x
	double Get() const
	{
		return mpf_get_d(m_x);
	}

	//	return the value of the internal variable

	mpf_t* GetReference() {return &m_x;}

};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Specialization of the cpWrap template class to complex<mpf_t> variables
//!	
//!	NOTE this is not a proper type, but it acts as a place holder to describe complex mpf variables
//!	The type is implemented by an array containing two mpf_t types
//!	Note that mpf_t typically a lot faster at arithmetic than mpfr_t or mpc_t,
//!	However the gmp library does not contain special functions e.g log,
//!	and it's rounding method does not confirm to IEEE.
//////////////////////////////////////////////////////////////////////////////////

template <int P> class HpWrap<mpfCmplx,P>
{
	private:

	//	internal data type variable
	mpf_t m_x[2];
	//	m_x[0] is the real part and m_x[1] is the imaginary part

	public:

	//	constructor initializes with precision P
	HpWrap(){mpf_init2(m_x[0],P);mpf_init2(m_x[1],P);}

	//	destructor clears variable to prevent all memory leaks
	~HpWrap(){mpf_clear(m_x[0]);mpf_clear(m_x[1]);}

	//	overload for the += operator
	HpWrap& operator+= (const HpWrap& rhs)
	{
		mpf_add(m_x[0],m_x[0],rhs.m_x[0]);
		mpf_add(m_x[1],m_x[1],rhs.m_x[1]);
		return *this;
	}

	//	overload for the -= operator
	HpWrap& operator-= (const HpWrap& rhs)
	{
		mpf_sub(m_x[0],m_x[0],rhs.m_x[0]);
		mpf_sub(m_x[1],m_x[1],rhs.m_x[1]);

		return *this;
	}

	//	overload for the *= operator
	HpWrap& operator*= (const HpWrap& rhs)
	{
		mpf_t tmp,tmp1;
		mpf_init2(tmp,P);
		mpf_init2(tmp1,P);

		mpf_mul(tmp,m_x[0],rhs.m_x[0]);
		mpf_mul(tmp1,m_x[1],rhs.m_x[1]);
		mpf_sub(tmp,tmp,tmp1);

		mpf_mul(m_x[1],m_x[1],rhs.m_x[0]);
		mpf_mul(tmp1,m_x[0],rhs.m_x[1]);
		mpf_add(m_x[1],m_x[1],tmp1);
		
		mpf_set(m_x[0],tmp);

		mpf_clear(tmp);
		mpf_clear(tmp1);

		return *this;
	}

	//	overload for the /= operator
	HpWrap& operator/= (const HpWrap& rhs)
	{
		mpf_t tmp,tmp1,tmp2;
		mpf_init2(tmp,P);
		mpf_init2(tmp1,P);
		mpf_init2(tmp2,P);
		
		mpf_mul(tmp,rhs.m_x[0],rhs.m_x[0]);
		mpf_mul(tmp1,rhs.m_x[1],rhs.m_x[1]);
		mpf_add(tmp2,tmp,tmp1);

		mpf_mul(tmp,m_x[0],rhs.m_x[0]);
		mpf_mul(tmp1,m_x[1],rhs.m_x[1]);
		mpf_add(tmp,tmp,tmp1);

		mpf_mul(m_x[1],m_x[1],rhs.m_x[0]);
		mpf_mul(tmp1,m_x[0],rhs.m_x[1]);
		mpf_sub(m_x[1],m_x[1],tmp1);
		
		mpf_set(m_x[0],tmp);

		mpf_div(m_x[0],m_x[0],tmp2);
		mpf_div(m_x[1],m_x[1],tmp2);

		mpf_clear(tmp);
		mpf_clear(tmp1);
		mpf_clear(tmp2);

		return *this;
	}

	//	overload for the = operator
	HpWrap& operator= (const HpWrap& rhs)
	{
		mpf_set(m_x[0],rhs.m_x[0]);
		mpf_set(m_x[1],rhs.m_x[1]);
		return *this;
	}

	//	defined an abs function
	HpWrap& AbsOf (const HpWrap& rhs)
	{
		mpf_t tmp,tmp2;
		mpf_init2(tmp,P);
		mpf_init2(tmp2,P);

		mpf_mul(tmp,rhs.m_x[1],rhs.m_x[1]);

		mpf_mul(tmp2,rhs.m_x[0],rhs.m_x[0]);

		mpf_add(tmp,tmp,tmp2);

		mpf_sqrt(m_x[0],tmp);

		mpf_set_d(m_x[1],0.0);

		mpf_clear(tmp);
		mpf_clear(tmp2);

		return *this;
	}

	//	define a log function
	HpWrap& LogOf (const HpWrap& rhs)
	{
		mpc_t tmp;
		mpfr_t rePart,imPart;

		mpfr_init2(rePart,P);
		mpfr_init2(imPart,P);
		mpc_init2(tmp,P);

		// convert to mpc type
		mpc_set_f_f(tmp,m_x[0],m_x[1],MPC_RNDNN);

		//	evaluate log
		mpc_log(tmp,tmp,MPC_RNDNN);

		mpc_real(rePart,tmp,MPFR_RNDN);
		mpc_imag(imPart,tmp,MPFR_RNDN);

		//	conert back to mpf type
		mpfr_get_f(m_x[0],rePart,MPFR_RNDN);
		mpfr_get_f(m_x[1],imPart,MPFR_RNDN);

		mpfr_clear(rePart);
		mpfr_clear(imPart);
		mpc_clear(tmp);

		return *this;
	}

	//	overload the < operator
	bool operator< (const HpWrap& rhs)  const
	{
		return (mpf_cmp(m_x[0],rhs.m_x[0])<0);
	}

	//	overload function MulBy for integer type
	HpWrap& MulBy (int lhs,const HpWrap& rhs)
	{
		mpf_mul_ui(m_x[0],rhs.m_x[0],(int) abs(lhs));
		mpf_mul_ui(m_x[1],rhs.m_x[1],(int) abs(lhs));
		
		if(lhs<0)
		{
			mpf_neg(m_x[0],m_x[0]);
			mpf_neg(m_x[1],m_x[1]);
		}
		return *this;
	}

	//	overload function MulBy for integer type
	HpWrap& MulBy (long int lhs,const HpWrap& rhs)
	{
		mpf_mul_ui(m_x[0],rhs.m_x[0],lhs);
		mpf_mul_ui(m_x[1],rhs.m_x[1],lhs);

		return *this;
	}

	//	multiplication by double type
	HpWrap& MulBy (double lhs,const HpWrap& rhs)
	{
		mpf_t tmp;
		mpf_init2(tmp,P);
		mpf_set_d(tmp,lhs);
		mpf_mul(m_x[0],rhs.m_x[0],tmp);
		mpf_mul(m_x[1],rhs.m_x[1],tmp);
		mpf_clear(tmp);
		return *this;
	}

	//	overload function MulBy for mpf_t* type
	HpWrap& MulBy (mpf_t* lhs,const HpWrap& rhs)
	{
		mpf_mul(m_x[0],rhs.m_x[0],*lhs);
		mpf_mul(m_x[1],rhs.m_x[1],*lhs);
		return *this;
	}

	//	divide by an integer
	HpWrap& DivBy (int lhs,const HpWrap& rhs)
	{
		mpf_div_ui(m_x[0],rhs.m_x[0],lhs);
		mpf_div_ui(m_x[1],rhs.m_x[1],lhs);
		return *this;
	}

	//	set the internal complex mpf_t variable equal to a complex double
	void Set(dcmplx in)
	{
		mpf_set_d(m_x[0],real(in));
		mpf_set_d(m_x[1],imag(in));
	}

	//	return a complex double variable
	dcmplx Get()  const
	{
		return dcmplx(mpf_get_d(m_x[0]),mpf_get_d(m_x[1]));
	}

	//	return the value of the internal variable

	mpf_t* GetReference() {return &m_x;}

};

}   //  End namespace utilities

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif
