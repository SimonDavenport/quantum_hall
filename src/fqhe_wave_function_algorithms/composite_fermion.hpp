////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This is the header file for the composite fermion wave function class
//!		This file contains template algorithms to evaluate different types of
//!		composite fermion wave function.
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

#ifndef _COMPOSITE_FERMION_HPP_INCLUDED_
#define _COMPOSITE_FERMION_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "fqhe_wave_function.hpp"
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
#include <utilities/wrappers/high_precision_wrapper.hpp>
#include <utilities/mathematics/dense_linear_algebra.hpp>
#include <utilities/mathematics/combinatorics.hpp>
#ifdef _ENABLE_MPI_
#include <utilities/wrappers/mpi_wrapper.hpp>
#endif
#ifdef _DEBUG_
#include <utilities/general/debug.hpp>
#endif

namespace FQHE
{
///////     STATIC CONST DEFINITIONS      //////////////////////////////////////
//!
//!	Define the 4 precision levels for the high precision switching 
//!
enum precisionLevel {_PREC_LEVEL1_=128, _PREC_LEVEL2_=256, _PREC_LEVEL3_=512, _PREC_LEVEL4_=1024};
static const double tol=0.0000000000001;	//!<	Tolerance in wave function value
											//!
static const double condTolL2=200;			//!<	Value of conditionNumber required
											//!  	to switch to use high precision algorithms
											//!     in negative effective field case
static const double extraNormFactor=10;     //!<    Extra factor to divide each monopole Harmonic
                                            //!     by, to avoid overflow errors at low precision
static const int maxSwitch=42;		        //!<	Maximum electron number to use with
											//!<	+ve effective field case before
											//!		switching to high precision

//!
//! Generate a list of composite fermion wave function program options
//!
inline boost::program_options::options_description GetCompositeFermionOptions()
{
    boost::program_options::options_description cfWaveFunctionOpt("CF Wave Function Options");
    cfWaveFunctionOpt.add_options()
    ("llup", boost::program_options::value<int>()->default_value(2),
     "Number of spin-up composite fermion Landau levels\n")
    ("lldown", boost::program_options::value<int>()->default_value(0),
     "Number of spin-down composite fermion Landau levels\n\n"
     "llup should be >=lldown and both signs should be the same\n\n"
     "  +ve values -> Positive effective field CF wave function\n"
     "  -ve values -> Negative effective field CF wave function\n")
    ("altProjection",boost::program_options::value<bool>()->default_value(false),
     "Set the alternative LLL projection for CF spin unpolarized states "
     "i.e write the wave function as polarised_up *polarised_down *(z_up-z_down)^2\n")
    ("un-proj",boost::program_options::value<bool>()->default_value(false),
     "Set to use un-projected CF wave functions");
     return cfWaveFunctionOpt;
}										
											
//!
//!	A data structure to contain a basic set of composite fermion
//! wave function data
//!
struct CompositeFermionData
{
	int halfFluxAttatched;	//!<		Amount of flux attached to CF wave function
							//!
	double effectiveMonopole;//!<		The effective monopole strength for reduced
							//!			magnetic field. 
	int LLup,LLdown;		//!<		Number of filled spin-up and spin-down CF
							//!		    effective Landau Levels
	bool NEF;				//!<		Option to use negative or positive effective field
							//!
	bool spinPolarised;		//!<		Option to specify spin polarised electrons or not
							//!
	bool altProjection;		//!<		Use alternative projection of spinful CF wave funcs
							//!
	bool bsMode;			//!<		Instead of CF wave functions, evaluate Bonderson-
							//! 		Slingerland wave functions
	bool unprojMode;		//!<		Instead of standard CF wave functions, evaluate
							//! 		unprojected versions (IQHE*jastrow factor)
	CompositeFermionData();			
    void InitFromCommandLine(
    #if _ENABLE_MPI_
        boost::program_options::variables_map* vm,
    const utilities::MpiWrapper& mpi);
    #else
        boost::program_options::variables_map* vm);
    #endif
    #if _ENABLE_MPI_
        void MpiSync(const int syncId, const utilities::MpiWrapper& mpi);
    #endif
};

//!
//!	Contains all data that can only needs to be calculated once	
//!	before any subsequent call to the composite fermion wave function algorithm.
//!
struct PreData
{
	double* llOccupationNumbers;//!<	Allow specified occupation of CF LLS
								//!
	double* cfLz;				//!<	Allow specified z-component of angular
								//!		momentum states to be occupied
	int maxBinom;				//!<	Maximum size of binomial coefficient stored	
								//!
	long int *binomList;        //!<	A list of binomial coefficients k Choose j
								//!
	int* nbrEachLL;	            //!<	Number of terms in the wave function for each nbrLLs
								//!
	double effectiveMonopole;	//!<	Monopole strength
};

double MonopoleHarmonicNorm(const double q, const int l, const double m);   

//////////////////////////////////////////////////////////////////////////////
//!	\brief	This template class contains arbitrary precision algorithms
//!	with which to evaluate various kinds of composite fermion wave functions
//!
//!	The template algorithms will be compiled at arbitrary precision level set 
//!	by the template parameter P.
//!
//!	The template parameter U sets the complex data type used in the template
//!
//!	The template parameter V sets the real data type used in the template
//////////////////////////////////////////////////////////////////////////////
template <typename U, typename V, int P>
class ArbitraryPrecisionCfWavefunctions
{
	private:
	friend class utilities::HpWrap<U,P>;	//!<	The high precision complex type 
								            //!		to be implemented
	friend class utilities::HpWrap<V,P>;	//!<	The high precision real type 
								            //!		to be implemented
	////////////////////////////////////////////////////////////////////////////////
	#if _CF_BENCHMARK_MODE_
		int m_noCalls;		
	#endif
	////////////////////////////////////////////////////////////////////////////////
	private:
	//////		Type independent data		////////////////////////////////////////
	CompositeFermionData *m_cfData;		//!<		Data structure to store CF wave 
										//! 		function data
	PreData *m_preData;					//!<		Data structure to store pre-calculated data	
										//!
	//////		Type V data	(real arbitrary precision)		////////////////////////
	utilities::HpWrap<V,P> *m_facList;  //!<	A list of combinatoric and factorial coefficients
								        //!<	associated with each term in the calculation
	//////		Type U data (complex arbitrary precision)		////////////////////
	utilities::HpWrap<U,P> *m_z;		//!<	Array to store z=v/u
								        //!
	utilities::HpWrap<U,P> *m_zPow;		//!<	Array to store powers of z e.g the z^2 etc.
								        //!
	utilities::HpWrap<U,P> *m_esp;	    //!<	Array to store elementary symmetric polynomials
								        //!
	utilities::HpWrap<U,P> *m_uProd;	//!<	An array to store a list of u*u*u....being looped
								        //!		over. Used to calculate the symmetric polynomials
	utilities::HpWrap<U,P> *m_symmPolyArray;//!<	Array of elementary symmetric polynomials in n-1 variables
								        //!		
	utilities::HpWrap<U,P> *m_zPowStore;//!<	Array storing different powers of z_i
								        //!
	utilities::HpWrap<U,P> *m_slaterArray;  //!<	Slater matrix (which we take the determinant of)
								        //!
	utilities::HpWrap<U,P> *m_slaterArrayUp;//!<	Slater matrix for spin-up particles
								        //!
	utilities::HpWrap<U,P> *m_slaterArrayDown;//!<	Slater matrix for spin-down particles
								        //!
	utilities::HpWrap<U,P> *m_fArray;	//!<	Array to store the values of the function f_j(a,b)
								        //!
	utilities::HpWrap<U,P> *m_pArray;	//!<	Array to store the values of the function P_j(a,b)
								        //!
	//////      FUNCTION DECLARATIONS		////////////////////////////////////////
	public:

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Constructor for ArbitraryPrecisionCfWavefunctions template class.
	//!
	//!	Used to assign memory required for internal calculations of composite 
	//!	fermion wave functions and to pre-calculate ratios of factorials and binomials.
	//!
	//!	Template U corresponds to a complex type. Template V corresponds to a real type.
	//!	Template P corresponds to a precision (in bits).
	////////////////////////////////////////////////////////////////////////////////	
	ArbitraryPrecisionCfWavefunctions(
		WaveFunctionData *wfData,		//!<	The address of an array containing 
										//!		the basic wave function data
		CompositeFermionData *cfData,	//!<	The address of an array containing
										//!		the basic composite fermion data
		PreData *preData)				//!<	The address of an array containing
										//!		pre-calculated data e.g. binomials
	:
	m_cfData(cfData),
	m_preData(preData),
	m_facList(0),
	m_z(0),
	m_zPow(0),
	m_esp(0),
	m_uProd(0),
	m_symmPolyArray(0),
	m_zPowStore(0),
	m_slaterArray(0),
	m_slaterArrayUp(0),
	m_slaterArrayDown(0),
	m_fArray(0),
	m_pArray(0)
	{
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		m_noCalls=0;
		#endif
		////////////////////////////////////////////////////////////////////////////////
		m_z = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr];
		if(m_cfData->spinPolarised || m_cfData->altProjection)
		{
			m_slaterArray = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr*wfData->nbr];
		}
		else
		{
			m_slaterArrayUp   = new utilities::HpWrap<U, P>[wfData->nbrUp*wfData->nbrUp];
			m_slaterArrayDown = new utilities::HpWrap<U, P>[wfData->nbrDown*wfData->nbrDown];
		}
		m_esp = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr];
		if(m_cfData->NEF)
		{
			m_zPow = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr];
			m_zPowStore = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr*wfData->nbr];
			m_uProd = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr];
			m_symmPolyArray = new (std::nothrow) utilities::HpWrap<U, P>[wfData->nbr*wfData->nbr];
		}
		else
		{
			int maxLL = std::max(abs(m_cfData->LLup), abs(m_cfData->LLdown));
			m_fArray = new (std::nothrow) utilities::HpWrap<U, P>[maxLL*maxLL*wfData->nbr];
			m_pArray = new (std::nothrow) utilities::HpWrap<U, P>[maxLL*maxLL*wfData->nbr];
		}
		if(m_cfData->NEF)
		{
			//	To save processing time in the CF wave function algorithm 
			//	(for the negative effective field case),
			//	store a list of the combinatoric factor associated with each term
			long int totalNumber=0;
			totalNumber = m_preData->nbrEachLL[abs(m_cfData->LLup)-1];
			if(m_cfData->LLdown>0)
			{
				totalNumber += m_preData->nbrEachLL[abs(m_cfData->LLdown)-1];
			}
			utilities::cout.DebuggingInfo() << "\tLength of m_facList allocated: " << totalNumber << std::endl;
			m_facList = new utilities::HpWrap<V, P>[totalNumber];
			//	pre-calculate values of all combinatoric factors
			//	spin-up LLs
			this->CalcualtePreFactor(wfData->nbr, wfData->flux, m_cfData->LLup, m_facList);
			if(m_cfData->LLdown>0)
			{
				//	spin-down LLs
				totalNumber = m_preData->nbrEachLL[abs(m_cfData->LLup)-1];
				this->CalcualtePreFactor(wfData->nbr, wfData->flux, m_cfData->LLdown, m_facList+totalNumber);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Destructor for ArbitraryPrecisionCfWavefunctions class
	//!
	//! Used to deallocate memory required for internal calculations.
	//! Called automatically when the object is deleted.
	////////////////////////////////////////////////////////////////////////////////
	~ArbitraryPrecisionCfWavefunctions()
	{	
		////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
			utilities::cout.SecondaryOutput() << "Number of calls to this precision: " 
			                                  << m_noCalls << std::endl;
		#endif
		////////////////////////////////////////////////////////////////////////////
		utilities::cout.DebuggingInfo() << "delete m_facList" << std::endl;
		if(m_facList!=0)			delete[] m_facList;
		utilities::cout.DebuggingInfo()<<"delete m_z"<<std::endl;
		if(m_z!=0)				    delete[] m_z;
		utilities::cout.DebuggingInfo() << "delete m_slaterArray" << std::endl;
		if(m_slaterArray!=0)		delete[] m_slaterArray;
		utilities::cout.DebuggingInfo() << "delete m_slaterArrayUp" << std::endl;
		if(m_slaterArrayUp!=0)	    delete[] m_slaterArrayUp;
		utilities::cout.DebuggingInfo() << "delete m_slaterArrayDown" << std::endl;
		if(m_slaterArrayDown!=0)	delete[] m_slaterArrayDown;
		utilities::cout.DebuggingInfo() << "delete m_zPow" << std::endl;
		if(m_zPow!=0)		        delete[] m_zPow;
		utilities::cout.DebuggingInfo() << "delete m_zPowStore" << std::endl;
		if(m_zPowStore!=0)		    delete[] m_zPowStore;
		utilities::cout.DebuggingInfo() << "delete m_uProd" << std::endl;
		if(m_uProd!=0)			    delete[] m_uProd;
		utilities::cout.DebuggingInfo() << "delete m_esp" << std::endl;
		if(m_esp!=0)				delete[] m_esp;
		utilities::cout.DebuggingInfo() << "delete m_symmPolyArray" << std::endl;
		if(m_symmPolyArray!=0)	    delete[] m_symmPolyArray;
		utilities::cout.DebuggingInfo() << "delete m_fArray" << std::endl;
		if(m_fArray!=0)			    delete[] m_fArray;
		utilities::cout.DebuggingInfo() << "delete m_pArray" << std::endl;
		if(m_pArray!=0)			    delete[] m_pArray;
	}

	private:

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Pre calculation of factorial and binomial factors for each term
	//!	occurring the the CF wave function. 
	//!
	//! Calculates the contents of the facList. This calculation only needs to
	//!	be done in the negative effective case, which requires additional
	//!	factorial factors.
	////////////////////////////////////////////////////////////////////////////////
	void CalcualtePreFactor(
		const int n,		        //!<	Number of particles
		const double flux,			//!<	Flux
		const int nbrLLs,			//!<	Number of Landau levels filled
		utilities::HpWrap<V, P>* facList)//!<	An arbitrary precision list of factorial 
									//!		ratios to be populated
	{
        {
		    utilities::HpWrap<V, P> realTmp;
		    utilities::HpWrap<V, P> *p_list = facList;
		    for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
		    {
			    double m = m_preData->cfLz[d];
                //  Calculate normalization of single particle orbitals
                double norm = FQHE::MonopoleHarmonicNorm(m_preData->effectiveMonopole, 0, 
                    m_preData->cfLz[d])*std::pow(FQHE::extraNormFactor, n/2);
			    for(int t=(int)(m_preData->effectiveMonopole+m); t<=(int)(n-1-m_preData->effectiveMonopole+m); ++t, ++p_list)
			    {
				    utilities::DivFactorial<V, P>(n-1-t,n-1-t-m_preData->effectiveMonopole+m, &realTmp);
				    *p_list = realTmp;
				    utilities::DivFactorial<V,P>(t, t-m_preData->effectiveMonopole-m, &realTmp);
				    *p_list *= realTmp;
				    (*p_list).MulBy(m_preData->binomList[(int)(2*m_preData->effectiveMonopole*m_preData->maxBinom+m_preData->effectiveMonopole+m)], (*p_list));
				    (*p_list).MulBy(1.0/norm, *p_list);
				    //	This line implements (-1)^t
					if(t & 1)   (*p_list).MulBy(-1.0, *p_list);
				    ////////////////////////////////////////////////////////////////////////////////
				    #if _DEBUG_
				    utilities::cout.DebuggingInfo() << "p_list = " << (*p_list).Get() << std::endl;
				    #endif
				    ////////////////////////////////////////////////////////////////////////////////
			    }
		    }
		    //	Now we do the same thing for higher Landau levels
            for(int ll=1, degen=0, cumulativeDegen=m_preData->llOccupationNumbers[0]; ll<nbrLLs; ++ll, cumulativeDegen+=degen)
		    {
			    degen = m_preData->llOccupationNumbers[ll];
			    for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
			    {
				    double m = m_preData->cfLz[d];
                    //  Calculate normalization of single particle orbitals
                    double norm = FQHE::MonopoleHarmonicNorm(m_preData->effectiveMonopole, ll, 
                                    m_preData->cfLz[d])*std::pow(FQHE::extraNormFactor, n/2);
				    for(int s=0; s<=ll; ++s)
				    {
					    //	skip if we are multiplying by a 0 binomial coefficient
					    if(m_preData->effectiveMonopole+m+s>=0 && m_preData->effectiveMonopole+m+s<=2*m_preData->effectiveMonopole+ll)
					    {
						    for(int t=(int)(m_preData->effectiveMonopole+m+s); t<=(int)(n-1-m_preData->effectiveMonopole+m-ll+s); ++t)
						    {
							    if(t>=0 && t<n)
							    {
								    utilities::DivFactorial<V, P>(n-1-t,n-1-t-m_preData->effectiveMonopole+m-ll+s, &realTmp);
								    *p_list = realTmp;
								    utilities::DivFactorial<V, P>(t, t-m_preData->effectiveMonopole-m-s, &realTmp);
								    (*p_list) *= realTmp;
								    (*p_list).MulBy(m_preData->binomList[(int)((
								        2*m_preData->effectiveMonopole+ll)
								        *m_preData->maxBinom+m_preData->effectiveMonopole+m+s)]
								        *m_preData->binomList[ll*m_preData->maxBinom+s], (*p_list));
								    (*p_list).MulBy(1.0/norm,(*p_list));
								    //	This line implements (-1)^t
					                if(t & 1)   (*p_list).MulBy(-1.0, (*p_list));
					                //	This line implements (-1)^s
					                if(s & 1)   (*p_list).MulBy(-1.0, (*p_list));
								    ////////////////////////////////////////////////////////////////////////////////
								    #if _DEBUG_
								    utilities::cout.DebuggingInfo() << "p_list = " << (*p_list).Get() << std::endl;
								    #endif
								    ////////////////////////////////////////////////////////////////////////////////
								    ++p_list;
							    }
						    }
					    }
				    }
			    }
		    }
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Performs pre-calculation of elementary symmetric polynomials and writes
	//!	a table of values to symmPolyArray. In the process, also populates zPowStore 
	//!	(a list of powers of v/u) and uProd (a list of products u,u*u etc.).
	////////////////////////////////////////////////////////////////////////////////
	void GenerateSymmetricPolys(
		const int n,		                    //!<	Number of particles
		dcmplx* u,			                    //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v,			                    //!<	Pointer to the current array of current spinor v coordinates
		utilities::HpWrap<U,P> *zPowStore,	    //!<	Pointer to a table to store powers of z = v/u
		utilities::HpWrap<U,P> *symmPolyArray,  //!<	Pointed to a table to store elementary symmetric polynomials
		utilities::HpWrap<U,P> *uProd)		    //!<	Pointer to a table to store values of u, u*u etc. 
	{
		//	Declare local variables	
		utilities::HpWrap<U,P>	tempU[n];	//	Temporary wrapper value for the u co-ordinate
		utilities::HpWrap<U,P>	temp;		//	Temporary value used in internal calculations
		// To save memory access time, it is convenient to set the ratios of v[j]/u[j] to be the variable z[j].
		for(int j=0;j<n;j++)
		{
			m_z[j].Set(v[j]/u[j]);
			tempU[j].Set(u[j]);
		}
		//	Need to multiply all of the Y matrix elements by prod u_i for i!=j and by u_j^n-1-2*m_preData->effectiveMonopole;
		uProd[0].Set(1.0);
		for(int j=0; j<n; ++j)
		{
			uProd[0] *= tempU[j];
		}
		for(int j=1; j<n; ++j)
		{
			uProd[j] = uProd[0];
		}
		for(int j=0; j<n; ++j)
		{
			for(int k=0; k<n-2-2*m_preData->effectiveMonopole; ++k)
			{
				uProd[j] *= tempU[j];
			}
			//	divide by factors of u if n-2-2*m_preData->effectiveMonopole<0
			for(int k=n-2-2*m_preData->effectiveMonopole; k<0; ++k)
			{
				uProd[j] /= tempU[j];
			}
		}
		//	Pre-calculate and store powers of v[j]/u[j] for later use
		for(int j=0; j<n; ++j)
		{
			zPowStore[j*n].Set(1.0);
		}
		for(int j=0; j<n; ++j)
		{
			m_zPow[j] = m_z[j];
		}
		for(int j=0; j<n; ++j)
		{
			if(j<(int)(n-2*m_preData->effectiveMonopole-1))
			{
				for(int i=0; i<n; ++i)
				{
					zPowStore[i*n+j+1] = m_zPow[i];
				}
			}
			for(int k=0; k<n; ++k)
			{
				m_zPow[k] *= m_z[k];
			}
		}
		//	Recursively evaluate elementary symmetric polynomials in n variables
		//	m_esp[k] contains the kth degree elementary symmetric polynomial
		m_esp[0].Set(1.0);
		m_esp[1]=m_z[0];
		for(int k=2; k<n; ++k)
		{
			m_esp[k].Set(0.0);
		}
		for(int k=2; k<=n; ++k)
		{
			for(int j=std::min((int)k, n-1); j>=1; --j)
			{
				temp = m_z[k-1];
				temp *= m_esp[j-1];
				m_esp[j] += temp;
			}
		}
		//	Finally, use e_k(n-1)=e_k(n) - z_j e_k-1(n-1) to convert to an array symmPolyArray[j][i] of esps in n-1
		//	variables. symmPolyArray[j][i] contains the elementary symm. poly. of order i, in n-1 variables
		//	excluding the jth variable.
		for(int j=0; j<n; ++j)
		{
			symmPolyArray[j*n].Set(1.0);
		}
		for(int j=0; j<n; ++j)
		{
			for(int k=1; k<n; ++k)
			{
				temp = m_z[j];
				temp *= symmPolyArray[j*n+k-1];
				symmPolyArray[j*n+k] = m_esp[k];
				symmPolyArray[j*n+k] -= temp;
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		utilities::cout.DebuggingInfo() << "\n\tSymm poly array:" << std::endl;
		for(int j=0; j<n; ++j)
		{
			for(int k=1; k<n; ++k)
			{
				utilities::cout.DebuggingInfo() << "\t[" << j << "][" << k << "]=" 
				                                << symmPolyArray[j*n+k].Get() << std::endl;
			}
		}
		#endif
		////////////////////////////////////////////////////////////////////////////////
		return;
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Performs pre-calculation of Jorgen's polynomials and writes
	//!	a table of values to symmPolyArray. In the process, also populates zPowStore 
	//!	(a list of powers of v/u).
	//!
	////////////////////////////////////////////////////////////////////////////////
    void GenerateJorgenPolynomials(
        const int n,		                    //!<	Number of particles
		dcmplx* u,			                    //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v,			                    //!<	Pointer to the current array of current spinor v coordinates
		utilities::HpWrap<U,P> *zPowStore,	    //!<	Pointer to a table to store powers of z = v/u
		utilities::HpWrap<U,P> *symmPolyArray,  //!<	Pointed to a table to store elementary symmetric polynomials
		utilities::HpWrap<U,P> *uProd)		    //!<	Pointer to a table to store values of u to some power
    {
		utilities::HpWrap<U,P>	tempU[n];
        utilities::HpWrap<U,P>	tempV[n];
		for(int j=0; j<n; ++j)
		{
			m_z[j].Set(v[j]/u[j]);
			tempU[j].Set(u[j]);
            tempV[j].Set(v[j]);
		}
        //	Pre-calculate and store powers of v[j]/u[j] for later use
		for(int j=0; j<n; ++j)
		{
			zPowStore[j*n].Set(1.0);
		}
		for(int j=0; j<n; ++j)
		{
			m_zPow[j] = m_z[j];
		}
		for(int j=0; j<n; ++j)
		{
			if(j<(int)(n-2*m_preData->effectiveMonopole-1))
			{
				for(int i=0; i<n; ++i)
				{
					zPowStore[i*n+j+1] = m_zPow[i];
				}
			}
			for(int k=0; k<n; ++k)
			{
				m_zPow[k] *= m_z[k];
			}
		}
        //	Need to multiply all of the Y matrix elements by u_j^n-1-2*m_preData->effectiveMonopole;
		for(int j=0; j<n; ++j)
		{
			uProd[j].Set(1.0);
		}
		for(int j=0; j<n; ++j)
		{
			for(int k=0; k<n-1-2*m_preData->effectiveMonopole; ++k)
			{
				uProd[j] *= tempU[j];
			}
			//	divide by factors of u if n-1-2*m_preData->effectiveMonopole<0
			for(int k=n-1-2*m_preData->effectiveMonopole; k<0; ++k)
			{
				uProd[j] /= tempU[j];
			}
		}
        //	Recursively evaluate the Jorgen polynomials
        //  Define f^i_t,0 = { 1: t=0
        //                   { 0: t>=1
        for(int i=0; i<n; ++i)
		{   
            symmPolyArray[i*n].Set(1.0);
            for(int t=1; t<n; ++t)
            {
                symmPolyArray[i*n+t].Set(0.0);
            }
		}
        //  Now go through from j=1 to N and set f^i_t,j = {  f^i_{t,j-1}         j=i
        //                                                 {  f^i_{t,j-1}u_j      j!=i and t=0
        //                                                 {  f^i_{t,j-1}u_j + f^i_{t-1,j-1} v_n  j!=i and t>0
        utilities::HpWrap<U,P>	temp;
        utilities::HpWrap<U,P>	temp2;
        for(int j=0; j<n; ++j)
        {
            for(int i=0; i<n; ++i)
            {
                if(i!=j)
                {   
                    //  i!=j and t>0 case
                    for(int t=n-1; t>0; --t)
                    {
                        temp = symmPolyArray[i*n+t];
                        temp *= tempU[j];
                        temp2 = symmPolyArray[i*n+t-1];
                        temp2 *= tempV[j];
                        temp2 += temp;
                        symmPolyArray[i*n+t] = temp2;
                    }
                    //  i!=j and t=0 case
                    symmPolyArray[i*n] *= tempU[j];
                }
            }
        }        
        ////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		utilities::cout.DebuggingInfo() << "\n\tJorgen poly array:" << std::endl;
		for(int j=0; j<n; ++j)
		{
			for(int k=0; k<n; ++k)
			{
				utilities::cout.DebuggingInfo() << "\t[" << j << "][" << k << "]=" 
				                                << symmPolyArray[j*n+k].Get() << std::endl;
			}
		}
		#endif
		////////////////////////////////////////////////////////////////////////////////
    }
    
	////////////////////////////////////////////////////////////////////////////////
	//!	\brief Pre calculation of the "p-functions" pArray and also z=v/u.
	//!
	//!	See the relevant appendix of J.K. Jain's book on composite fermions
	//!	for the definition of these functions.
	////////////////////////////////////////////////////////////////////////////////
	void GeneratePfunctions(
		const int n,	                //!<	Number of particles
		const int nbrLLs,		        //!<	Number of Landau levels
		dcmplx* u,		                //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v,		                //!<	Pointer to the current array of current spinor v coordinates
		utilities::HpWrap<U,P> *pArray,	//!<	Pointer to a table of "p-functions" to be populated
		utilities::HpWrap<U,P> *z)		//!<	Pointer to an array containing v/u to be populated
	{
		utilities::HpWrap<U,P> tmp1,tmp2,tmp3; 
		utilities::HpWrap<U,P> runTot,runTot2;
		utilities::HpWrap<U,P> *p_fArray;
		utilities::HpWrap<U,P> *p_pArray;
		for(int j=0; j<n; ++j)
		{
			z[j].Set(v[j]/u[j]);
		}
		//	First, determine f_j(alpha, beta) functions
		for(int a=0; a<nbrLLs; ++a)
		{
			for(int b=0; b<nbrLLs; ++b)
			{
				if((a==0 && b==0) || (a+b>=nbrLLs))
				{}
				else
				{
					for(int j=0; j<n; ++j)
					{
						runTot.Set(0.0);
						//	Store z_j
						tmp1=*(z+j);
						for(int k=0;k<n;k++)
						{
							if(k!=j)
							{
								//	Store z_k
								tmp2 = *(z+k);
								//	Store z_k - z_j
								tmp3 = tmp2;
								tmp3 -= tmp1;
								runTot2.Set(1.0);
								for(int i=0; i<(a+b); ++i)
								{
									if(i<a)
									{
										//	calculate z_k^a
										runTot2 *= tmp2;
									}
									//	Calculate 1/(z_k-z_j)^(a+b)
									runTot2 /= tmp3;
								}
								//	Sum over k
								runTot += runTot2;
							}
						}
						//	Implement (-1)^b
						if(b&1)
						{
							runTot.MulBy(-1,runTot);
						}
						tmp1.Set(*(u+j));
						//	Divide by u_j^(a+b)
						for(int i=0; i<(a+b); ++i)
						{
							runTot /= tmp1;
						}
						//	Set pointer to the correct row of m_fArray
						p_fArray = m_fArray+(a*nbrLLs+b)*n;
						//	Store array of values
						*(p_fArray+j) = runTot;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//output matrix elements
		for(int a=0; a<nbrLLs; ++a){
			for(int b=0; b<nbrLLs; ++b){
				if((a==0 && b==0) ||(a+b>=nbrLLs))
				{}
				else
				{
					for(int j=0; j<n; ++j){
						utilities::cout.DebuggingInfo() << "m_fArray[" << a << "][" << b << "][" 
						    << j << "]=" << (*(m_fArray+(a*nbrLLs+b)*n+j)).Get()<<std::endl;
					}}}} getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////

		//	Second, determine P_j(alpha, beta) functions. These are given e.g in J.K. Jain 
		//  "Composite Fermions" appendix J.2.2. pp. 496--497
		//	Implement P_j(0,0)=1;
		p_pArray=pArray;
		for(int j=0; j<n; ++j)
		{
			(*(p_pArray+j)).Set(1.0);
		}
		//	Implement P_j(1,0)=m_cfData->halfFluxAttatched*f_j(1,0);
		p_pArray=pArray+1*nbrLLs*n;
		p_fArray=m_fArray+1*nbrLLs*n;
		for(int j=0; j<n; ++j)
		{
			*(p_pArray+j) = *(p_fArray+j);
			(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
		}
		//	Implement P_j(0,1)=m_cfData->halfFluxAttatched*f_j(0,1);
		p_pArray = pArray+1*n;
		p_fArray = m_fArray+1*n;
		for(int j=0; j<n; ++j)
		{
			*(p_pArray+j) = *(p_fArray+j);
			(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
		}
		if(nbrLLs>2)
		{
			//	Implement P_j(2,0)=m_cfData->halfFluxAttatched^2*f_j(1,0)^2-m_cfData->halfFluxAttatched*f(2,0);
			p_pArray = pArray+2*nbrLLs*n;
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				(*(p_pArray+j)) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				(*(p_pArray+j)) -= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			//	Implement P_j(0,2)=m_cfData->halfFluxAttatched^2*f_j(0,1)^2-m_cfData->halfFluxAttatched*f(0,2);
			p_pArray = pArray+2*n;
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				(*(p_pArray+j)) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) -= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			//	Implement P_j(1,1)=m_cfData->halfFluxAttatched^2*f_j(0,1)*f_j(1,0)
			//  -m_cfData->halfFluxAttatched*f(1,1);
			p_pArray = pArray+1*nbrLLs*n+1*n;
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
			}
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+1*nbrLLs*n+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) -= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
		}
		if(nbrLLs>3)
		{
			//	Implement P_j(3,0)=m_cfData->halfFluxAttatched*f_j(1,0)*[m_cfData->halfFluxAttatched^2*f_j(1,0)^2-3*m_cfData->halfFluxAttatched*f_j(2,0)]+2*m_cfData->halfFluxAttatched*f_j(3,0);
			p_pArray = pArray+3*nbrLLs*n;
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				*(p_pArray+j) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(3, *(p_fArray+j));
				*(p_pArray+j) -= tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) *= *(p_fArray+j);
			}
			p_fArray = m_fArray+3*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(2, *(p_fArray+j));
				*(p_pArray+j) += tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			//	Implement P_j(0,3)=m_cfData->halfFluxAttatched*f_j(0,1)*[m_cfData->halfFluxAttatched^2*f_j(0,1)^2
			//  -3*m_cfData->halfFluxAttatched*f_j(0,2)]+2*m_cfData->halfFluxAttatched*f_j(0,3);
			p_pArray = pArray+3*n;
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				*(p_pArray+j) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(3, *(p_fArray+j));
				*(p_pArray+j) -= tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) *= *(p_fArray+j);
			}
			p_fArray = m_fArray+3*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(2, *(p_fArray+j));
				*(p_pArray+j) += tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			//	Implement P_j(2,1)=m_cfData->halfFluxAttatched*f_j(0,1)[m_cfData->halfFluxAttatched^2
			//  *f_j(1,0)^2-m_cfData->halfFluxAttatched*f_j(2,0)]-2
			//  *m_cfData->halfFluxAttatched^2*f_j(1,0)*f_j(1,1)+2*m_cfData->halfFluxAttatched*f_(2,1);
			p_pArray = pArray+2*nbrLLs*n+1*n;
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				*(p_pArray+j) *= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) -= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) *= *(p_fArray+j);
			}
			for(int j=0; j<n; ++j)
			{
				p_fArray = m_fArray+1*nbrLLs*n;
				tmp1.MulBy(2*m_cfData->halfFluxAttatched, *(p_fArray+j));
				p_fArray = m_fArray+1*nbrLLs*n+1*n;
				tmp1 *= *(p_fArray+j);
				*(p_pArray+j) -= tmp1;
			}
			p_fArray = m_fArray+2*nbrLLs*n+1*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(2, *(p_fArray+j));
				*(p_pArray+j) += tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			//	Implement P_j(1,2)=f_j(1,0)[f_j(0,1)^2-f_j(0,2)]-2*f_j(0,1)*f_j(1,1)+2*f_(1,2);
			p_pArray = pArray+1*nbrLLs*n+2*n;
			p_fArray = m_fArray+1*n;
			for(int j=0; j<n; ++j)
			{
				*(p_pArray+j) = *(p_fArray+j);
				(*(p_pArray+j)) *= (*(p_fArray+j));
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+2*n;
			for(int j=0; j<n; ++j)
			{
				(*(p_pArray+j)) -= *(p_fArray+j);
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
			p_fArray = m_fArray+1*nbrLLs*n;
			for(int j=0; j<n; ++j)
			{
				(*(p_pArray+j)) *= *(p_fArray+j);
			}
			for(int j=0; j<n; ++j)
			{
				p_fArray = m_fArray+1*n;
				tmp1.MulBy(2*m_cfData->halfFluxAttatched, (*(p_fArray+j)));
				p_fArray = m_fArray+1*nbrLLs*n+1*n;
				tmp1 *= *(p_fArray+j);
				*(p_pArray+j) -= tmp1;
			}
			p_fArray = m_fArray+1*nbrLLs*n+2*n;
			for(int j=0; j<n; ++j)
			{
				tmp1.MulBy(2, *(p_fArray+j));
				*(p_pArray+j) += tmp1;
				(*(p_pArray+j)).MulBy(m_cfData->halfFluxAttatched, (*(p_pArray+j)));
			}
		}
		if(nbrLLs>4)
		{
			std::cerr << "\n\tERROR: cannot do nbrLLs>4 at the moment!" << std::endl;
			exit(EXIT_FAILURE);
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//  Print matrix elements
		for(int a=0; a<nbrLLs; ++a){
		for(int b=0; b<nbrLLs; ++b){
		if(a+b<nbrLLs)
		{
		for(int j=0; j<n; ++j){
		utilities::cout.DebuggingInfo() << "pArray[" << a << "][" << b << "][" << j << "]=" 
		                                << (*(pArray+(a*nbrLLs+b)*n+j)).Get() << std::endl;
		}}}}
		getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////
		return;
	}

	public:

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	Calculates a condition number for the nth elementary symmetric 
	//! polynomial evaluation. This will act as an effective condition number
	//! for the full CF wave function calculation (negative effective field)
	//! Given by: \f[ k e_k (|z|)/|e_k| \f]
	//!
	//!	\return		The value of the condition number.
	////////////////////////////////////////////////////////////////////////////////
	double CalculateConditionNumber(
		const int n,	        //!<	Number of particles
		dcmplx* u,		        //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v)		        //!<	Pointer to the current array of current spinor v coordinates	
	{
	    //  Determine the closest particle to the north pole
	    double minU = 100000;
	    for(int i=0; i<n; ++i)
        {
            double absU = abs(u[i]);
            if(absU<minU)   minU = absU;
        }
        ////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		utilities::cout.DebuggingInfo() << "\tCondition Number: " << 1.0/minU << std::endl;
        #endif
		////////////////////////////////////////////////////////////////////////////////
		return 1.0/minU;
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	Composite fermion wave function algorithm (negative effective field) 
	//!	for the sphere	geometry and spin/valley polarized particles.	
	//!
	//! This function returns the log of the mod squared of the CF wave function
	//! for a given set of co-ordinates.
	//!	
	//! There are two steps: first we populate a n x n Slater matrix with 
	//! lowest Landau level projected monopole harmonics, then we calculate
	//! the log of its Slater determinant.
	//!	
	//!	\return		log of wave function (NOT NORMALIZED)
	////////////////////////////////////////////////////////////////////////////////
	dcmplx PolarisedNEF(
		const int n,	        //!<	Number of particles
		const int nbrLLs,		//!<	Number of Landau levels in the construction
		dcmplx* u,		        //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v)		        //!<	Pointer to the current array of current spinor v coordinates
	{
        ////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		++m_noCalls;
		#endif
		////////////////////////////////////////////////////////////////////////////////

		//	Generate elementary symmetric polynomials and prod u_i for i!=j * u_j^n-1-2*m_preData->effectiveMonopole
		//this->GenerateSymmetricPolys(n,u,v,m_zPowStore,m_symmPolyArray,m_uProd);
        this->GenerateJorgenPolynomials(n,u,v,m_zPowStore,m_symmPolyArray,m_uProd);
		//	Populate the Slater matrix up with projected monopole functions
		utilities::HpWrap<U,P> term;
		utilities::HighPrecKahanAccumulation<utilities::HpWrap<U,P> > runTot;
		//  (see highPrecisionWrapper.h for definition of HighPrecKahanAccumulation)
		utilities::HpWrap<V,P>* p_list = 0;
		//	First we work out states filling the lowest effective Landau level.
		for(int j=0; j<n; ++j)
		{
			p_list = m_facList;
			for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
			{
				double m = m_preData->cfLz[d];
				runTot.Init();
				for(int t=(int)(m_preData->effectiveMonopole+m); t<=(int)(n-1-m_preData->effectiveMonopole+m);
				++t, ++p_list)
				{
					term.MulBy(p_list->GetReference(), m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
					term *= m_symmPolyArray[j*n+t];
					runTot += term;
				}
				runTot *= m_uProd[j];
				m_slaterArray[d*n+j] = runTot.GetTotal();
			}
		}
		//	Now we do the same thing for higher Landau levels
        for(int ll=1, degen=0, cumulativeDegen=m_preData->llOccupationNumbers[0]; ll<nbrLLs; ++ll, cumulativeDegen+=degen)
		{
            degen = m_preData->llOccupationNumbers[ll];
            utilities::HpWrap<V,P> *p_next = p_list;
			for(int j=0; j<n; ++j)
			{
				p_list = p_next;
				for(int d=(int)cumulativeDegen; d<(int)(cumulativeDegen+degen); ++d)
				{
					double m = m_preData->cfLz[d];
					runTot.Init();
					for(int s=0; s<=ll; ++s)
					{
						//	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+m+s>=0 && m_preData->effectiveMonopole+m+s<=2*m_preData->effectiveMonopole+ll)
						{
							for(int t=(int)(m_preData->effectiveMonopole+m+s); t<=(int)(n-1-m_preData->effectiveMonopole+m-ll+s); ++t)
							{
								if(t>=0 && t<n)
								{
									term.MulBy(p_list->GetReference(), m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
									term *= m_symmPolyArray[j*n+t];
									runTot += term;
									++p_list;
								}
							}
						}
					}
					runTot *= m_uProd[j];
					m_slaterArray[d*n+j] = runTot.GetTotal();
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//	print out matrix elements
		for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
		utilities::cout.DebuggingInfo() << "Y[" << i << "][" << j << "]=" 
		                                << m_slaterArray[i*n+j].Get() << std::endl;
		}}
		getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////
		//	The final step is to calculate the determinant
		return utilities::linearAlgebra::LogDeterminant<U, P>(m_slaterArray,n);
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	Composite fermion wave function algorithm (negative effective field) 
	//!	for the sphere geometry and non spin/valley polarized particles.
	//!
	//! This function returns the log of the mod squared of the CF wave function
	//! for a given set of co-ordinates.
	//! There are two steps: first we populate a nbrUp x nbrUp Slater matrix and
	//! a nbrDown x nbrDown Slater matrix with lowest Landau level projected 
	//! monopole harmonics, then we calculate the sum of logs of the Slater 
	//! determinants.
	//!	
	//!	\return	log of wave function (NOT NORMALIZED)
	////////////////////////////////////////////////////////////////////////////////
	dcmplx NonPolarisedNEF(
		const int nbrUp,		//!<	Number of spin/valley up particles
		const int nbrDown,		//!<	Number of spin/valley down particles
		const int LLup,		    //!<	Number of spin/valley up Landau levels
		const int LLdown,	    //!<	Number of spin/valley down Landau levels
		dcmplx* u,	            //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v)	            //!<	Pointer to the current array of current spinor v coordinates
	{
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		++m_noCalls;
		#endif
		////////////////////////////////////////////////////////////////////////////////
		int n = nbrUp+nbrDown;
		//	Generate elementary symmetric polynomials and prod u_i for i!=j * u_j^n-1-2*m_preData->effectiveMonopole
		//this->GenerateSymmetricPolys(n,u,v,m_zPowStore,m_symmPolyArray,m_uProd);
		this->GenerateJorgenPolynomials(n,u,v,m_zPowStore,m_symmPolyArray,m_uProd);
		//	Populate the Slater matrix up with projected monopole functions
        utilities::HpWrap<U,P> term;
        utilities::HighPrecKahanAccumulation<utilities::HpWrap<U,P> > runTot;
        //  (see highPrecisionWrapper.h for definition of HighPrecKahanAccumulation)
		utilities::HpWrap<V,P>* p_list=0;
		for(int j=0; j<nbrUp; ++j)
		{
			p_list = m_facList;
			for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
			{
				double m = m_preData->cfLz[d];
				runTot.Init();
				for(int t=(int)(m_preData->effectiveMonopole+m); t<=(int)(n-1-m_preData->effectiveMonopole+m); ++t, ++p_list)
				{
					term.MulBy(p_list->GetReference(), m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
					term *= m_symmPolyArray[j*n+t];
					runTot += term;
				}
				runTot *= m_uProd[j];
				m_slaterArrayUp[d*nbrUp+j] = runTot.GetTotal();
			}
		}
		// Now we do the same thing for higher spin up Landau levels
		for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<LLup; ++ll, cumulativeDegen+=degen)
		{
            degen = m_preData->llOccupationNumbers[ll];
            utilities::HpWrap<V, P> *p_next = p_list;
			for(int j=0; j<nbrUp; ++j)
			{
				p_list = p_next;
				for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
				{
					double m = m_preData->cfLz[d];
					runTot.Init();
					for(int s=0; s<=ll; ++s)
					{
						//	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+m+s>=0 && m_preData->effectiveMonopole+m+s<=2*m_preData->effectiveMonopole+ll)
						{
							for(int t=(int)(m_preData->effectiveMonopole+m+s); t<=(int)(n-1-m_preData->effectiveMonopole+m-ll+s); ++t)
							{
								if(t>=0 && t<n)
								{
									term.MulBy(p_list->GetReference(), m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
									term *= m_symmPolyArray[j*n+t];
									runTot += term;
									++p_list;
								}
							}
						}
					}
					runTot *= m_uProd[j];
					m_slaterArrayUp[d*nbrUp+j] = runTot.GetTotal();
				}
			}
		}
		//	Now follow the same procedure with the spin down Landau Levels
		for(int j=nbrUp; j<n; ++j)
		{
			p_list = m_facList;
			for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
			{
				double m = m_preData->cfLz[d];
				runTot.Init();
				for(int t=(int)(m_preData->effectiveMonopole+m); t<=(int)(n-1-m_preData->effectiveMonopole+m); ++t, ++p_list)
				{
					term.MulBy(p_list->GetReference(), m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
					term *= m_symmPolyArray[j*n+t];
					runTot += term;
				}
				runTot *= m_uProd[j];
				m_slaterArrayDown[d*nbrDown+j-nbrUp] = runTot.GetTotal();
			}
		}
		// Now we do the same thing for higher spin up Landau levels
		for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<LLdown; ++ll, cumulativeDegen += degen)
		{
            degen = m_preData->llOccupationNumbers[ll];
            utilities::HpWrap<V, P> *p_next = p_list;
			for(int j=nbrUp; j<n; ++j)
			{
				p_list = p_next;
				for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
				{
					double m = m_preData->cfLz[d];
					runTot.Init();
					for(int s=0; s<=ll; ++s)
					{
					    //	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+m+s>=0 && m_preData->effectiveMonopole+m+s<=2*m_preData->effectiveMonopole+ll)
						{
							for(int t=(int)(m_preData->effectiveMonopole+m+s); t<=(int)(n-1-m_preData->effectiveMonopole+m-ll+s); ++t)
							{
								if(t>=0 && t<n)
								{
									term.MulBy(p_list->GetReference(),m_zPowStore[j*n+(int)(n-1-t-m_preData->effectiveMonopole+m)]);
									term *= m_symmPolyArray[j*n+t];
									runTot += term;
									++p_list;
								}
							}
						}
					}
					runTot *= m_uProd[j];
					m_slaterArrayDown[d*nbrDown+j-nbrUp] = runTot.GetTotal();
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//	print out matrix elements
		for(int i=0; i<nbrUp; ++i){
		for(int j=0; j<nbrUp; ++j){
		utilities::cout.DebuggingInfo() << "Y[" << i << "][" << j << "]=" 
		                                << m_slaterArrayUp[i*nbrUp+j].Get() << std::endl;
		}}
		getchar();
		for(int i=0; i<nbrDown; ++i){
		for(int j=0; j<nbrDown; ++j){
		utilities::cout.DebuggingInfo() << "Y[" << i << "][" << j << "]="
		                                << m_slaterArrayDown[i*nbrDown+j].Get() << std::endl;
		}}
		getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////
		//	The final step is to calculate the determinants
		return utilities::linearAlgebra::LogDeterminant<U, P>(m_slaterArrayUp, nbrUp)      
		    +utilities::linearAlgebra::LogDeterminant<U, P>(m_slaterArrayDown, nbrDown);
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	Composite fermion wave function algorithm (positive effective field) 
	//!	for the sphere geometry and spin/valley polarized particles.	
	//!
	//! This function returns the log of the mod squared of the CF wave function
	//! for a given set of co-ordinates.
	//!	
	//! There are two steps: first we populate a n x n Slater matrix with 
	//! lowest Landau level projected monopole harmonics, then we calculate
	//! the log of its Slater determinant.
	//!	
	//!	\return		log of wave function (NOT NORMALIZED)
	////////////////////////////////////////////////////////////////////////////////
	dcmplx PolarisedPEF(
		const int n,	        //!<	Number of particles
		const int nbrLLs,		//!<	Number of Landau levels in the construction
		dcmplx* u,		        //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v)		        //!<	Pointer to the current array of current spinor v coordinates
	{
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		++m_noCalls;
		#endif
		////////////////////////////////////////////////////////////////////////////////
		if(nbrLLs==1)
		{
			//	Use Laughlin 1/3 wave function
            dcmplx wfValue = 0.0;
			for (int i=0; i<n-1; ++i)
			{
				for (int j=i+1; j<n; ++j)
				{
					wfValue += log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
			wfValue *= (1.0+2.0*m_cfData->halfFluxAttatched);
			return wfValue;
		}
		//	Generate P_j(alpha, beta) functions.
		this->GeneratePfunctions(n, nbrLLs, u, v, m_pArray, m_z);
		//	Populate the Slater matrix with projected monopole functions
		for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
		{
			double m = m_preData->cfLz[d];
			for(int j=0; j<n; ++j)
			{
				(*(m_slaterArray+d*n+j)).Set((double)m_preData->binomList[(int)(2*m_preData->effectiveMonopole*m_preData->maxBinom+m_preData->effectiveMonopole-m)]);
				int imax = m_preData->effectiveMonopole-m;
				for(int i=0; i<imax; ++i)
				{
					*(m_slaterArray+d*n+j) *= *(m_z+j);
				}
			}
		}
        utilities::HpWrap<U,P> term,runTot;
		utilities::HpWrap<U,P> tmpCmplx1,tmpCmplx2;
        for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<nbrLLs; ++ll, cumulativeDegen+=degen)
		{
			degen = m_preData->llOccupationNumbers[ll];
			for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
			{
				double m = m_preData->cfLz[d];
				for(int j=0; j<n; ++j)
				{
					runTot.Set(0.0);
					tmpCmplx1.Set(*(v+j));
					tmpCmplx2.Set(*(u+j));
					for(int s=0; s<=ll; ++s)
					{
						//	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+ll-m-s>=0 && (m_preData->effectiveMonopole+ll-m-s)<=(2*m_preData->effectiveMonopole+ll))
						{
							term.MulBy(m_preData->binomList[ll*m_preData->maxBinom+s]*m_preData->binomList[(int)((2*m_preData->effectiveMonopole+ll)*m_preData->maxBinom+m_preData->effectiveMonopole+ll-m-s)],m_pArray[(s*nbrLLs+ll-s)*n+j]);
							int imax = ll-s;
							for(int i=0; i<imax; ++i)
							{
								term *= tmpCmplx1;
							}
							for(int i=0; i<s; ++i)
							{
								term *= tmpCmplx2;
							}
							//This line implements (-1)^s
							(s & 1)	? (runTot-=term) : (runTot+=term);
						}
					}
					int imax = (int)(m_preData->effectiveMonopole-m);
					for(int i=0; i<imax; ++i)
					{
						runTot *= *(m_z+j);
					}
					for(int i=imax; i<0; ++i)
					{
						runTot /= *(m_z+j);
					}
					*(m_slaterArray+d*n+j) = runTot;
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//	Print out matrix elements
		for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
		utilities::cout.DebuggingInfo() << "Y[" << i << "][" << j << "]=" 
		                                << (*(m_slaterArray+i*n+j)).Get() << std::endl;
		}}
		getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////
		//	The final step is to calculate the determinant
		dcmplx wfValue = utilities::linearAlgebra::LogDeterminant<U,P>(m_slaterArray,n);
		//	Multiply by a Jastrow factor
		for (int i=0; i<n-1; ++i)
		{
			for (int j=i+1; j<n; ++j)
			{
				wfValue += 2.0*m_cfData->halfFluxAttatched*log(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}
		//	Multiply by a Product_j u_j^(2Q) factor, which we removed from the determinant
		for(int i=0; i<n; ++i)
		{
			wfValue += 2.0*m_preData->effectiveMonopole*(log(*(u+i)));
		}
		return wfValue;
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	Composite fermion wave function algorithm (positive effective field) 
	//!	for the sphere geometry and non spin/valley polarized particles.
	//!
	//! This function returns the log of the mod squared of the CF wave function
	//! for a given set of co-ordinates.
	//! There are two steps: first we populate a nbrUp x nbrUp Slater matrix and
	//! a nbrDown x nbrDown Slater matrix with lowest Landau level projected 
	//! monopole harmonics, then we calculate the sum of logs of the Slater 
	//! determinants.
	//!	
	//!	\return	log of wave function (NOT NORMALIZED)
	////////////////////////////////////////////////////////////////////////////////
	dcmplx NonPolarisedPEF(
		const int nbrUp,		//!<	Number of spin/valley up particles
		const int nbrDown,		//!<	Number of spin/valley down particles
		const int LLup,		    //!<	Number of spin/valley up Landau levels
		const int LLdown,	    //!<	Number of spin/valley down Landau levels
		dcmplx* u,	            //!<	Pointer to the current array of current spinor u coordinates
		dcmplx* v)	            //!<	Pointer to the current array of current spinor v coordinates
	{
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		++m_noCalls;
		#endif
		////////////////////////////////////////////////////////////////////////////////
        int n = nbrUp+nbrDown;    //  total number of particles
		if(LLup==1 && LLdown==1)
		{
		    dcmplx wfValue = 0.0;
			//	Use Halperin 2p+1 2p+1 2p wavefunction
			for (int i=0; i<nbrUp-1; ++i)
			{
				for (int j=i+1; j<nbrUp; ++j)
				{
					wfValue += (1.0+2.0*m_cfData->halfFluxAttatched)*log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
			for (int i=nbrUp; i<n-1; ++i)
			{
				for (int j=i+1; j<n; ++j)
				{
					wfValue += (1.0+2.0*m_cfData->halfFluxAttatched)*log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
			for (int i=0; i<nbrUp; ++i)
			{
				for (int j=nbrUp; j<n; ++j)
				{
					wfValue += (2.0*m_cfData->halfFluxAttatched)*log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
			return wfValue;
		}
		//	Generate P_j(alpha, beta) functions.
		this->GeneratePfunctions(n,std::max(LLup,LLdown), u, v, m_pArray, m_z);
		//	Populate the spin-up Slater matrix with projected monopole functions for the spin up particles
		for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
		{
			double m = m_preData->cfLz[d];
			for(int j=0; j<nbrUp; ++j)
			{
				(*(m_slaterArrayUp+d*nbrUp+j)).Set((double)m_preData->binomList[(int)(2*m_preData->effectiveMonopole*m_preData->maxBinom+m_preData->effectiveMonopole-m)]);
				int imax = m_preData->effectiveMonopole-m;
				for(int i=0; i<imax; ++i)
				{
					*(m_slaterArrayUp+d*nbrUp+j) *= *(m_z+j);
				}
			}
		}
        utilities::HpWrap<U,P> term, runTot;
		utilities::HpWrap<U,P> tmpCmplx1, tmpCmplx2;
        for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<LLup; ++ll, cumulativeDegen+=degen)
        {
			degen = m_preData->llOccupationNumbers[ll];
			for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
			{
				double m = m_preData->cfLz[d];
				for(int j=0; j<nbrUp; ++j)
				{
					runTot.Set(0.0);
					tmpCmplx1.Set(*(v+j));
					tmpCmplx2.Set(*(u+j));
					for(int s=0; s<=ll; ++s)
					{
						//	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+ll-m-s>=0 && m_preData->effectiveMonopole+ll-m-s<=2*m_preData->effectiveMonopole+ll)
						{
							term.MulBy(m_preData->binomList[ll*m_preData->maxBinom+s]*m_preData->binomList[(int)((2*m_preData->effectiveMonopole+ll)*m_preData->maxBinom+m_preData->effectiveMonopole+ll-m-s)],*(m_pArray+(s*LLup+ll-s)*n+j));
							int imax = ll-s;
							for(int i=0; i<imax; ++i)
							{
								term *= tmpCmplx1;
							}
							for(int i=0; i<s; ++i)
							{
								term *= tmpCmplx2;
							}
							//This line implements (-1)^s
							(s & 1)	? (runTot-=term) : (runTot+=term);
						}
					}
					int imax = (int)(m_preData->effectiveMonopole-m);
					for(int i=0; i<imax; ++i)
					{
						runTot *= *(m_z+j);
					}
					for(int i=imax; i<0; ++i)
					{
						runTot /= *(m_z+j);
					}
					*(m_slaterArrayUp+d*nbrUp+j) = runTot;
				}
			}
		}
		//	Populate the spin-down Slater matrix with projected monopole functions for the spin down particles
		if(LLdown!=1)
		{
			for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
			{
				double m = m_preData->cfLz[d];
				for(int j=nbrUp; j<n; ++j)
				{
					(*(m_slaterArrayDown+d*nbrDown+j-nbrUp)).Set((double)m_preData->binomList[(int)(2*m_preData->effectiveMonopole*m_preData->maxBinom+m_preData->effectiveMonopole-m)]);
					int imax = m_preData->effectiveMonopole-m;
					for(int i=0; i<imax; ++i)
					{
						*(m_slaterArrayDown+d*nbrDown+j-nbrUp) *= *(m_z+j);
					}
				}
			}
	    for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<LLdown; ++ll, cumulativeDegen+=degen)
        {
			degen = m_preData->llOccupationNumbers[ll];
				for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
				{
					double m = m_preData->cfLz[d];
					for(int j=nbrUp; j<n; ++j)
					{
						runTot.Set(0.0);
						tmpCmplx1.Set(*(v+j));
						tmpCmplx2.Set(*(u+j));
						for(int s=0; s<=ll; ++s)
						{
							//	skip if we are multiplying by a 0 binomial coefficient
							if(m_preData->effectiveMonopole+ll-m-s>=0 && m_preData->effectiveMonopole+ll-m-s<=2*m_preData->effectiveMonopole+ll)
							{
								term.MulBy(m_preData->binomList[ll*m_preData->maxBinom+s]*m_preData->binomList[(int)((2*m_preData->effectiveMonopole+ll)*m_preData->maxBinom+m_preData->effectiveMonopole+ll-m-s)], *(m_pArray+(s*LLdown+ll-s)*n+j));
								int imax = ll-s;
								for(int i=0; i<imax; ++i)
								{
									term *= tmpCmplx1;
								}
								for(int i=0; i<s; ++i)
								{
									term *= tmpCmplx2;
								}
								//This line implements (-1)^s
								(s & 1)	? (runTot-=term) : (runTot+=term);
							}
						}
						int imax = m_preData->effectiveMonopole-m;
						for(int i=0; i<imax; ++i)
						{
							runTot *= *(m_z+j);
						}
						for(int i=imax; i<0; ++i)
						{
							runTot /= *(m_z+j);
						}
						*(m_slaterArrayDown+d*nbrDown+j-nbrUp) = runTot;
					}
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//	Print out matrix elements
		for(int i=0; i<nbrUp; ++i){
		for(int j=0; j<nbrUp; ++j){
		utilities::cout.DebuggingInfo() << "Yup[" << i << "][" << j << "]=" 
		                                << (*(m_slaterArrayUp+i*nbrUp+j)).Get() << std::endl;
		}}
		getchar();
		if(LLdown!=1){
		for(int i=0; i<nbrDown; ++i){
		for(int j=0; j<nbrDown; ++j){
		utilities::cout.DebuggingInfo() << "Ydown[" << i << "][" << j << "]=" 
		                                << (*(m_slaterArrayDown+i*nbrDown+j)).Get() << std::endl;
		}}getchar();}
		#endif
		////////////////////////////////////////////////////////////////////////////////
		//	The final step is to calculate the determinants
		dcmplx wfValue = utilities::linearAlgebra::LogDeterminant<U>(m_slaterArrayUp, nbrUp);
		if(LLdown!=1)
		{
			wfValue += utilities::linearAlgebra::LogDeterminant<U, P>(m_slaterArrayDown, nbrDown);
		}
		//	Multiply by a Jastrow factor
		for (int i=0; i<n-1; ++i)
		{
			for (int j=i+1; j<n; ++j)
			{
				wfValue += 2.0*m_cfData->halfFluxAttatched*log(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}
		//	Multiply by a Product_j u_j^(2Q) factor, which we removed from the determinant
		for(int i=0; i<nbrUp; ++i)
		{
			wfValue += 2.0*m_preData->effectiveMonopole*(log(*(u+i)));
		}
		if(LLdown!=1)
		{
			for(int i=nbrUp; i<n; ++i)
			{
				wfValue += 2.0*m_preData->effectiveMonopole*(log(*(u+i)));
			}
		}
		return wfValue;
	}

	////////////////////////////////////////////////////////////////////////////////
	//!	\brief	This function evaluates the unprojected CF wave functions i.e. IQHE	
	//!	 wave functions times Jastrow factors		
	//!
	//! There are two steps: first we populate a n x n Slater matrix with
	//! unprojected monopole harmonics, then we calculate the log of its
	//! Slater determinant.
	//!	
	//! \return log of wave function (NOT NORMALIZED)
	////////////////////////////////////////////////////////////////////////////////
	dcmplx PolarisedUnprojected(
		const int n,	        //!<	Number of particles
		const int nbrLLs,		//!<	Number of full Landau levels in construction
		dcmplx* u,		        //!<	Address of spinor u co-ordinate array	
		dcmplx* v)		        //!<	Address of spinor v co-ordinate array	
	{
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		++m_noCalls;
		#endif
		////////////////////////////////////////////////////////////////////////////////
		if(nbrLLs==1)
		{
			//	Use Laughlin 1/3 wave function
            dcmplx wfValue = 0.0;
			for(int i=0; i<n-1; ++i)
			{
				for(int j=i+1; j<n; ++j)
				{
					wfValue += log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
			return wfValue;
		}
		dcmplx z[n];
		for(int j=0; j<n; ++j)
		{
			z[j]=v[j]/u[j];
		}
		//	Populate the Slater matrix with unprojected monopole functions
		dcmplx slaterArray[n*n];
		for(int d=0; d<m_preData->llOccupationNumbers[0]; ++d)
		{
			double m = m_preData->cfLz[d];
			for(int j=0; j<n; ++j)
			{
				*(slaterArray+d*n+j) = std::pow(*(z+j),(m_preData->effectiveMonopole-m))*(double)m_preData->binomList[(int)(2*m_preData->effectiveMonopole*m_preData->maxBinom+m_preData->effectiveMonopole-m)];
			}
		}
        for(int ll=1, degen=0, cumulativeDegen = m_preData->llOccupationNumbers[0]; ll<nbrLLs; ++ll, cumulativeDegen+=degen)
        {
			degen = m_preData->llOccupationNumbers[ll];
			for(int d=cumulativeDegen; d<cumulativeDegen+degen; ++d)
			{
				double m = m_preData->cfLz[d];
				for(int j=0; j<n; ++j)
				{
					dcmplx runTot = 0.0;
					for(int s=0; s<=ll; ++s)
					{
						//	skip if we are multiplying by a 0 binomial coefficient
						if(m_preData->effectiveMonopole+ll-m-s>=0 && (m_preData->effectiveMonopole+ll-m-s)<=(2*m_preData->effectiveMonopole+ll))
						{
							dcmplx term = m_preData->binomList[(int)((2*m_preData->effectiveMonopole+ll)*m_preData->maxBinom+m_preData->effectiveMonopole+ll-m-s)];
							term *= pow(abs(*(z+j)), -2.0*s)*m_preData->binomList[ll*m_preData->maxBinom+s];
							//This line implements (-1)^s
							(s & 1)	? (runTot-=term) : (runTot+=term);
						}
					}
					*(slaterArray+d*n+j) = runTot*pow(*(z+j),(int)(m_preData->effectiveMonopole-m))*pow(abs(v[j]),2.0*ll);
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		#if _DEBUG_
		//	Print out matrix elements
		for(int i=0; i<n; ++i){
		for(int j=0; j<n; ++j){
		utilities::cout.DebuggingInfo() << "Y[" << i << "][" << j << "]=" 
		                                << *(slaterArray+i*n+j) << std::endl;
		}}
		getchar();
		#endif
		////////////////////////////////////////////////////////////////////////////////
		//	The final step is to calculate the determinant
		dcmplx wfValue = utilities::linearAlgebra::LogDeterminant<dcmplx>(slaterArray, n);
		//	Multiply by a Product_j u_j^(2Q) factor, which we removed from the determinant
		for(int i=0; i<n; ++i)
		{
			wfValue += 2.0*m_preData->effectiveMonopole*(log(*(u+i)));
		}
		return wfValue;
	}
};

////////////////////////////////////////////////////////////////////////////////
//!	\brief	This class defines the basics data associated with, and functions 
//!	to evaluate, composite fermion type FQHE wave functions.
//!	
//!	For more information see "Composite Fermions" by J.K. Jain.
////////////////////////////////////////////////////////////////////////////////
class CompositeFermion : public WaveFunction
{
	private:
	CompositeFermionData *m_cfData;	//!<	Class internal instance of a data structure 
									//!		to store CF wave function data
	PreData *m_preData;				//!<	Class internal instance of a data structure 
									//!		to store pre-calculated data
	ArbitraryPrecisionCfWavefunctions<dcmplx, double, 64> *m_lowPrec;   
	                                //!<	complex<double> implementation
	ArbitraryPrecisionCfWavefunctions<mpfCmplx, mpf_t, _PREC_LEVEL1_> *m_precLevel1;	
	                                //!<	mpf complex _PREC_LEVEL1_ implementation
	ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL2_> *m_precLevel2;	
	                                //!<	mpf complex _PREC_LEVEL2_ implementation
	ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL3_> *m_precLevel3;	
	                                //!<	mpf complex _PREC_LEVEL3_ implementation
	ArbitraryPrecisionCfWavefunctions<mpfCmplx, mpf_t, _PREC_LEVEL4_> *m_precLevel4;	
	                                //!<	mpf complex _PREC_LEVEL4_ implementation
	//	For most cases, do not need to consider higher precision
	//	(a warning message is produced instead)
    public:
	CompositeFermion(WaveFunctionData*, CompositeFermionData*);
	~CompositeFermion();	
	////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_SPHERE_GEOMETRY_
	bool IsCloseToPole(const int n, dcmplx* u, dcmplx* v) const;
	void RotateAwayFromPole(const int n, dcmplx* u, dcmplx* v) const;
	dcmplx EvaluateWfSphere(const int n, dcmplx* u, dcmplx* v) const;
	#endif
	////////////////////////////////////////////////////////////////////////////////
							
	////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_DISC_GEOMETRY_					
	dcmplx EvaluateWfDisc(const int n,dcmplx*) const;			
	#endif
	////////////////////////////////////////////////////////////////////////////////
	bool IncreasePrecCondition(const dcmplx lowPrecResult, const dcmplx highPrecResult, 
	                           int& maxPrecLevel) const;
	void GeneratePreFactorNumber(const int n, const int nbrLLs, int* nbrEachLL);
	void SetLLdata(const int n, const int nbrLL, const int* cfLz, const double* llOccupationNumbers);
};
}   //  End FQHE namespace
#endif
