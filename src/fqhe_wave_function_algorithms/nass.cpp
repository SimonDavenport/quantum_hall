////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This cpp file contains functions to evaluate the NASS wave function
//!		For information about the NASS state see Nucl. Phys. B 607, 549--576
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

#include "nass.hpp"

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the NonAbelianSpinSinglet class. 
//!
//!	 Assigns basic wave function characteristics. 
//!	 Reads in a file containing the unique permutations of the sets of particles
//!	 forming the NASS state and generates a list of minimal transpositions
//!	 relating these permutations. This information is used to optimize the NASS
//!	 wave function algorithm.
////////////////////////////////////////////////////////////////////////////////
FQHE::NonAbelianSpinSinglet::NonAbelianSpinSinglet(
	FQHE::WaveFunctionData *wavefunctionData)		//!<	Address of a struct of wave function data
	:  FQHE::WaveFunction(wavefunctionData)
{
	m_wfData->fillNumerator = 4;
	m_wfData->fillDenominator = 4*m_wfData->jastrowExponent+3;
	m_wfData->shift = 2+m_wfData->jastrowExponent;
	if(m_wfData->jastrowExponent%2==0)			m_wfData->statistics=_BOSONS_;
	else										m_wfData->statistics=_FERMIONS_;
	m_wfData->fillingFactor = (double)m_wfData->fillNumerator/m_wfData->fillDenominator;
	m_wfData->flux = (double)m_wfData->nbr/(m_wfData->fillingFactor)-m_wfData->shift;
	m_wfData->monopoleStrength = (m_wfData->nbr*m_wfData->fillDenominator/m_wfData->fillNumerator-m_wfData->shift)/2.0;
	m_wfData->wfName = "non-abelian spin singlet (NASS)";
	if(m_wfData->nbr%4==1)
	{
		std::cerr << "\n\tERROR WITH NASS WAVE FUNCTION: n must be divisible by 4." << std::endl;
		exit(EXIT_FAILURE);
	} 
	m_wfData->GenerateFileName();
    m_wfData->InitRadius();
	m_wfData->CheckAndPrint();
	//	For the NASS state we need to pre-process the sets of permutations
	//	Read in file of unique permutations of the NASS state (stored as e.g 12
	//	indicating the transposition 1<->2)
	std::stringstream permsFileName;	//	Name of file containing list of permutations
	std::ifstream f_perms;
	permsFileName.str("");
	permsFileName << "nassPermutationData/permutations" << m_wfData->nbr << ".dat";
	f_perms.open(permsFileName.str().c_str());
	if(f_perms.is_open()==0)
	{
		std::cout << "\tERROR: Cannot open file: " << permsFileName.str() << std::endl;
		exit(EXIT_FAILURE);
	}
	int permsNumber[10]={2, 6, 20, 70, 252, 924, 3432, 12870, 48620, 184756};
	m_nbrPerms = permsNumber[m_wfData->nbr/4-1];
	m_nOverTwo = m_wfData->nbr/2;
	m_permsList = new int[m_nbrPerms*m_nOverTwo];
	for(int i=0; i<m_nbrPerms; ++i)
	{
		for(int j=0; j<m_nOverTwo; ++j)
		{
			f_perms >> m_permsList[i*m_nOverTwo+j];
			m_permsList[i*m_nOverTwo+j] -= 1;
		}
	}
	f_perms.close();
	//	Calculate the lists of index shifts from one term to the next in the NASS state.
	//	The a and b clusters are done separately at the moment.
	//	Index shifts for the a cluster:
	m_sizeListA = new int[m_nbrPerms-1];
	int ktot = 0;
	for(int i=0; i<m_nbrPerms-1; ++i)
	{
		int k = 0;
		for(int j=0; j<m_wfData->nbr/4; ++j)
		{
			int dperms = m_permsList[i*m_nOverTwo+j]-m_permsList[(i+1)*m_nOverTwo+j];
			if(dperms!=0)	++k;
		}
		ktot += k;
		m_sizeListA[i] = k;
	}
	m_indexA1 = new int[ktot];
	m_indexA2 = new int[ktot];
	int k = 0;
	for(int i=0; i<m_nbrPerms-1; ++i)
	{
		for(int j=0; j<m_wfData->nbr/4; ++j)
		{
			int dperms = m_permsList[i*m_nOverTwo+j]-m_permsList[(i+1)*m_nOverTwo+j];
			if (dperms!=0)
			{
				m_indexA1[k] = m_permsList[i*m_nOverTwo+j];
				m_indexA2[k] = m_permsList[(i+1)*m_nOverTwo+j];
				++k;
			}
		}
	}
	//	Index shifts for the b cluster:
	m_sizeListB = new int[m_nbrPerms-1];
	ktot = 0;
	for(int i=0; i<m_nbrPerms-1; ++i)
	{	
		k = 0;
		for(int j=m_wfData->nbr/4; j<m_wfData->nbr/2; ++j)
		{
			int dperms = m_permsList[i*m_nOverTwo+j]-m_permsList[(i+1)*m_nOverTwo+j];
			if (dperms!=0)	++k;
		}
		ktot = ktot+k;
		m_sizeListB[i] = k;
	}
	m_indexB1 = new int[ktot];
	m_indexB2 = new int[ktot];
	k = 0;
	for(int i=0; i<m_nbrPerms-1; ++i)
	{
		for(int j=m_wfData->nbr/4; j<m_wfData->nbr/2; ++j)
		{
			int dperms = m_permsList[i*m_nOverTwo+j]-m_permsList[(i+1)*m_nOverTwo+j];
			if(dperms!=0)
			{
				m_indexB1[k] = m_permsList[i*m_nOverTwo+j];
				m_indexB2[k] = m_permsList[(i+1)*m_nOverTwo+j];
				++k;
			}
		}
	}
	//	Initialise variables which depend on n
	m_cmplxCtr1 = new dcmplx[m_nbrPerms/2];
	m_cmplxCtr2 = new dcmplx[m_nbrPerms];
	m_zwDiff = new dcmplx[m_nOverTwo*m_nOverTwo];
	//	in fact, we will only need access to the first m_nbrPerms/2 permutations.
	delete[] m_permsList;
	f_perms.open(permsFileName.str().c_str());
	m_permsList = new int[m_nbrPerms*m_nOverTwo/2];
	for(int i=0; i<m_nbrPerms/2; ++i)
	{
		for(int j=0; j<m_nOverTwo; ++j)
		{
			f_perms >> m_permsList[i*m_nOverTwo+j];
			//	N.B for c++ arrays, we now have to reduce the indices by 1.
			m_permsList[i*m_nOverTwo+j] -= 1;
		}
	}
	f_perms.close();
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the NonAbelianSpinSinglet class. 
//!
//!	 Deallocates memory assigned to store the various sets of permutations
//!	 and transpositions that were used to optimize the NASS wave function
//!	 algorithm.
////////////////////////////////////////////////////////////////////////////////
FQHE::NonAbelianSpinSinglet::~NonAbelianSpinSinglet()
{
	delete[] m_permsList;
	delete[] m_sizeListA;
	delete[] m_sizeListB;
	delete[] m_indexA1;
	delete[] m_indexA2;
	delete[] m_indexB1;
	delete[] m_indexB2;
	delete[] m_cmplxCtr1;
	delete[] m_cmplxCtr2;
	delete[] m_zwDiff;
}

#if _ENABLE_SPHERE_GEOMETRY_
////////////////////////////////////////////////////////////////////////////////
//! \brief 	NASS Wave function algorithm in the sphere geometry
//!	
//!	(C) Simon C Davenport 
//!
//!	First complete form appeared in June 2011
////////////////////////////////////////////////////////////////////////////////
dcmplx FQHE::NonAbelianSpinSinglet::EvaluateWfSphere(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* u,	    //!<	The memory address of the start of an array of spinor u coordinates
	dcmplx* v)	    //!<	The memory address of the start of an array of spinor v coordinates
	const
	{
	    dcmplx wfValue = dcmplx(0, 0);
	    //	Multiply by the Jastrow factor (z_i-z_j)^m_wfData->jastrowExponent
	    if(m_wfData->jastrowExponent>0)
	    {
		    for(int i=0; i<n-1; ++i)
		    {
			    for(int j=i+1; j<n; ++j)
			    {
				    wfValue += log(*(u+i)**(v+j)-*(u+j)**(v+i));
			    }
		    }
		    wfValue *= m_wfData->jastrowExponent;
	    }
	    //	The NASS wavefunction is a sum over permutations of some Halperin 221 factors
	    //	Spin up co-ordinates are z and spin dosn are w
	    //	First, we generate an array of the parts of each term which change under 
	    //	a permutation of z, keeping w fixed or a permutation of w keeping z fixed
	    int* p_currPerm = m_permsList;
	    for(int k=0; k<m_nbrPerms/2; ++k, p_currPerm+=m_nOverTwo)
	    {
		    dcmplx cmplxCtrA = 0.0;
		    dcmplx cmplxCtrB = 0.0;
		    for(int i=0; i<(n/4-1); ++i)
		    {
			    int a1 = *(p_currPerm+i);
			    int a2 = *(p_currPerm+i+n/4);
			    for(int j=i+1; j<n/4; ++j)
			    {
				    int b1 = *(p_currPerm+j);
				    int b2 = *(p_currPerm+j+n/4);
				    cmplxCtrA += log(*(u+a1)**(v+b1)-*(u+b1)**(v+a1))+log(*(u+a2)**(v+b2)-*(u+b2)**(v+a2));
				    b1 += n/2;
				    b2 += n/2;
				    cmplxCtrB += log(*(u+a1+n/2)**(v+b1)-*(u+b1)**(v+a1+n/2))+log(*(u+a2+n/2)**(v+b2)-*(u+b2)**(v+a2+n/2));
			    }
		    }
		    *(m_cmplxCtr1+k) = 2.0*cmplxCtrA;
		    *(m_cmplxCtr2+k) = 2.0*cmplxCtrB;
		    *(m_cmplxCtr2+(k+m_nbrPerms/2)) = 2.0*cmplxCtrB;
	    }
	    //	Now we generate an array of the parts of each terms which change under 
	    //	permutations of both z and w
	    //	Work out a matrix of differences between w and z, to reduce the number of times
	    //	we access the memory.
	    for(int i=0; i<m_nOverTwo; ++i)
	    {
		    for(int j=0; j<m_nOverTwo; ++j)
		    {
			    m_zwDiff[i*m_nOverTwo+j] = log(*(u+m_nOverTwo+i)**(v+j)-*(u+j)**(v+m_nOverTwo+i));
			    //std::cout << "m_zwDiff[" << i << "][" << j << "]=" << m_zwDiff[i*m_nOverTwo+j] << std::endl;
		    }
	    }
	    //	Start by fully calculating the first term:
	    dcmplx cmplxCtrA = 0.0;
	    dcmplx* currDiff1 = m_zwDiff;
	    dcmplx* currDiff2 = m_zwDiff+n/4*m_nOverTwo;
	    for(int i=0; i<n/4; ++i, currDiff1 += m_nOverTwo, currDiff2 += m_nOverTwo)
	    {
            for(int j=0; j<n/4; ++j)
		    {
			    cmplxCtrA += *(currDiff1+j)+*(currDiff2+j+n/4);
		    }
	    }
	    dcmplx cmplxCtrE = cmplxCtrA;
	    dcmplx cmplxCtrB = cmplxCtrE;
	    p_currPerm = m_permsList;
	    dcmplx* p_cmplxCtr1 = m_cmplxCtr1;
	    dcmplx*p_zPerm_sizeA = m_sizeListA;
	    dcmplx*p_zPerm_sizeB = m_sizeListB;
	    dcmplx*p_zPerm_indexA1 = m_indexA1;
	    dcmplx*p_zPerm_indexA2 = m_indexA2;
	    dcmplx*p_zPerm_indexB1 = m_indexB1;
	    dcmplx*p_zPerm_indexB2 = m_indexB2;
	    for(int zPerm=0; zPerm<m_nbrPerms/2; ++zPerm, p_currPerm+=m_nOverTwo, ++p_zPerm_sizeA, ++p_zPerm_sizeB, ++p_cmplxCtr1)
	    {
		    dcmplx cmplxCtrC = exp(cmplxCtrB+*(m_cmplxCtr2));
		    p_cmplxCtr2 = m_cmplxCtr2+1;
		    p_wPerm_sizeA = m_sizeListA;
		    p_wPerm_sizeB = m_sizeListB;
		    p_wPerm_indexA1 = m_indexA1;
		    p_wPerm_indexA2 = m_indexA2;
		    p_wPerm_indexB1 = m_indexB1;
		    p_wPerm_indexB2 = m_indexB2;
		    int lmax;
		    dcmplx dcmplx;
		    for(int wPerm=0; wPerm<m_nbrPerms-1; ++wPerm, ++p_wPerm_sizeA, ++p_wPerm_sizeB, ++p_cmplxCtr2)
		    {
			    cmplxCtrD = 0.0;
			    //	Multiply by the ratios of adjacent terms,
			    //	for fixed z perm and changing w perm
			    lmax = *p_wPerm_sizeA;
			    for(int l=0; l<lmax; ++l)
			    {
				    currDiff1 = m_zwDiff+*(p_wPerm_indexA1)*m_nOverTwo;
				    currDiff2 = m_zwDiff+*(p_wPerm_indexA2)*m_nOverTwo;
				    ++p_wPerm_indexA1;
				    ++p_wPerm_indexA2;
				    dcmplx* p_tmp = p_currPerm;
				    for(int i=0; i<n/4; ++i, ++p_tmp)
				    {
					    cmplxCtrD += *(currDiff2+*(p_tmp))-*(currDiff1+*(p_tmp));
				    }
			    }						
			    lmax = *p_wPerm_sizeB;
			    for(int l=0; l<lmax; ++l)
			    {
				    currDiff1 = m_zwDiff+*(p_wPerm_indexB1)*m_nOverTwo;
				    currDiff2 = m_zwDiff+*(p_wPerm_indexB2)*m_nOverTwo;
				    ++p_wPerm_indexB1;
				    ++p_wPerm_indexB2;
				    dcmplx* p_tmp = p_currPerm;
				    for(int i=n/4; i<n/2; ++i, ++p_tmp)
				    {
					    cmplxCtrD += *(currDiff2+*(p_tmp))-*(currDiff1+*(p_tmp));
				    }
			    }
			    //	Update the value of the current term:
			    cmplxCtrB += cmplxCtrD;
			    //	Udpate cmplxCtrC, the value of the term times the w-only part:
			    cmplxCtrC += exp(cmplxCtrB+*p_cmplxCtr2);
		    }
		    //	Calculate running total using logarithmic addition:
		    dcmplx dotProd = *p_cmplxCtr1+log(cmplxCtrC);
		    if(zPerm == 0)
		    {
			    runTot = dotProd;
		    }
		    else
		    {
			    runTot = log(exp(runTot-dotProd)+dcmplx(1, 0))+dotProd;
		    }
		    //	Multiply by the ratios of adjacent terms,
		    //	for fixed w perm and changing z perm:
		    lmax = *p_zPerm_sizeA;
		    cmplxCtrD = 0.0;
		    for(int l=0; l<lmax; ++l)
		    {
			    int a1 = *p_zPerm_indexA1;
			    int a2 = *p_zPerm_indexA2;
			    ++p_zPerm_indexA1;
			    ++p_zPerm_indexA2;
			    currDiff1 = m_zwDiff;
			    for(int i=0; i<n/4; ++i)
			    {
				    cmplxCtrD += *(currDiff1+a2)-*(currDiff1+a1);
				    currDiff1 += m_nOverTwo;
			    }
		    }
		    lmax = *p_zPerm_sizeB;
		    for(int l=0; l<lmax; ++l)
		    {
			    b1 = *p_zPerm_indexB1;
			    b2 = *p_zPerm_indexB2;
			    ++p_zPerm_indexB1;
			    ++p_zPerm_indexB2;
			    currDiff2 = m_zwDiff+n/4*m_nOverTwo;
			    for(int i=0; i<n/4; ++i)
			    {
				    cmplxCtrD += *(currDiff2+b2)-*(currDiff2+b1);
				    currDiff2 += m_nOverTwo;
			    }
		    }	
		    //	Update the value of the current term:
		    cmplxCtrE += cmplxCtrD;
		    //	Reset CmplxCtrB for the next set of terms to be calculated:
		    cmplxCtrB = cmplxCtrE;
	    }
	    wfValue += runTot;
	return wfValue;
}
#endif

#if _ENABLE_DISC_GEOMETRY_
////////////////////////////////////////////////////////////////////////////////
//! \brief 	NASS Wave function algorithm in the disc geometry
//!	
//!	(C) Simon C Davenport 
//!
//!	First complete form appeared in June 2011
////////////////////////////////////////////////////////////////////////////////
dcmplx FQHE::NonAbelianSpinSinglet::EvaluateWfDisc(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* z)	    //!<	The memory address of the start of an array of disc coordinates
	const
{
	//////////////////////////////////////////////////////////////////////////////
	//	declare local variables													//
	//////////////////////////////////////////////////////////////////////////////
	int i,j,k,l;		//	generic loop counters
	int zPerm,wPerm;	//	index of the z or w permutation
	dcmplx cmplxCtrA,cmplxCtrB,cmplxCtrC,cmplxCtrD,cmplxCtrE;
	//		Some Counters to organise the data in the NASS algorithm
	dcmplx	runTot;		//		Store running total	
	dcmplx	dotProd;	//		Store dot product between vector of z only
						//		permutations with matrix of z	
	int lmax;
	int a1,a2,b1,b2;
	int *p_tmp;
	int *p_currPerm;
	dcmplx	*p_cmplxCtr1,*p_cmplxCtr2;
	int *p_zPerm_sizeA,*p_wPerm_sizeA;
	int *p_zPerm_sizeB,*p_wPerm_sizeB;
	int *p_zPerm_indexA1,*p_wPerm_indexA1;
	int *p_zPerm_indexA2,*p_wPerm_indexA2;
	int *p_zPerm_indexB1,*p_wPerm_indexB1;
	int *p_zPerm_indexB2,*p_wPerm_indexB2;
	dcmplx *currDiff1,*currDiff2;
	dcmplx wfValue;
	//////////////////////////////////////////////////////////////////////////////
	
	dcmplx wfValue = dcmplx(0, 0);
	//	Multiply by the Jastrow factor (z_i-z_j)^m_wfData->jastrowExponent
	if(m_wfData->jastrowExponent>0)
	{
		for(int i=0; i<n-1; ++i)
		{
			for(int j=i+1; j<n; ++j)
			{
				wfValue += log(*(z+i)-*(z+j));
			}
		}
		wfValue *= m_wfData->jastrowExponent;
	}
	//	The NASS wavefunction is a sum over permutations of some Halperin 221 factors
	//	Spin up co-ordinates are z and spin dosn are w
	//	First, we generate an array of the parts of each term which change under 
	//	a permutation of z, keeping w fixed or a permutation of w keeping z fixed
	int* p_currPerm=m_permsList;
	for(int k=0; k<m_nbrPerms/2; ++k, p_currPerm+=m_nOverTwo)
	{
		dcmplx cmplxCtrA = 0.0;
		dcmplx cmplxCtrB = 0.0;
		for(int i=0; i<(n/4-1); ++i)
		{
			int a1 = *(p_currPerm+i);
			int a2 = *(p_currPerm+i+n/4);
			for(int j=i+1; j<n/4; ++j)
			{
				int b1 = *(p_currPerm+j);
				int b2 = *(p_currPerm+j+n/4);
				cmplxCtrA += log(*(z+a1)-*(z+b1))+log(*(z+a2)-*(z+b2));
				b1 += n/2;
				b2 += n/2;
				cmplxCtrB += log(*(z+a1+n/2)-*(z+b1))+log(*(z+a2+n/2)-*(z+b2));
			}
		}
		*(m_cmplxCtr1+k) = 2.0*cmplxCtrA;
		*(m_cmplxCtr2+k) = 2.0*cmplxCtrB;
		*(m_cmplxCtr2+(k+m_nbrPerms/2)) = 2.0*cmplxCtrB;
	}
	//	Now we generate an array of the parts of each terms which change under 
	//	permutations of both z and w
	//	Work out a matrix of differences between w and z, to reduce the number of times
	//	we access the memory.
	for(int i=0; i<m_nOverTwo; ++i)
	{
		for(int j=0; j<m_nOverTwo; ++j)
		{
			m_zwDiff[i*m_nOverTwo+j] = log(*(z+m_nOverTwo+i)-*(z+j));
			//std::cout<<"m_zwDiff["<<i<<"]["<<j<<"]="<<m_zwDiff[i*m_nOverTwo+j]<<std::endl;
		}
	}
	//	Start by fully calculating the first term:
	cmplxCtrA = 0.0;
	dcmplx* currDiff1 = m_zwDiff;
	dcmplx* currDiff2 = m_zwDiff+n/4*m_nOverTwo
	for(int i=0; i<n/4; ++i, currDiff1+=m_nOverTwo, currDiff2+=m_nOverTwo)
	{
        for(int j=0; j<n/4; ++j)
		{
			cmplxCtrA += *(currDiff1+j)+*(currDiff2+j+n/4);
		}
	}
	dcmplx cmplxCtrE = cmplxCtrA;
	dcmplx cmplxCtrB = cmplxCtrE;
	p_currPerm = m_permsList;
	dcmplx* p_cmplxCtr1 = m_cmplxCtr1;
	dcmplx* p_zPerm_sizeA = m_sizeListA;
	dcmplx* p_zPerm_sizeB = m_sizeListB;
	dcmplx* p_zPerm_indexA1 = m_indexA1;
	dcmplx* p_zPerm_indexA2 = m_indexA2;
	dcmplx* p_zPerm_indexB1 = m_indexB1;
	dcmplx* p_zPerm_indexB2 = m_indexB2;
	for(int zPerm=0; zPerm<m_nbrPerms/2; ++zPerm, p_currPerm += m_nOverTwo, ++p_zPerm_sizeA, ++p_zPerm_sizeB, ++p_cmplxCtr1)
	{
		dcmplx cmplxCtrC = exp(cmplxCtrB+*(m_cmplxCtr2));
		dcmplx* p_cmplxCtr2 = m_cmplxCtr2+1;
		dcmplx* p_wPerm_sizeA = m_sizeListA;
		dcmplx* p_wPerm_sizeB = m_sizeListB;
		dcmplx* p_wPerm_indexA1 = m_indexA1;
		dcmplx* p_wPerm_indexA2 = m_indexA2;
		dcmplx* p_wPerm_indexB1 = m_indexB1;
		dcmplx* p_wPerm_indexB2 = m_indexB2;
		int lmax;
		for(int wPerm=0; wPerm<m_nbrPerms-1; ++wPerm, ++p_wPerm_sizeA, ++p_wPerm_sizeB, ++p_cmplxCtr2)
		{
			dcmplx cmplxCtrD = 0.0;
			//	Multiply by the ratios of adjacent terms,
			//	for fixed z perm and changing w perm
			lmax = *p_wPerm_sizeA;
			for(int l=0; l<lmax; ++l)
			{
				currDiff1 = m_zwDiff+*(p_wPerm_indexA1)*m_nOverTwo;
				currDiff2 = m_zwDiff+*(p_wPerm_indexA2)*m_nOverTwo;
				++p_wPerm_indexA1;
				++p_wPerm_indexA2;
				dcmplx* p_tmp = p_currPerm;
				for(int i=0; i<n/4; ++i, ++p_tmp)
				{
					cmplxCtrD += *(currDiff2+*(p_tmp))-*(currDiff1+*(p_tmp));
				}
			}						
			lmax = *p_wPerm_sizeB;
			for(int l=0; l<lmax; ++l)
			{
				currDiff1 = m_zwDiff+*(p_wPerm_indexB1)*m_nOverTwo;
				currDiff2 = m_zwDiff+*(p_wPerm_indexB2)*m_nOverTwo;
				++p_wPerm_indexB1;
				++p_wPerm_indexB2;
				dcmplx* p_tmp = p_currPerm;
				for(int i=n/4; i<n/2; ++i, ++p_tmp)
				{
					cmplxCtrD += *(currDiff2+*(p_tmp))-*(currDiff1+*(p_tmp));
				}
			}
			//	Update the value of the current term:
			cmplxCtrB += cmplxCtrD;
			//	Udpate cmplxCtrC, the value of the term times the w-only part:
			cmplxCtrC += exp(cmplxCtrB+*p_cmplxCtr2);
		}
		//	Calculate running total using logarithmic addition:
		dcmplx dotProd = *p_cmplxCtr1+log(cmplxCtrC);
		if(zPerm==0)
		{
			runTot = dotProd;
		}
		else
		{
			runTot = log(exp(runTot-dotProd)+dcmplx(1, 0))+dotProd;
		}
		//	Multiply by the ratios of adjacent terms,
		//	for fixed w perm and changing z perm:
		lmax = *p_zPerm_sizeA;
		cmplxCtrD = 0.0;
		for(int l=0; l<lmax; ++l)
		{
			int a1 = *p_zPerm_indexA1;
			int a2 = *p_zPerm_indexA2;
			++p_zPerm_indexA1;
			++p_zPerm_indexA2;
			currDiff1 = m_zwDiff;
			for(int i=0; i<n/4; ++i)
			{
				cmplxCtrD += *(currDiff1+a2)-*(currDiff1+a1);
				currDiff1 += m_nOverTwo;
			}
		}
		lmax = *p_zPerm_sizeB;
		for(int l=0; l<lmax; ++l)
		{
			int b1 = *p_zPerm_indexB1;
			int b2 = *p_zPerm_indexB2;
			++p_zPerm_indexB1;
			++p_zPerm_indexB2;
			currDiff2 = m_zwDiff+n/4*m_nOverTwo;
			for(int i=0; i<n/4; ++i)
			{
				cmplxCtrD += *(currDiff2+b2)-*(currDiff2+b1);
				currDiff2 += m_nOverTwo;
			}
		}
		//	Update the value of the current term:
		cmplxCtrE += cmplxCtrD;
		//	Reset CmplxCtrB for the next set of terms to be calculated:
		cmplxCtrB = cmplxCtrE;
	}
	wfValue += runTot;
	//	Gaussian exponential factors
	wfValue *= m_wfData->jastrowExponent;
	for(int i=0; i<n; ++i)
	{
		double tmp = abs(*(z+i));
		wfValue -= tmp*tmp/4;	
	}
	return wfValue;
}
#endif
