////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 04/12/2013
//!
//!  \file 
//!		This cpp file contains functions to evaluate the Moore-Read (Pfaffian)
//! 	wave function.
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

//!
//!	\todo Fix what happens for the filling factor 1 wave function (jastrowExponent=0)
//!	in the sphere geometry case.
//!

#include "moore_read.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the MooreRead class. 
//!
//!	 Assigns basic wave function characteristics. 
//!	 Allocates memory to store the Pfaffian matrix.
//!
//!	 In the torus geometry case, we populate an array of possible electron-electron
//!	 theta functions. By the end of the constructor we have populated the following 
//!	 theta function table:
//!   - m_thetaFuncs To store a table of theta function values on the lattice
//!	    with no offset
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::MooreRead::MooreRead(
	FQHE::WaveFunctionData *wavefunctionData) 		//!<	Address of a struct of wave function data
	: FQHE::WaveFunction(wavefunctionData)
{
	//	Generate basic data 

	m_wfData->fillNumerator=1;
	m_wfData->fillDenominator=m_wfData->jastrowExponent;
	m_wfData->shift=m_wfData->jastrowExponent+1;
	
	//	Specify whether we are dealing with bosons or fermions
	if(m_wfData->jastrowExponent%2==1)	m_wfData->statistics=_BOSONS_;
	else								m_wfData->statistics=_FERMIONS_;
	
	//	Specify a default folder name
	m_wfData->wfName="Moore-Read";
	
	if(m_wfData->nbr%2==1)
	{
		std::cout<<"\n\tERROR: n must be even."<<std::endl;
		exit (1);
	}
	
	m_wfData->fillingFactor =(double)m_wfData->fillNumerator/m_wfData->fillDenominator;

	if(m_wfData->geometry==_SPHERE_ || m_wfData->geometry==_DISC_)				        		
	{
		m_wfData->flux=(double)m_wfData->nbr/(m_wfData->fillingFactor)-m_wfData->shift;
		m_wfData->monopoleStrength = (m_wfData->nbr*m_wfData->fillDenominator/m_wfData->fillNumerator-m_wfData->shift)/2.0;
	}
	
	m_wfData->flux -= m_wfData->nbrQh;

	//	Allocate memory to store Pfaffian matrix
	
	m_pfaff = new dcmplx[m_wfData->nbr*m_wfData->nbr];
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_
	
	if(m_wfData->geometry==_TORUS_)
	{
	
		//	Generate basic data 
		
		m_wfData->flux=(double)m_wfData->nbr/(m_wfData->fillingFactor);
		m_wfData->torusDegeneracy=3;
		m_wfData->fusionChannels=1;
		m_wfData->maxOccupation=3;
		
		if(m_wfData->nbrQh==4)
		{
			m_wfData->fusionChannels=2;
		}
		
		//  First generate theta function objects with no offsets.
		//	These will be used for evaluations of electron-electron Jastrow factors
		//	and therefore they are needed for both the ground state and quasi hole states
		
		m_thetaFuncs = new utilities::ThetaLookUp*[4];
	 
		m_thetaFuncs[0]=new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,-0.5,dcmplx(-smallShift,0.0),false,1);
		m_thetaFuncs[1]=new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,0.0,dcmplx(-smallShift,0.0),false,1);
		m_thetaFuncs[2]=new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.0,dcmplx(-smallShift,0.0),false,1);
		m_thetaFuncs[3]=new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.5,dcmplx(-smallShift,0.0),false,1);

		//  Set pointers to torus pre-calculation arrays to be null

		m_thetaFuncsQh        = 0;
		m_thetaFuncsW1W2      = 0;
		m_thetaFuncsW1W2W3W4A = 0;
		m_thetaFuncsW1W2W3W4B = 0;
		m_thetaFuncsCom       = 0;
		
		m_uniqueQhPositions   = 0;
		m_uniqueW1W2          = 0;
		m_uniqueW1W2W3W4A     = 0;
		m_uniqueW1W2W3W4B     = 0;
		m_uniqueCom           = 0;
	}

	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//  Now we can generate a file name
	m_wfData->GenerateFileName();

    //  Generate sphere/disc radius value
    m_wfData->InitRadius();

	//  Print a summary
	m_wfData->CheckAndPrint();

	//	FOR BENCHMARKING

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

	m_timeJastrow=0;
	m_timePfaffian=0;
	m_timeGauss=0;
	m_timeQuasihole=0;
	m_timeCM=0;

	#endif
	//////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	std::cout.precision(15);
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the MooreRead class. 
//!
//!	 Deallocates memory as required
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::MooreRead::~MooreRead()
{
	delete[] m_pfaff;

    //////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_

	if(m_wfData->geometry==_TORUS_)
	{
	
		delete[] m_thetaFuncs;

		if(m_thetaFuncsQh!=0)	            delete[] m_thetaFuncsQh;
		//std::cout<<" deleted m_thetaFuncsQh"<<std::endl;
		if(m_thetaFuncsW1W2!=0)	        	delete[] m_thetaFuncsW1W2;
		//std::cout<<" deleted m_thetaFuncsW1W2"<<std::endl;
		if(m_thetaFuncsW1W2W3W4A!=0)  		delete[] m_thetaFuncsW1W2W3W4A;
		//std::cout<<" deleted m_thetaFuncsW1W2W3W4A"<<std::endl;
		if(m_thetaFuncsW1W2W3W4B!=0)	    delete[] m_thetaFuncsW1W2W3W4B;
		//std::cout<<" deleted m_thetaFuncsW1W2W3W4B"<<std::endl;
		if(m_uniqueQhPositions!=0)	    	delete[] m_uniqueQhPositions;
		//std::cout<<" deleted m_uniqueQhPositions"<<std::endl;
		if(m_uniqueW1W2!=0)	            	delete[] m_uniqueW1W2;
		//std::cout<<" deleted m_uniqueW1W2"<<std::endl;
		if(m_uniqueW1W2W3W4A!=0)	        delete[] m_uniqueW1W2W3W4A;
		//std::cout<<" deleted m_uniqueW1W2W3W4A"<<std::endl;
		if(m_uniqueW1W2W3W4B!=0)	        delete[] m_uniqueW1W2W3W4B;
		//std::cout<<" deleted m_uniqueW1W2W3W4B"<<std::endl;
		if(m_thetaFuncsCom!=0)	        	delete[] m_thetaFuncsCom;
		//std::cout<<" deleted m_thetaFuncsCom"<<std::endl;
		if(m_uniqueCom!=0)	            	delete[] m_uniqueCom;
		//std::cout<<" deleted m_uniqueCom"<<std::endl;

	}
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

	std::cout<<"\t\t\t-----BENCHMARKING-----\n\n";
	std::cout<<"\tJASTROW TIME:\t"<<m_timeJastrow<<std::endl;
	std::cout<<"\tPFAFFIAN TIME:\t"<<m_timePfaffian<<std::endl;
	std::cout<<"\tGAUSS TIME:\t"<<m_timeGauss<<std::endl;
	std::cout<<"\tQUASIHOLE TIME:\t"<<m_timeQuasihole<<std::endl;
	std::cout<<"\tCM TERM TIME:\t"<<m_timeCM<<std::endl<<std::endl;

	#endif
	//////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Moore--Read wave function algorithm in the sphere geometry 
//!
//!	For more information see e.g. Nuclear Physics B360 (1991) 362-396 and do 
//!	a forward citation search!
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::MooreRead::EvaluateWfSphere(
	const int n,		//!<	The number of particles to perform the evaluation for
	dcmplx* u,	        //!<	The memory address of the start of an array of spinor u coordinates
	dcmplx* v)	        //!<	The memory address of the start of an array of spinor v coordinates
	const
{
	//	Local variables
	
	int i,j;			//	generic loop counters
	dcmplx *p_row;		//	Pointer to a row of the Pfaffian matrix
	dcmplx wfValue;
	
	wfValue=dcmplx(0,0);

	//	Calculate Jastrow factor part
	
	if(m_wfData->jastrowExponent>0)
	{
		for (i=0;i<n-1;i++)
		{
			for (j=i+1;j<n;j++)
			{
				wfValue+=log(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}
		wfValue*=m_wfData->jastrowExponent;
	}

	//	Calculate Pfaffian part

	//	Set upper triangular parts to be 1/(u[i]v[j]-u[j]v[i])
	
	for(i=0,p_row=m_pfaff;i<n;i++,p_row+=n)
	{
		for(j=i+1;j<n;j++)
		{
			*(p_row+j)=dcmplx(1,0)/(*(u+i)**(v+j)-*(u+j)**(v+i));
		}
	}
	
	//PrintPfaffian(n);

	wfValue+=utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,n);
	
	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////
#if _ENABLE_DISC_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Moore--Read wave function algorithm in the disc geometry 
//!
//!	For more information see e.g. Nuclear Physics B360 (1991) 362-396 and do 
//!	a forward citation search!
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::MooreRead::EvaluateWfDisc(
	const int n,		//!<	The number of particles to perform the evaluation for
	dcmplx* z)			//!<	The memory address of the start of an array of disc coordinates
	const
{
	//	Local variables
	
	int i,j;				//	generic loop counters
	dcmplx *p_row;			//	Pointer to a row of the Pfaffian matrix	
	dcmplx wfValue;
	
	wfValue=dcmplx(0,0);
	
	//	Calculate Jastrow factor part
	
	if(m_wfData->jastrowExponent>0)
	{
		for (i=0;i<n-1;i++)
		{
			for (j=i+1;j<n;j++)
			{
				wfValue+=log(*(z+i)-*(z+j));
			}
		}
		
		wfValue*=m_wfData->jastrowExponent;
	}
	
	//Calculate Pfaffian part
	
	//	Set upper triangular parts to be 1/(z[i]-z[j])
	
	for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
	{
		for(j=i+1;j<n;j++)
		{
			*(p_row+j)=dcmplx(1,0)/(*(z+i)-*(z+j));
		}
	}
	
	 //PrintPfaffian(n,);
	 
	wfValue+=utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,n);
	
	//	Gaussian exponential part
	for(i=0;i<n;i++)
	{
		double tmp=abs(*(z+i));
		wfValue-=tmp*tmp/4;	
	}
	
	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
#if _TEST_MODE_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	A function to print the values stored in the Pfaffian matrix
//!
//!	Used for testing only.
////////////////////////////////////////////////////////////////////////////////

void FQHE::MooreRead::PrintPfaffian(const int dim) //!<	Dimension of the matrix
const
{	
	std::cout<<"print Pfaffian matrix:"<<std::endl<<std::endl;
	
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			std::cout<<m_pfaff[i*dim+j]<<"\t";
		}
		std::cout<<std::endl;
	}
	
	getchar();
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_TORUS_GEOMETRY_

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function pre-calculates all the possible values of theta
//!	 functions that can be called to evaluate the wave function for a given set
//!	 of quasi hole positions.														
//! 
//!  At the end of the function we have populated various arrays:
//!
//!	  - m_thetaFuncsQh To store a list of tables of theta function values 
//!		offset by each unique quasi hole position
//!	  - m_thetaFuncsW1W2 To store a list of tables of theta function values 
//!     offset by values of (W1-W2)/2
//!	  - m_thetaFuncsW1W2W3W4A To store a list of tables of theta function values 
//!     offset by values of W1+W2-W3-W4)/2
//!	  - m_thetaFuncsW1W2W3W4B To store a list of tables of theta function values 
//!     offset by values of (W1-W2+W3-W4)/2
//!	  -	m_thetaFuncsCom To store a list of tables of theta function values
//!		offset by the centre of mass coordinate (minus the COM zero depending
//!		on whether we use the COM zeros formalism or not)
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::MooreRead::LatticeDataPrecalculation(
	const int nbrPos,		//!<	Number of quasi hole positions along the braid
	const int nbrQh,		//!<	Number of quasi holes   
	dcmplx* listQh)			//!<	The starting address of a list of all quasi-hole
							//!		positions that are visited along the path 
							//!		(can include duplicates). The size of this list
							//!		must equal nbrPos.
{
	//	Local variables

    int i,j,counter;
    int *keepIndexes;
	
    //////////////////////////////////////////////////////////////////////////////////
    //  Generate theta function objects offset by every unique possible value 
    //  of a quasi-hole position z_i      
            
	//	Determine a list of all quasi hole positions without repeated values

	keepIndexes = new int[nbrPos];

	for(i=0;i<nbrPos;i++)
	{
		keepIndexes[i]=i;
	}

	for(i=0;i<nbrPos;i++)
	{
		for(j=i+1;j<nbrPos;j++)
		{
			if(listQh[j]==listQh[i])
			{
				keepIndexes[j]=-1;
			}
		}
	}

	counter=0;

	for(i=0;i<nbrPos;i++)
	{
		if(keepIndexes[i]!=-1)	counter++;
	}

	m_nbrQhPositions=counter;

	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
		std::cout<<"\treduced number of positions is "<<m_nbrQhPositions<<std::endl;
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	m_uniqueQhPositions = new (std::nothrow) dcmplx[m_nbrQhPositions];

	counter=0;

	for(i=0;i<nbrPos;i++)
	{
		if(keepIndexes[i]!=-1)
		{
			m_uniqueQhPositions[counter]=listQh[i];
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_
			std::cout<<"\t filtered position "<<counter<<" is "<<m_uniqueQhPositions[counter]<<std::endl;
			#endif
			//////////////////////////////////////////////////////////////////////////////////
			
			counter++;
		}
	}

	delete[] keepIndexes;

    //  Generate theta function objects with the corresponding offsets

	m_thetaFuncsQh=new (std::nothrow) utilities::ThetaLookUp*[m_nbrQhPositions];

	for(i=0;i<m_nbrQhPositions;i++)
	{
		m_thetaFuncsQh[i] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,-0.5,m_uniqueQhPositions[i]-smallShift,true,1);
	}

    //////////////////////////////////////////////////////////////////////////////////
    //  For the 2 quasi-hole case determine all unique values of (w1-w2)/2 along the
    //  path and generate theta function objects with corresponding offsets
  
    if(nbrQh==2)
    {
        int origNbr = (int)nbrPos/nbrQh;
        
        //  First make a unique list of (w1-w2)/2 values
        
        keepIndexes = new int[origNbr];
        
        for(i=0;i<origNbr;i++)
	    {
		    keepIndexes[i]=i;
	    }
        
        for(i=0;i<origNbr;i++)
        {
            for(j=i+1;j<origNbr;j++)
            {

                if( abs((listQh[i*nbrQh]-listQh[i*nbrQh+1])/2.0 - (listQh[j*nbrQh]-listQh[j*nbrQh+1])/2.0) < sameTol )
                {
                    keepIndexes[j]=-1;
                }
            }
        }
        
        counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)	counter++;
	    }

	    m_nbrW1W2=counter;
        	    
	    //////////////////////////////////////////////////////////////////////////////////
	    #if _TEST_MODE_
	    std::cout<<"\treduced number of (w1-w2)/2 is "<<m_nbrW1W2<<std::endl;
		#endif
	    //////////////////////////////////////////////////////////////////////////////////
        
        m_uniqueW1W2 = new (std::nothrow) dcmplx[m_nbrW1W2];

	    counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)
		    {
			    m_uniqueW1W2[counter]=(listQh[i*nbrQh]-listQh[i*nbrQh+1])/2.0;
			    
			    //////////////////////////////////////////////////////////////////////////////////
			    #if _TEST_MODE_
			    std::cout<<"\t filtered (w1-w2)/2 "<<counter<<" is "<<m_uniqueW1W2[counter]<<std::endl;
				#endif
			    //////////////////////////////////////////////////////////////////////////////////
			    counter++;
		    }
	    }
        
        delete[] keepIndexes;
        
        //  Generate theta function objects with the corresponding offsets
        //  The first m_nbrW1W2 are for the Theta_2 functions, the next m_nbrW1W2 are for Theta_3
        //  and the third set are for Theta_4 functions

	    m_thetaFuncsW1W2 = new (std::nothrow) utilities::ThetaLookUp*[m_wfData->torusDegeneracy*m_nbrW1W2];

        counter=0;

	    for(i=0;i<m_nbrW1W2;i++)    //  Theta_2 objects
	    {
		    m_thetaFuncsW1W2[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,0.0,-m_uniqueW1W2[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2;i++)    //  Theta_3 objects
	    {
		    m_thetaFuncsW1W2[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.0,-m_uniqueW1W2[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2;i++)    //  Theta_4 objects
	    {
		    m_thetaFuncsW1W2[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.5,-m_uniqueW1W2[i]-smallShift,false,1);
		    
		    counter++;
	    }
    }
    
    //////////////////////////////////////////////////////////////////////////////////
    //  For the 4 quasi-hole case determine all unique values of (w1+w2-w3-w4)/2 
    //  and (w1-w2+w3-w4)/2 along the path
    
    if(nbrQh==4)
    {
        int origNbr = (int)nbrPos/nbrQh;
        
        //  Case A

        //  First make a unique list of (w1+w2-w3-w4)/2 
        
        keepIndexes = new int[origNbr];
        
        for(i=0;i<origNbr;i++)
	    {
		    keepIndexes[i]=i;
	    } 
        
        for(i=0;i<origNbr;i++)
        {
            for(j=i+1;j<origNbr;j++)
            {
                if( abs((listQh[i*nbrQh]+listQh[i*nbrQh+1]-listQh[i*nbrQh+2]-listQh[i*nbrQh+3])/2.0 - (listQh[j*nbrQh]+listQh[j*nbrQh+1]-listQh[j*nbrQh+2]-listQh[j*nbrQh+3])/2.0) < sameTol )
                {
                    keepIndexes[j]=-1;
                }
            }
        }
        
        counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)	counter++;
	    }

	    m_nbrW1W2W3W4A=counter;
        
	    //////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_
        std::cout<<"\treduced number of (w1+w2-w3-w4)/2 is "<<m_nbrW1W2W3W4A<<std::endl;
		#endif
        //////////////////////////////////////////////////////////////////////////////////
        
        m_uniqueW1W2W3W4A = new (std::nothrow) dcmplx[m_nbrW1W2W3W4A];

	    counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)
		    {
			    m_uniqueW1W2W3W4A[counter]=(listQh[i*nbrQh]+listQh[i*nbrQh+1]-listQh[i*nbrQh+2]-listQh[i*nbrQh+3])/2.0;
			    
			    //////////////////////////////////////////////////////////////////////////////////
			    #if _TEST_MODE_
			    std::cout<<"\t filtered (w1+w2-w3-w4)/2 "<<counter<<" is "<<m_uniqueW1W2W3W4A[counter]<<std::endl;
				#endif
			    
			    counter++;
		    }
	    }
        
        delete[] keepIndexes;
        
        //  Generate theta function objects with the corresponding offsets
        //  The first m_nbrW1W2W3W4A are for the Theta_2 functions, the next m_nbrW1W2W3W4A are for Theta_3
        //  and the third set are for Theta_4 functions

	    m_thetaFuncsW1W2W3W4A=new (std::nothrow) utilities::ThetaLookUp*[m_wfData->torusDegeneracy*m_nbrW1W2W3W4A];

        counter=0;

	    for(i=0;i<m_nbrW1W2W3W4A;i++)    //  Theta_2 objects
	    {
		    m_thetaFuncsW1W2W3W4A[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,0.0,-m_uniqueW1W2W3W4A[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2W3W4A;i++)    //  Theta_3 objects
	    {
		    m_thetaFuncsW1W2W3W4A[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.0,-m_uniqueW1W2W3W4A[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2W3W4A;i++)    //  Theta_4 objects
	    {
		    m_thetaFuncsW1W2W3W4A[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.5,-m_uniqueW1W2W3W4A[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    //  Case B

        //  First make a unique list of (w1-w2+w3-w4)/2 
        
        keepIndexes = new int[origNbr];
        
        for(i=0;i<origNbr;i++)
	    {
		    keepIndexes[i]=i;
	    } 
        
        for(i=0;i<origNbr;i++)
        {
            for(j=i+1;j<origNbr;j++)
            {
                if( abs((listQh[i*nbrQh]-listQh[i*nbrQh+1]+listQh[i*nbrQh+2]-listQh[i*nbrQh+3])/2.0 - (listQh[j*nbrQh]-listQh[j*nbrQh+1]+listQh[j*nbrQh+2]-listQh[j*nbrQh+3])/2.0) < sameTol )
                {
                    keepIndexes[j]=-1;
                }
            }
        }
        
        counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)	counter++;
	    }

	    m_nbrW1W2W3W4B=counter;
        
	    //////////////////////////////////////////////////////////////////////////////////
	    #if _TEST_MODE_
        std::cout<<"\treduced number of (w1-w2+w3-w4)/2 is "<<m_nbrW1W2W3W4B<<std::endl;
		#endif
        //////////////////////////////////////////////////////////////////////////////////

        m_uniqueW1W2W3W4B = new (std::nothrow) dcmplx[m_nbrW1W2W3W4B];

	    counter=0;

	    for(i=0;i<origNbr;i++)
	    {
		    if(keepIndexes[i]!=-1)
		    {
			    m_uniqueW1W2W3W4B[counter]=(listQh[i*nbrQh]-listQh[i*nbrQh+1]+listQh[i*nbrQh+2]-listQh[i*nbrQh+3])/2.0;
			    
			    //////////////////////////////////////////////////////////////////////////////////
			    #if _TEST_MODE_
			    std::cout<<"\t filtered (w1-w2+w3-w4)/2 "<<counter<<" is "<<m_uniqueW1W2W3W4B[counter]<<std::endl;
				#endif
			    //////////////////////////////////////////////////////////////////////////////////
			    
			    counter++;
		    }
	    }

        delete[] keepIndexes;
        
        //  Generate theta function objects with the corresponding offsets
        //  The first m_nbrW1W2W3W4B are for the Theta_2 functions, the next m_nbrW1W2W3W4B are for Theta_3
        //  and the third set are for Theta_4 functions

	    m_thetaFuncsW1W2W3W4B=new (std::nothrow) utilities::ThetaLookUp*[m_wfData->torusDegeneracy*m_nbrW1W2W3W4B];

        counter=0;

	    for(i=0;i<m_nbrW1W2W3W4B;i++)    //  Theta_2 objects
	    {
		    m_thetaFuncsW1W2W3W4B[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,0.0,-m_uniqueW1W2W3W4B[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2W3W4B;i++)    //  Theta_3 objects
	    {
		    m_thetaFuncsW1W2W3W4B[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.0,-m_uniqueW1W2W3W4B[i]-smallShift,false,1);
		    
		    counter++;
	    }
	    
	    for(i=0;i<m_nbrW1W2W3W4B;i++)    //  Theta_4 objects
	    {
		    m_thetaFuncsW1W2W3W4B[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.0,0.5,-m_uniqueW1W2W3W4B[i]-smallShift,false,1);
		    
		    counter++;
	    }
    }
    
    //////////////////////////////////////////////////////////////////////////////////
	//  Tabulate COM theta functions (for calculation of Theta_1(zcm - com0)
	//  Where zcm = sum_i z_i + (w1+w2+...)/2.0
	
	int origNbr = (int)nbrPos/nbrQh;
	int k;
	dcmplx tot1,tot2;
		
	//  First make a unique list of (w1+w2+...) values
	
	keepIndexes = new int[origNbr];
	
	for(i=0;i<origNbr;i++)
	{
		keepIndexes[i]=i;
	}
	
	for(i=0;i<origNbr;i++)
	{
		tot1=0.0;
		
		for(k=0;k<nbrQh;k++)
		{
			tot1+=listQh[i*nbrQh+k];
		}
	
		for(j=i+1;j<origNbr;j++)
		{
			tot2=0.0;
					
			for(k=0;k<nbrQh;k++)
			{
				tot2+=listQh[j*nbrQh+k];
			}

			if( abs(tot1-tot2) < sameTol )
			{
				keepIndexes[j]=-1;
			}
		}
	}
	
	counter=0;

	for(i=0;i<origNbr;i++)
	{

		if(keepIndexes[i]!=-1)	counter++;
	}

	m_nbrCom=counter;
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	std::cout<<"\treduced number of (w1+w2+...)/2.0 is "<<m_nbrCom<<std::endl;
	#endif
	//////////////////////////////////////////////////////////////////////////////////
	
	m_uniqueCom = new (std::nothrow) dcmplx[m_nbrCom];

	counter=0;

	for(i=0;i<origNbr;i++)
	{
		if(keepIndexes[i]!=-1)
		{
			tot1=0.0;
		
			for(k=0;k<nbrQh;k++)
			{
				tot1+=listQh[i*nbrQh+k];
			}
		
			m_uniqueCom[counter]=tot1/2.0;
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_
			std::cout<<"\t filtered (w1+w2+...)/2.0 "<<counter<<" is "<<m_uniqueCom[counter]<<std::endl;
			#endif
			//////////////////////////////////////////////////////////////////////////////////
			
			counter++;
		}
	}
	
	delete[] keepIndexes;
	
	//  Generate theta function objects with the corresponding offsets

	m_thetaFuncsCom = new (std::nothrow) utilities::ThetaLookUp*[m_nbrCom];

	counter=0;

	for(i=0;i<m_nbrCom;i++)
	{
		m_thetaFuncsCom[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,-0.5,-m_uniqueCom[i]+mooreReadComZero-smallShift,true,1);
		counter++;
	}
	
	return;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_TORUS_GEOMETRY_

//////////////////////////////////////////////////////////////////////////////////
//! \brief 	Moore--Read wave function algorithm for a lattice system, in the torus
//!	geometry. This is just the ground state. 
//! 
//!  The look-up table (thetaLookUpTable) should contain a 1D array such that
//!  the log of the theta function which would be evaluated for a given
//!  co-ordinate separation dx and dy lattice spacings 
//!	
//!	 See PRB 54, 16864 (1996) for more details about these wave functions.
//!  See also PRB 77, 125321 (2008).
//!	
//!	Note: the Pfaffian matrix is given by Pf[Theta_a(z_i-z_j)/Theta_1(z_i-z_j)]
//!
//!	In some cases the numerators vanish, making taking the Pfaffian numerically
//!	tricky. Nevertheless it should be OK.
//! 
//////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::MooreRead::EvaluateWfTorus(
	const int n,					    //!<	Number of particles to perform 
									    //!		the evaluation for
	std::complex<int>* latticeConfig,	//!<	An array of complex integers representing 
										//!		the electron configuration on the lattice
	const int whichState)			    //!<	An integer denoting which of the degenerate 
										//!		torus states to evaluate for
	const
{
    //	Local variables
	
    int i,j;
    double fluxFactor;
    dcmplx tmpDcmplx,zcm;
    std::complex<int> tmpIntCmplx,diff;
    int Lx=m_wfData->Lx,Ly=m_wfData->Ly;
    int diffx,diffy;
	dcmplx *p_row;
	dcmplx wfValue;

    //////////////////////////////////////////////////////////////////////////////////
    //	Print co-ordinates such that they can be imported into Mathematica
    #if _TEST_MODE_
    
    std::cout.precision(15);
    
    std::cout<<"\n\t-----------START EVALUATION-----------"<<std::endl;
    
    std::cout<<"\n\tElectron co-ordinates:"<<std::endl;

	std::cout<<"\t{";
	for (i=0;i<n;i++)
	{
		std::cout<<real(latticeConfig[i])<<"+I*"<<imag(latticeConfig[i]);
		if(i<n-1)
		{
			std::cout<<",";
		}
	}
	std::cout<<"}"<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//	Error and warning messages	

    //	Three wave function types are available in the pfaffian term

    if(whichState<0 || whichState>2)		
    {
        std::cout<<"\nERROR: MooreRead::evaluate_wf(int n, dcmplx *z,int whichState) whichState argument must be between 0 and 2"<<std::endl;
        exit(-1);
    }

	if(m_wfData->jastrowExponent>1)
	{
		std::cout<<"WARNING: MooreRead::evaluate_wf() is only currently programmed with whichState nu=1 COM term"<<std::endl;
		getchar();
	}

	wfValue=dcmplx(0,0);
	
	//////////////////////////////////////////////////////////////////////////////////
	//  Implement Jastrow factor part
	
	for (i=0;i<n;i++)
	{
		tmpIntCmplx=latticeConfig[i];
		
		for (j=i+1;j<n;j++)
		{
			diff=tmpIntCmplx-latticeConfig[j];
			
			//	Map the domain [-Lx to Lx ] onto the domain [0 to 2 Lx] and similarly for Ly
			
			diffx=real(diff)+Lx;
			diffy=imag(diff)+Ly;
			
			wfValue+=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			
				m_testValA=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			
			#endif
			//////////////////////////////////////////////////////////////////////////////////

			/*
			
			if(diffx>=0 && diffy>=0)
			{
				wfValue+=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx<0 && diffy>=0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly+diffy);
				wfValue+=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx>=0 && diffy<0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly-diffy);
				wfValue+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx<0 && diffy<0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly-diffy);
				wfValue+=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			
			*/
			
			//////////////////////////////////////////////////////////////////////////////////
			//  Check all values 
			#if _TEST_MODE_ == 2
			
			std::cout<<"\tCheck z"<<i<<"-z"<<j<<" table value (should match)"<<std::endl;
			std::cout<<"\t"<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
			(dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
			
			#endif
			//////////////////////////////////////////////////////////////////////////////////
		}
	}
	
	wfValue*=(double)m_wfData->jastrowExponent;
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
		std::cout<<"\n\tWave function value after Jastrow term: "<<wfValue<<std::endl;
		
		m_prevWfValue=wfValue;
		
		std::cout<<"\n\tCalculating Pfaffian term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//  Implement Pfaffian part

    //////////////////////////////////////////////////////////////////////////////////
	//	Set upper triangular parts to be theta_a(z[i]-z[j])/theta_1(z[i]-z[j])

	for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
	{
		tmpIntCmplx=latticeConfig[i];

		for(j=i+1;j<n;j++)
		{
			diff=tmpIntCmplx-latticeConfig[j];

			//	Map the domain [-Lx to Lx ] onto the domain [0 to 2 Lx] and similarly for Ly
						
			diffx=real(diff)+Lx;
			diffy=imag(diff)+Ly;

			//	There are 4 cases to consider, depending on the signs of diffx and diffy
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			
			    std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<" table value (should match)"<<std::endl;
			
			#endif
            //////////////////////////////////////////////////////////////////////////////////

			*(p_row+j)=*(m_thetaFuncs[whichState+1]->thetaLookUpTable+diffx*(2*Ly)+diffy) - *(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			
				m_testValA=*(m_thetaFuncs[whichState+1]->thetaLookUpTable+diffx*Ly+diffy);
				m_testValB=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
			
			#endif
			//////////////////////////////////////////////////////////////////////////////////
			    
			/*    
			    
			if(diffx>=0 && diffy>=0)
			{
				*(p_row+j)=*(m_thetaFuncs[whichState+1]->thetaLookUpTable+diffx*Ly+diffy) - *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=*(m_thetaFuncs[whichState+1]->thetaLookUpTable+diffx*Ly+diffy);
					m_testValB=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////

			}
			else if(diffx<0 && diffy>=0)
			{
				*(p_row+j)=0;

				//	numerator theta_a
				tmpDcmplx=*(m_thetaFuncs[whichState+1]->thetaLookUpTable-diffx*Ly+diffy);
				*(p_row+j)+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				#endif
	            //////////////////////////////////////////////////////////////////////////////////

				//	denominator theta_1
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly+diffy);
				*(p_row+j)-=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValB=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				#endif
	            //////////////////////////////////////////////////////////////////////////////////

			}
			else if(diffx>=0 && diffy<0)
			{
				*(p_row+j)=0;

				//	numerator theta_a
				tmpDcmplx=*(m_thetaFuncs[whichState+1]->thetaLookUpTable+diffx*Ly-diffy);
				*(p_row+j)+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				    
				#endif
	            //////////////////////////////////////////////////////////////////////////////////

				//	denominator theta_1
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly-diffy);
				*(p_row+j)-=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValB=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				#endif
	            //////////////////////////////////////////////////////////////////////////////////

			}
			else if(diffx<0 && diffy<0)
			{
				*(p_row+j)=0;

				//	numerator theta_a
				tmpDcmplx=*(m_thetaFuncs[whichState+1]->thetaLookUpTable-diffx*Ly-diffy);
				*(p_row+j)+=dcmplx(real(tmpDcmplx),imag(tmpDcmplx));

                //////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValA=dcmplx(real(tmpDcmplx),imag(tmpDcmplx));
				
				#endif
	            //////////////////////////////////////////////////////////////////////////////////

				//	denominator theta_1
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly-diffy);
				*(p_row+j)-=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
				    m_testValB=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				#endif
	            //////////////////////////////////////////////////////////////////////////////////
			}
			
			*/
			
			*(p_row+j)=exp(*(p_row+j));
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			
			switch(whichState)
			{
			    case 0:
			    
			        std::cout<<"\n\tNumerator: "<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
	                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx);
			    
			        break;
			        
			    case 1:
			    
			        std::cout<<"\n\tNumerator: "<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
	                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx);
			    
			        break;
			        
			    case 2:
			    
			        std::cout<<"\n\tNumerator: "<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
	                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx);
			    
			        break;
			}        
			
			std::cout<<"\n\tDenominator: "<<m_testValB<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
	        (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl; 
			   
			#endif
            //////////////////////////////////////////////////////////////////////////////////
		}
	}

    //////////////////////////////////////////////////////////////////////////////////
    //  Print Pfaffian before processing
    #if _TEST_MODE_ == 2

    PrintPfaffian(n);
    
    #endif
	//////////////////////////////////////////////////////////////////////////////////

    wfValue+=utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,n);

    //////////////////////////////////////////////////////////////////////////////////
    //  Print Pfaffian after processing
    #if _TEST_MODE_ == 2

    PrintPfaffian(n);
    
    #endif
	//////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Pfaffian part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating Gaussian term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
    //	Implement Gaussian exponential part

	fluxFactor=(PI/2.0)*(double)(m_wfData->flux)/(Lx*Ly);

	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];

		wfValue-=fluxFactor*abs(tmpDcmplx)*abs(tmpDcmplx);
		wfValue+=fluxFactor*tmpDcmplx*tmpDcmplx;
	}

    //////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Gaussian part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating COM term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//  Implement centre of mass part

	zcm=0.0;

	for(i=0;i<n;i++)
	{
		zcm+=latticeConfig[i];
	}

	//	In the nu=1 case there is only 1, associated with the nu=1 Jastrow factor

	wfValue+=utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(zcm-mooreReadComZero+smallShift)/(double)Lx,(double)Lx/Ly);

    //////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after COM term: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\t----------COMPLETED EVALUATION----------\n\n"<<std::endl;
	    getchar();
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_TORUS_GEOMETRY_

//////////////////////////////////////////////////////////////////////////////////
//! \brief 	Moore--Read wave function algorithm for a lattice system, in the torus
//!	geometry and with some number of quasi-holes.
//!
//!	This function consults various look up tables for efficient evaluation
//!	of theta functions. In order to use the look up tables you must call
//!	LatticeDataPrecalculation() first. This function will still be evaluated 
//!	even if the look up tables are not provided. 
//! 
//////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::MooreRead::EvaluateQuasiholeWfTorus(
	const int n,						        //!<	Number of particles to perform 
												//!		the evaluation for
	std::complex<int>* latticeConfig,		    //!<	An array of complex integers representing
												//!	    the electron configuration on the lattice
	const int whichState,				        //!<	An integer to specify which degenerate 
												//!		torus state we're in:
												//!		- 0,2,4 is fusion channel A, 
												//!		  (for the 3 degenerate states)	
												//!		- 1,2,5 is the same for fusion channel B.
	const int nbrQh,					        //!<	Number of quasi holes (NOTE the quasi holes are not 
												//!		necessarily confined to the lattice
	dcmplx* zQh)							    //!<	An array of quasi hole coordinates
	const
{
	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

		m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////

	//	Local variables
	
	int i,j;
	int flux=m_wfData->flux+nbrQh;
	double fluxFactor;
	dcmplx tmpDcmplx,zcm;
	std::complex<int> tmpIntCmplx,diff;
	int Lx=m_wfData->Lx,Ly=m_wfData->Ly;
	int diffx,diffy;
	dcmplx *p_row;	
	dcmplx wfValue;

    //////////////////////////////////////////////////////////////////////////////////
    //	Print co-ordinates such that they can be imported into Mathematica
    #if _TEST_MODE_
	
	std::cout.precision(15);

    std::cout<<"\n\t-----------START EVALUATION-----------"<<std::endl;
    
    std::cout<<"\n\tElectron co-ordinates:"<<std::endl;

	std::cout<<"\t{";
	for (i=0;i<n;i++)
	{
		std::cout<<real(latticeConfig[i])<<"+I*"<<imag(latticeConfig[i]);
		if(i<n-1)
		{
			std::cout<<",";
		}
	}
	std::cout<<"}"<<std::endl;
	
	std::cout<<"\n\tQuasi-hole co-ordinates:"<<std::endl;
	
	std::cout<<"\t{";
	for (i=0;i<nbrQh;i++)
	{
		std::cout<<real(zQh[i])<<"+I*"<<imag(zQh[i]);
		if(i<nbrQh-1)
		{
			std::cout<<",";
		}
	}
	std::cout<<"}"<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//	Error and warning messages	
    
    if(nbrQh==2)       //	Three wave function types are available in the Pfaffian term
    {
        if(whichState<0 || whichState>2)		
        {
            std::cout<<"ERROR: MooreRead::evaluate_wf(int n, dcmplx *z,int whichState) whichState argument must be between 0 and 2"<<std::endl;
            exit(-1);
        }
    }
    
    if(nbrQh==4)       //	Four wave function types are available in the Pfaffian term
    {
        if(whichState<0 || whichState>5)		
        {
            std::cout<<"ERROR: MooreRead::evaluate_wf(int n, dcmplx *z,int whichState) whichState argument must be between 0 and 5"<<std::endl;
            exit(-1);
        }
    }

	if(m_wfData->jastrowExponent>1)
	{
		std::cout<<"\nWARNING: MooreRead::evaluate_wf() is only currently programmed with whichState nu=1 COM term"<<std::endl;
	}

	wfValue=dcmplx(0,0);

	//////////////////////////////////////////////////////////////////////////////////
	//  Implement Jastrow factor part
	
	for (i=0;i<n;i++)
	{
		tmpIntCmplx=latticeConfig[i];
		
		for (j=i+1;j<n;j++)
		{
			diff=tmpIntCmplx-latticeConfig[j];
			
			//	Map the domain [-Lx to Lx ] onto the domain [0 to 2 Lx] and similarly for Ly
			
			diffx=real(diff)+Lx;
			diffy=imag(diff)+Ly;
			
			wfValue+=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			
				m_testValA=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			
			#endif
			//////////////////////////////////////////////////////////////////////////////////

			/*
			
			if(diffx>=0 && diffy>=0)
			{
				wfValue+=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx<0 && diffy>=0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly+diffy);
				wfValue+=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx>=0 && diffy<0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly-diffy);
				wfValue+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			else if(diffx<0 && diffy<0)
			{
				tmpDcmplx=*(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly-diffy);
				wfValue+=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
				
					m_testValA=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
				
				#endif
				//////////////////////////////////////////////////////////////////////////////////
			}
			
			*/
			
			//////////////////////////////////////////////////////////////////////////////////
			//  Check all values 
			#if _TEST_MODE_ == 2
			
			std::cout<<"\tCheck z"<<i<<"-z"<<j<<" table value (should match)"<<std::endl;
			std::cout<<"\t"<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
			(dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
			
			#endif
			//////////////////////////////////////////////////////////////////////////////////
		}
	}
	
	wfValue*=(double)m_wfData->jastrowExponent;
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
		std::cout<<"\n\tWave function value after Jastrow term: "<<wfValue<<std::endl;
		
		m_prevWfValue=wfValue;
		
		std::cout<<"\n\tCalculating Pfaffian term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

		m_timeJastrow+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

		m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//  Implement Pfaffian part
	//	This part depends on the quasi hole number and position, as well as the number
	//  of quasi-hole in the system
	
	switch(nbrQh)
	{
	    case 2:
	    {
	        int currW1=0,currW2=0,currW1W2=0;
	        dcmplx qhDiff;
	        dcmplx termIJ,termJI;
	        dcmplx termW1,termW2;
	        dcmplx denominator;
	
	        if(m_uniqueQhPositions!=0)		
	        {
		        //	determine which quasi hole positions correspond to zQh[0] and zQh[1]
		        //  assuming that the pre-calculation has been carried out, otherwise ignore

		        for(int k=0;k<m_nbrQhPositions;k++)
		        {
			        if(m_uniqueQhPositions[k]==zQh[0])
			        {
				        currW1=k;
			        }
		        }
		        
		        for(int k=0;k<m_nbrQhPositions;k++)
		        {
			        if(m_uniqueQhPositions[k]==zQh[1])
			        {
				        currW2=k;
			        }
		        }
		        
		        //////////////////////////////////////////////////////////////////////////////////
                #if _TEST_MODE_ == 2

                   std::cout<<"\n\tLocation of quasi-hole(s) in list "<<currW1<<" and "<<currW2<<std::endl;

                #endif
                //////////////////////////////////////////////////////////////////////////////////
		        
		        //  determine also which value of (w1-w2)/2 we currently have
		        
		        for(int k=0;k<m_nbrW1W2;k++)
		        {
			        if( abs(m_uniqueW1W2[k]-(zQh[0]-zQh[1])/2.0) < sameTol )
			        {
				        currW1W2=k;
			        }
		        }
		        //////////////////////////////////////////////////////////////////////////////////
                #if _TEST_MODE_ == 2

                   std::cout<<"\n\tLocation of (w1-w2)/2 in list "<<currW1W2<<std::endl;

                #endif
                //////////////////////////////////////////////////////////////////////////////////
		        
	        }
	        
	        //////////////////////////////////////////////////////////////////////////////////
	        //	Set upper triangular parts to be:
		    //
		    //	theta_a(z[i]-z[j]+(w1-w2)/2)*theta_1(z[i]-w1)*theta_1(z[j]-w2) + theta_a(z[i]-z[j]+(w2-w1)/2) theta_1(z[i]-w2) theta_1(z[j]-w1)
		    //	/ theta_1(z[i]-z[j])
	    
	        if(m_uniqueQhPositions!=0)    //  use pre-calculated values
	        {
	            for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
	            {
		            tmpIntCmplx=latticeConfig[i];
		            
		            //  termW1 = theta_1(z[i]-w1)
		            
		            termW1 = *(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]));
		            
		            //  termW2 = theta_1(z[i]-w2)
		            
		            termW2 = *(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]));
		            
		            //////////////////////////////////////////////////////////////////////////////////
                    #if _TEST_MODE_ == 2
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w1 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<termW1<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
                    (dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w2 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<termW2<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
                    (dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    #endif
                    //////////////////////////////////////////////////////////////////////////////////

		            for(j=i+1;j<n;j++)
		            {
				        diff=tmpIntCmplx-latticeConfig[j];

                        //  termIJ = theta_a(z[i]-z[j]+(w1-w2)/2)

				        diffx=real(diff)+Lx;            //  map the index [-Lx,Lx] to [0,2Lx]
				        diffy=imag(diff)+Ly;            //  map the index [-Ly,Ly] to [0,2Ly]
                        
				        termIJ = *(m_thetaFuncsW1W2[whichState*m_nbrW1W2+currW1W2]->thetaLookUpTable+diffx*(2*Ly)+diffy);   
				       
				        //  termJI = theta_a(z[j]-z[i]+(w1-w2)/2)
				       
				        diffx=-real(diff)+Lx;            //  map the index [-Lx,Lx] to [0,2Lx]
				        diffy=-imag(diff)+Ly;            //  map the index [-Ly,Ly] to [0,2Ly]
				        
                        termJI = *(m_thetaFuncsW1W2[whichState*m_nbrW1W2+currW1W2]->thetaLookUpTable+diffx*(2*Ly)+diffy); 
                        
                        //////////////////////////////////////////////////////////////////////////////////
                        #if _TEST_MODE_ == 2
                        
                        switch(whichState)
                        {    
                            case 0:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"argument of theta: "<<(dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                            
                            case 1:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                            
                            case 2:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]-zQh[1])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;   
                        }
                        
                        #endif
                        //////////////////////////////////////////////////////////////////////////////////
				      
				        //  denominator = theta_1(z[i]-z[j])

				        //	Map the domain [-Lx to Lx ] onto the domain [0 to 2 Lx] and similarly for Ly
				        			
				        diffx=real(diff)+Lx;
				        diffy=imag(diff)+Ly;
				      
				        denominator = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*(2*Ly)+diffy);
				        
				        /*
				        
                        if(diffx>=0 && diffy>=0)
                        {
                            denominator = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
                        }
                        else if(diffx<0 && diffy>=0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly+diffy);
                            denominator = dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
                        }
                        else if(diffx>=0 && diffy<0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly-diffy);
                            denominator = dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
                        }
                        else if(diffx<0 && diffy<0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly-diffy);
                            denominator = dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
                        }
                        
                        */
                        
                        //////////////////////////////////////////////////////////////////////////////////
						#if _TEST_MODE_ == 2
						
						std::cout<<"\n\tCheck z"<<j<<"-w1 table value (should match)"<<std::endl;
						std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
						(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
						
						std::cout<<"\n\tCheck z"<<j<<"-w2 table value (should match)"<<std::endl;
						std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
						(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;

						#endif
						//////////////////////////////////////////////////////////////////////////////////
	
                        //  combine all of the terms
				        
				        tmpDcmplx = 
				        //  theta_a(z[i]-z[j]+(w1-w2)/2)*theta_1(z[i]-w1)*theta_1(z[j]-w2)
                        exp(termIJ+termW1+*(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j])))
                        +
                        //  theta_a(z[j]-z[i]+(w1-w2)/2) theta_1(z[i]-w2) theta_1(z[j]-w1)
                        exp(termJI+termW2+*(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j])))
                        ;

                        *(p_row+j) = exp(log(tmpDcmplx)  - denominator);
				        
                        //////////////////////////////////////////////////////////////////////////////////
                        #if _TEST_MODE_ == 2
                        
                        std::cout<<"\tCombined numerator "<<i<<" "<<j<<" "<<log(tmpDcmplx)<<std::endl;
                        std::cout<<"\tCombined denominator "<<i<<" "<<j<<" "<<denominator<<std::endl;
                        
						#endif
                        ///////////////////////////////////////////////////////////////////////
                        
				    }
				}    
	        }
	        else                        //  calculate directly from theta functions
	        {
	            qhDiff=(zQh[0]-zQh[1])/2.0;     //  gives (w1-w2)/2
	        
	            for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
	            {
		            for(j=i+1;j<n;j++)
		            {
			            diff = dcmplx(real(latticeConfig[i]),imag(latticeConfig[i])) - dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]));
			        
				        switch(whichState)
				        {
				            case 0:
				                
				                //  whichState=0 corresponds to Theta_2
				                termIJ = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				                termJI = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				                
				                break;
				            
				            case 1:
				            
				                //  whichState=1 corresponds to Theta_3
				                termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				                termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				            
				                break;
				            
				            case 2:
				            
				                //  whichState=2 corresponds to Theta_4
				                termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				                termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( -dcmplx(real(diff),imag(diff)) + qhDiff+smallShift)/(double)Lx,(double)Ly/Lx);
				            
				                break;
				        }
				        
				        
                        tmpDcmplx = 
                        exp(termIJ 
                        // theta_1 (z_i - w_1)  
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_2)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)
                        )+
                        exp(termJI
                        // theta_1 (z_i - w_2)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_1)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)
                        );
                        
                        *(p_row+j) = exp(log(tmpDcmplx)  - utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(diff),imag(diff))+smallShift)/(double)Lx,(double)Ly/Lx));
				    }
				}  
	        }
	    
	        break;
	    }
	    
	    case 4:
	    {
	        int k;
	        int b;
	        int currW1=0,currW2=0,currW3=0,currW4=0,currW1W2W3W4A=0,currW1W2W3W4B=0;
	        dcmplx qhDiffA,qhDiffB;
	        dcmplx termIJ,termJI;
	        dcmplx termW1W2,termW3W4;
	        dcmplx denominator;

	        if(m_uniqueQhPositions!=0)		
	        {
		        //	Determine which quasi hole positions correspond to zQh[0],...,zQh[3]
		        //  assuming that the pre-calculation has been carried out, otherwise ignore

		        for(k=0;k<m_nbrQhPositions;k++)
		        {
			        if(m_uniqueQhPositions[k]==zQh[0])
			        {
				        currW1=k;
			        }
		        }
		        
		        for(k=0;k<m_nbrQhPositions;k++)
		        {
			        if(m_uniqueQhPositions[k]==zQh[1])
			        {
				        currW2=k;
			        }
		        }
		        
		        for(k=0;k<m_nbrQhPositions;k++)
                {
                    if(m_uniqueQhPositions[k]==zQh[2])
                    {
                        currW3=k;
                    }
                }
		        
		        for(k=0;k<m_nbrQhPositions;k++)
		        {
			        if(m_uniqueQhPositions[k]==zQh[3])
			        {
				        currW4=k;
			        }
		        }
		        
		        //////////////////////////////////////////////////////////////////////////////////
                #if _TEST_MODE_ == 2

                   std::cout<<"\n\tLocation of quasi-hole(s) in list "<<currW1<<","<<currW2<<","<<currW3<<", and "<<currW4<<std::endl;

                #endif
                //////////////////////////////////////////////////////////////////////////////////
		        
		        //  Determine also which value of (W1+W2-W3-W4)/2 we currently have
		        
		        for(k=0;k<m_nbrW1W2W3W4A;k++)
		        {
			        if( abs(m_uniqueW1W2W3W4A[k]-(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0) < sameTol )
			        {
				        currW1W2W3W4A=k;
			        }
		        }
                
                //////////////////////////////////////////////////////////////////////////////////
                #if _TEST_MODE_ == 2
                		        
		            std::cout<<"\tcurrent (W1+W2-W3-W4)/2 "<<currW1W2W3W4A<<std::endl;
		            std::cout<<"\t correspondg to value = "<<(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0<<std::endl;
		            
		        #endif
                //////////////////////////////////////////////////////////////////////////////////    
		        
		        //  Determine also which value of (W1-W2+W3-W4)/2 we currently have
		        
		        for(k=0;k<m_nbrW1W2W3W4B;k++)
		        {
			        if( abs(m_uniqueW1W2W3W4B[k]-(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0)< sameTol )
			        {
				        currW1W2W3W4B=k;
			        }
		        }
		        
		        //////////////////////////////////////////////////////////////////////////////////
                #if _TEST_MODE_ == 2
		        
		            std::cout<<"\tcurrent (W1-W2+W3-W4)/2 "<<currW1W2W3W4B<<std::endl;
		            std::cout<<"\t correspondg to value = "<<(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0<<std::endl;
		            
		        #endif
                //////////////////////////////////////////////////////////////////////////////////
	        }

            //////////////////////////////////////////////////////////////////////////////////
            //	Set upper triangular parts to be:
	        //
	        //  A-type fusion channel
	        //
	        //	theta_a(z_i-z_j+(w1+w2-w3-w4)/2) * theta_1(z_i-w1) * theta_1(z_i-w2) * theta_1(z_j-w3) * theta_1(z_j-w4)
	        //	+theta_a(z_j-z_i+(w1+w2-w3-w4)/2) * theta_1(z_j-w1) * theta_1(z_j-w2) * theta_1(z_i-w3) * theta_1(z_i-w4)
	        //  / theta_1(z_i - z_j)
	        //
	        //  B-type fusion channel
	        //
	        //	theta_a(z_i-z_j+(w1-w2+w3-w4)/2) * theta_1(z_i-w1) * theta_1(z_i-w2) * theta_1(z_j-w3) * theta_1(z_j-w4)
	        //	+theta_a(z_j-z_i+(w1-w2+w3-w4)/2) * theta_1(z_j-w1) * theta_1(z_j-w2) * theta_1(z_i-w3) * theta_1(z_i-w4)
	        //  / theta_1(z_i - z_j)
	        //
        
            if(m_uniqueQhPositions!=0)    //  use pre-calculated values
            {
                for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
                {
	                tmpIntCmplx=latticeConfig[i];
	                
	                //  termW1W2 = theta_1(z[i]-w1) * theta_1(z[i]-w2)
	                
	                termW1W2 = *(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
	                         + *(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]));
	                
	                //  termW3W4 = theta_1(z[i]-w3) * theta_1(z[i]-w4)
	                
	                termW3W4 = *(m_thetaFuncsQh[currW3]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
	                         + *(m_thetaFuncsQh[currW4]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]));

                    //////////////////////////////////////////////////////////////////////////////////
                    #if _TEST_MODE_ == 2
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w1 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
                    <<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w2 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
                    <<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w3 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW3]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
                    <<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[2]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    std::cout<<"\n\tCheck z"<<i<<"-w4 table value (should match)"<<std::endl;
                    std::cout<<"\n\t"<<*(m_thetaFuncsQh[currW4]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))
                    <<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[3]+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                    
                    #endif
                    //////////////////////////////////////////////////////////////////////////////////

	                for(j=i+1;j<n;j++)
	                {
			            diff=tmpIntCmplx-latticeConfig[j];

                        //  type-A fusion : termIJ = theta_a(z[i]-z[j]+(w1+w2-w3-w4)/2)
                        //  type-B fusion : termIJ = theta_a(z[i]-z[j]+(w1-w2+w3-w4)/2)  

			            diffx=real(diff)+Lx;            //  map the index [-Lx,Lx] to [0,2Lx]
			            diffy=imag(diff)+Ly;            //  map the index [-Ly,Ly] to [0,2Ly]
                        
                        if(whichState%2==0)  //  type-A fusion
                        {
                            b=whichState/2;      // degeneracy index
                        
			                termIJ = *(m_thetaFuncsW1W2W3W4A[currW1W2W3W4A+b*m_nbrW1W2W3W4A]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			            }
			            else        //  type-B fusion
			            {
			                 b=(whichState-1)/2;      // degeneracy index
                        
			                 std::cout<<"b="<<b<<std::endl;
			                 
			                 termIJ = *(m_thetaFuncsW1W2W3W4B[currW1W2W3W4B+b*m_nbrW1W2W3W4B]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			            }   
			           
			            //  type-A fusion : termJI = theta_a(z[j]-z[i]+(w1+w2-w3-w4)/2)
                        //  type-B fusion : termJI = theta_a(z[j]-z[i]+(w1-w2+w3-w4)/2) 
			           
			            diffx=-real(diff)+Lx;            //  map the index [-Lx,Lx] to [0,2Lx]
			            diffy=-imag(diff)+Ly;            //  map the index [-Ly,Ly] to [0,2Ly]
			            
                        if(whichState%2==0)  //  type-A fusion
                        {
                            b=whichState/2;      // degeneracy index
                        
			                termJI = *(m_thetaFuncsW1W2W3W4A[currW1W2W3W4A+b*m_nbrW1W2W3W4A]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			            }
			            else        //  type-B fusion
			            {
			                 b=(whichState-1)/2;      // degeneracy index
                        
			                 std::cout<<"b="<<b<<std::endl;
			                 
			                termJI = *(m_thetaFuncsW1W2W3W4B[currW1W2W3W4B+b*m_nbrW1W2W3W4B]->thetaLookUpTable+diffx*(2*Ly)+diffy);
			            }   
			            
			            //////////////////////////////////////////////////////////////////////////////////
                        #if _TEST_MODE_ == 2
                        
                        switch(whichState)
                        {    
                            case 0:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                                
                            case 1:              

                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                            
                            case 2:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                                
                            case 3:              

                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                                
                            case 4:
                            
                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1+w2-w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[j]-latticeConfig[i]))+(zQh[0]+zQh[1]-zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                                
                            case 5:              

                                std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termIJ<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                std::cout<<"\n\tCheck z"<<j<<"-z"<<i<<"+(w1-w2+w3-w4)/2 table value (should match)"<<std::endl;
                                std::cout<<"\n\t"<<termJI<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,
                                (dcmplx(real(latticeConfig[j]-latticeConfig[i]),imag(latticeConfig[i]-latticeConfig[j]))+(zQh[0]-zQh[1]+zQh[2]-zQh[3])/2.0+smallShift)/(double)Lx,(double)Ly/Lx)<<std::endl;
                                
                                break;
                        }
                        
                        #endif
                        //////////////////////////////////////////////////////////////////////////////////
			          
			            //  Denominator = theta_1(z[i]-z[j])
			             
                        //	Map the domain [-Lx to Lx ] onto the domain [0 to 2 Lx] and similarly for Ly
                        				        			
                        diffx=real(diff)+Lx;
                        diffy=imag(diff)+Ly;
			            
			            denominator = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
			          
			            /*
			            
                        if(diffx>=0 && diffy>=0)
                        {
                            denominator = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly+diffy);
                        }
                        else if(diffx<0 && diffy>=0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly+diffy);
                            denominator = dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
                        }
                        else if(diffx>=0 && diffy<0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable+diffx*Ly-diffy);
                            denominator = dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
                        }
                        else if(diffx<0 && diffy<0)
                        {
                            tmpDcmplx   = *(m_thetaFuncs[0]->thetaLookUpTable-diffx*Ly-diffy);
                            denominator = dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
                        }
                        
                        */

                        //  combine all of the terms
			            
			            tmpDcmplx = 
			            //theta_a(z_i-z_j+(w1+w2-w3-w4)/2) * theta_1(z_i-w1) * theta_1(z_i-w2) * theta_1(z_j-w3) * theta_1(z_j-w4)
                        exp(
                        termIJ+termW1W2
                        +*(m_thetaFuncsQh[currW3]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))
                        +*(m_thetaFuncsQh[currW4]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))
                        )+
                        //theta_a(z_j-z_i+(w1+w2-w3-w4)/2) * theta_1(z_j-w1) * theta_1(z_j-w2) * theta_1(z_i-w3) * theta_1(z_i-w4)
                        exp(
                        termJI+termW3W4
                        +*(m_thetaFuncsQh[currW1]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))
                        +*(m_thetaFuncsQh[currW2]->thetaLookUpTable+real(latticeConfig[j])*Ly+imag(latticeConfig[j]))
                        );

			            //////////////////////////////////////////////////////////////////////////////////
						#if _TEST_MODE_ == 2
						
						std::cout<<"Combined numerator "<<i<<" "<<j<<" "<<log(tmpDcmplx)<<std::endl;
						std::cout<<"Combined denominator "<<i<<" "<<j<<" "<<denominator<<std::endl;
						
						#endif
						///////////////////////////////////////////////////////////////////////
			            
                        *(p_row+j) = exp(log(tmpDcmplx)  - denominator);
			            
			        }
			    }    
            }
            else                        //  calculate directly from theta functions
            {
                qhDiffA=(zQh[0]+zQh[1]-zQh[2]-zQh[4])/2.0;     //  gives (w1+w2-w3-w4)/2
                qhDiffB=(zQh[0]-zQh[1]+zQh[2]-zQh[4])/2.0;     //  gives (w1-w2+w3-w4)/2
            
                for(i=0,p_row=m_pfaff;i<n-1;i++,p_row+=n)
                {
	                for(j=i+1;j<n;j++)
	                {
		                diff = dcmplx(real(latticeConfig[i]),imag(latticeConfig[i])) - dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]));
		            
			            switch(whichState)
			            {
			                case 0:
			                    
			                    //  whichState=0 corresponds to fusion channel A, Theta_2
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( -dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                    
			                    break;
			                    
			                case 1:
			                    
			                    //  whichState=1 corresponds to fusion channel B, Theta_2
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.5,( -dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                    
			                    break;
			                
			                case 2:
			                
			                    //  whichState=2 corresponds to fusion channel A, Theta_3
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                
			                    break;
			                    
			                case 3:
			                
			                    //  whichState=3 corresponds to fusion channel B, Theta_3
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.0,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                
			                    break;
			                
			                case 4:
			                
			                    //  whichState=4 corresponds to fusion channel A, Theta_4
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiffA+smallShift)/(double)Lx,(double)Ly/Lx);
			                
			                    break;
			                    
			                case 5:
			                
			                    //  whichState=5 corresponds to fusion channel B, Theta_4
			                    termIJ = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                    termJI = utilities::thetaFunction::GeneralisedJacobi(0.5,0.0,( -dcmplx(real(diff),imag(diff)) + qhDiffB+smallShift)/(double)Lx,(double)Ly/Lx);
			                
			                    break; 
			            }
			            
			            
                        tmpDcmplx = 
                        exp(termIJ 
                        // theta_1 (z_i - w_1)  
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_i - w_2)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_3)  
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[2]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_4)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[3]+smallShift)/(double)Lx,(double)Ly/Lx)
                        )+
                        exp(termJI
                        // theta_1 (z_i - w_3)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[2]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_i - w_4)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[3]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_1)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[0]+smallShift)/(double)Lx,(double)Ly/Lx)
                        // theta_1 (z_j - w_2)
                        +utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[j]),imag(latticeConfig[j]))-zQh[1]+smallShift)/(double)Lx,(double)Ly/Lx)
                        );
                        
                        *(p_row+j) = exp(log(tmpDcmplx)  - utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(diff),imag(diff))+smallShift)/(double)Lx,(double)Ly/Lx));
			        }
			    }  
            }
	    
	        break;
	    }
	    
	    default:
	            std::cout<<"\nERROR: MooreRead::evaluate_wf() called with"<<nbrQh<<"quasi-hole(s). The algorithm currently only handles 2 or 4 quasiholes."<<std::endl;
	            exit(-1);
	        break;
	}

	//  Now evaluate the pfaffian

	//////////////////////////////////////////////////////////////////////////////////
    //  Print pfaffian before processing
    #if _TEST_MODE_ == 2

    PrintPfaffian(n);
    
    #endif
	//////////////////////////////////////////////////////////////////////////////////

    wfValue+=utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,n);

    //////////////////////////////////////////////////////////////////////////////////
    //  Print pfaffian after processing
    #if _TEST_MODE_ == 2

    PrintPfaffian(n);
    
    #endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Pfaffian part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating Gaussian term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

		m_timePfaffian+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

		m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//	Implement Gaussian exponential part

	fluxFactor=(PI/2.0)*(double)(flux)/(Lx*Ly);

	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];

		wfValue-=fluxFactor*abs(tmpDcmplx)*abs(tmpDcmplx);
		wfValue+=fluxFactor*tmpDcmplx*tmpDcmplx;
	}

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

		m_timeGauss+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

		m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Gaussian part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating COM term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    //  Implement centre of mass part

	if(m_thetaFuncsCom!=0)        //  if a table of value is stored, use that
	{
		//  determine current (w1+w2+...)/2.0 value
		
		int currCom=0;
		dcmplx tot1;
		
		for(i=0;i<m_nbrCom;i++)
		{
			tot1=0.0;
			
			for(j=0;j<nbrQh;j++)
			{
				tot1+=zQh[j]/2.0;
			}

			if( abs(m_uniqueCom[i]-tot1) < sameTol)
			{
				currCom=i;
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_ == 2
	
			std::cout<<"\n\tLocation of (w1+w2+...)/2.0 in list "<<currCom<<std::endl;
			
			if(nbrQh==2)
			{
				std::cout<<"\tcheck: "<<m_uniqueCom[currCom]<<" compare with "<<(zQh[0]+zQh[1])/2.0<<std::endl;
			}
			else if(nbrQh==4)
			{
				std::cout<<"\tcheck: "<<m_uniqueCom[currCom]<<" compare with "<<(zQh[0]+zQh[1]+zQh[2]+zQh[3])/2.0<<std::endl;
			}
	
		#endif
		//////////////////////////////////////////////////////////////////////////////////

		dcmplx thetaArgWithQuasiHole;
		dcmplx thetaArgNoQuasiHole;

		zcm=0.0;
						
		for(i=0;i<n;i++)
		{
			zcm+=latticeConfig[i];
		}

		thetaArgNoQuasiHole=(1.0/Lx)*zcm;
		
		for(j=0;j<nbrQh;j++)
		{
			zcm+=zQh[j]/2.0;
		}
		
		//	subtract off COM zero from argument
		
		zcm-=mooreReadComZero;
		
		thetaArgWithQuasiHole=((1.0/Lx)*(zcm+smallShift));
		
		//  Shift sum of z_i such that the theta function has an argument 
		//	within the range 0<sum_i z_i /Lx <1 (the allows us to use a 
		//	look-up table
		
		dcmplx xShiftFactor=0;
		dcmplx yShiftFactor=0;
		double t = (double)Ly/Lx;
		double arg1=0.5;
		double arg2=-0.5;
		
		//////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_ == 2
	
			std::cout<<"\n\tOriginal argument "<<thetaArgWithQuasiHole<<std::endl;
			
			std::cout<<"\tCorrect result: "<<utilities::thetaFunction::GeneralisedJacobi(
			arg1,arg2,thetaArgWithQuasiHole,t)<<std::endl;
	
		#endif
		//////////////////////////////////////////////////////////////////////////////////
		
		while(real(thetaArgNoQuasiHole)>=1.0)
		{
			thetaArgWithQuasiHole-=1.0;
			thetaArgNoQuasiHole-=1.0;
			xShiftFactor+=2.0*PI*I*arg1;
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			//	check that this correction is consistent
				
				std::cout<<"Result with shift: "<<xShiftFactor+utilities::thetaFunction::GeneralisedJacobi(
				arg1,arg2,thetaArgWithQuasiHole,t)<<std::endl;				
		
			#endif
			//////////////////////////////////////////////////////////////////////////////////
			
		}

		while(imag(thetaArgNoQuasiHole)>=t)
		{
			thetaArgWithQuasiHole-=I*t;
			thetaArgNoQuasiHole-=I*t;
			yShiftFactor+=(PI*t-2.0*PI*I*(thetaArgWithQuasiHole+arg2));
			
			//////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_ == 2
			//	check that this correction is consistent
				
				std::cout<<"Result with shift: "<<xShiftFactor+yShiftFactor+utilities::thetaFunction::GeneralisedJacobi(
				arg1,arg2,thetaArgWithQuasiHole,t)<<std::endl;				
		
			#endif
			//////////////////////////////////////////////////////////////////////////////////
		}

		//////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_ == 2
	
			std::cout<<"\n\tValue of (z1+z2+...) in table "<<thetaArgNoQuasiHole*(double)Lx<<std::endl;
	
		#endif
		//////////////////////////////////////////////////////////////////////////////////
		
		//  Pass this value to the lookup table

		wfValue+=*(m_thetaFuncsCom[currCom]->thetaLookUpTable+(int)real(thetaArgNoQuasiHole*(double)Lx)*Ly*m_wfData->jastrowExponent+(int)imag(thetaArgNoQuasiHole*(double)Lx));
		
		//  correct for the shift factors
		
		wfValue+=xShiftFactor+yShiftFactor;
		
		//////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_ == 2

		//for(int bla=0;bla<Lx*Ly;bla++)
		//{
		//	
		//	std::cout<<"\tTable value: "<<bla<<" "<<*(m_thetaFuncsCom[currCom]->thetaLookUpTable+bla)<<std::endl;
		//}

		std::cout<<"\n\tCheck COM term vs table value. "<<std::endl;
		std::cout<<"\tTable value: "<<wfValue-m_prevWfValue<<std::endl;
		std::cout<<"\tCorrect result: "<<utilities::thetaFunction::GeneralisedJacobi(
					arg1,arg2,((1.0/Lx)*(zcm+smallShift)),t)<<std::endl;

		#endif
		//////////////////////////////////////////////////////////////////////////////////
	}
	else                        //  otherwise calculate from scratch
	{
		wfValue+=utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(zcm-mooreReadComZero+smallShift)/(double)Lx,(double)Lx/Ly);
	}

	//////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

		m_timeCM+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

		m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after COM term: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\t----------COMPLETED EVALUATION----------\n\n"<<std::endl;
	    getchar();
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

    return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
