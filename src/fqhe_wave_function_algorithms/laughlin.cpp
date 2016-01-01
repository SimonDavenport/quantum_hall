////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 04/12/2013
//!
//!  \file 
//!		This cpp file contains functions to evaluate the Laughlin wave function
//!     in various geometries, both for the ground state and quasi-hole states
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

#include "laughlin.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the Laughlin class. 
//!
//!	 Assigns basic wave function characteristics. 
//!
//!	 In the torus geometry case, we populate an array of possible electron-electron
//!	 theta functions. By the end of the constructor we have populated the following 
//!	 theta function table:
//!   - m_thetaFuncs To store a table of theta function values on the lattice
//!	    with no offset
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Laughlin::Laughlin(
	FQHE::WaveFunctionData *wavefunctionData) 		//!<	Address of a struct of wave function data
	: FQHE::WaveFunction(wavefunctionData)
{
	//	Generate basic wave function data
    m_wfData->fillNumerator=1;
    m_wfData->fillDenominator=m_wfData->jastrowExponent;
    m_wfData->shift=m_wfData->jastrowExponent;
    m_wfData->wfName="Laughlin";
    m_wfData->statistics = (m_wfData->jastrowExponent%2==0 ? _BOSONS_ : _FERMIONS_);

	m_wfData->fillingFactor =(double)m_wfData->fillNumerator/m_wfData->fillDenominator;
	
	if(m_wfData->geometry==_SPHERE_ || m_wfData->geometry==_DISC_)
	{	
		m_wfData->flux=(double)m_wfData->nbr/(m_wfData->fillingFactor)-m_wfData->shift;
		m_wfData->monopoleStrength = (m_wfData->nbr*m_wfData->fillDenominator/m_wfData->fillNumerator-m_wfData->shift)/2.0;
	}

    //////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_
	
	if(m_wfData->geometry==_TORUS_)
	{
		//	Generate basic data 
		
		m_wfData->flux=(double)m_wfData->nbr/(m_wfData->fillingFactor);
		m_wfData->maxOccupation=m_wfData->jastrowExponent;
		m_wfData->torusDegeneracy=m_wfData->jastrowExponent;
		m_wfData->fusionChannels=1;
		
		m_thetaFuncs = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,-0.5,dcmplx(0.0,0.0),true,1);
		
		//  set pointers to be null
		
		m_thetaFuncsQh=0;
		m_thetaFuncsCom=0;
		m_uniqueQhPositions=0;
		m_uniqueCom=0;
	}

    //	FOR BENCHMARKING
    //////////////////////////////////////////////////////////////////////////////////
    #if	_BENCHMARK_MODE_

    m_timeJastrow=0;
    m_timeGauss=0;
    m_timeQuasihole=0;
    m_timeCM=0;

    #endif
    //////////////////////////////////////////////////////////////////////////////////
		
	#endif
	//////////////////////////////////////////////////////////////////////////////////
    
    //  Now we can generate a file name
	m_wfData->GenerateFileName();

    //  Generate sphere/disc radius value
    m_wfData->InitRadius();

	//  Print a summary
	m_wfData->CheckAndPrint();
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the Laughlin class. 
//!
//!	 Deallocates memory as required
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Laughlin::~Laughlin()
{	
    ////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_TORUS_GEOMETRY_	
	
	if(m_wfData->geometry==_TORUS_)
	{
		
		delete m_thetaFuncs;

		if(m_thetaFuncsQh!=0)	        delete[] m_thetaFuncsQh;
		if(m_thetaFuncsCom!=0)	        delete[] m_thetaFuncsCom;
		if(m_uniqueQhPositions!=0)	    delete[] m_uniqueQhPositions;
		if(m_uniqueCom!=0)	            delete[] m_uniqueCom;
	
	}

    ////////////////////////////////////////////////////////////////////////////////
    #if	_BENCHMARK_MODE_

    std::cout<<"\t\t\t-----BENCHMARKING-----\n\n";
    std::cout<<"\tJASTROW TIME:\t"<<m_timeJastrow<<std::endl;
    std::cout<<"\tGAUSS TIME:\t"<<m_timeGauss<<std::endl;
    std::cout<<"\tQUASIHOLE TIME:\t"<<m_timeQuasihole<<std::endl;
    std::cout<<"\tCM TERM TIME:\t"<<m_timeCM<<std::endl<<std::endl;

    #endif
    ////////////////////////////////////////////////////////////////////////////////

	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Laughlin wave function algorithm in the sphere geometry 
//!
//!	For more information see e.g. PRL 50 1395--1398 (1983).
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::Laughlin::EvaluateWfSphere(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* u,	    //!<	The memory address of the start of an array of spinor u coordinates
	dcmplx* v)	    //!<	The memory address of the start of an array of spinor v coordinates
	const
{
    int i,j;
	dcmplx wfValue;
    
	wfValue=dcmplx(0,0);
	
	for (i=0;i<n-1;i++)
	{
		for (j=i+1;j<n;j++)
		{
			wfValue+=log(*(u+i)**(v+j)-*(u+j)**(v+i));
		}
	}
	
	wfValue*=m_wfData->jastrowExponent;
	
	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////
#if _ENABLE_DISC_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Laughlin wave function algorithm in the disc geometry 
//!
//!	For more information see e.g. PRL 50 1395--1398 (1983).
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::Laughlin::EvaluateWfDisc(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* z)	    //!<	The memory address of the start of an array of disc coordinates
	const
{
    int i,j;
    double tmp;
	dcmplx wfValue;
    
	wfValue=dcmplx(0,0);
	
	for (i=0;i<n;i++)
	{
		for (j=i+1;j<n;j++)
		{
			wfValue+=log(*(z+i)-*(z+j));
		}
	}
	
	wfValue*=m_wfData->jastrowExponent;
	
	for(i=0;i<n;i++)
	{
		//	Implement Gaussian exponential part
		tmp=abs(*(z+i));
		wfValue-=tmp*tmp/4;	
	}
	
	return wfValue;	
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Laughlin wave function algorithm with quasi holes in the sphere geometry 
//!
//!	For more information see e.g. PRL 50 1395--1398 (1983).
//!
////////////////////////////////////////////////////////////////////////////////
dcmplx FQHE::Laughlin::EvaluateQuasiholeWfSphere(
    const int n,		//!<	The number of particles to perform the evaluation for
    dcmplx* u,	        //!<	The memory address of the start of an array of spinor u coordinates
    dcmplx* v,	        //!<	The memory address of the start of an array of spinor v coordinates
    const int nbrQh,    //!<    The number of quasi holes to perform the evaluation for
    dcmplx* uQh,        //!<    The memory address of the start of an array 
                        //!     of quasi hole spinor u coordinates
    dcmplx* vQh)        //!<    The memory address of the start of an array 
                        //!     of quasi hole spinor v coordinates
	const
{
	dcmplx wfValue;

	wfValue=dcmplx(0,0);
    
	for (int i=0;i<n-1;i++)
	{
		for (int j=i+1;j<n;j++)
		{
			wfValue+=log(*(u+i)**(v+j)-*(u+j)**(v+i));
		}
	}
	
	wfValue*=m_wfData->jastrowExponent;
    
    for(int j=0;j<nbrQh;j++)
    {
        for(int i=0;i<n;i++)
        {
			wfValue+=log(*(u+i)**(vQh+j)-*(uQh+j)**(v+i));
		}
	}
	
	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////
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
//!	  -	m_thetaFuncsCom To store a list of tables of theta function values
//!		offset by the centre of mass coordinate (minus the COM zero depending
//!		on whether we use the COM zeros formalism or not)
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Laughlin::LatticeDataPrecalculation(
	const int nbrPos,       //!<	Number of quasi hole positions along the braid
	const int nbrQh,        //!<	Number of quasi holes   
	dcmplx* listQh)	        //!<	The starting address of a list of all quasi-hole
						    //!		positions that are visited along the path 
						    //!		(can include duplicates). The size of this list
						    //!		must equal nbrPos.
{
    //  Local variables

	int i,j,counter;
	int *keepIndexes;
	dcmplx wfValue;

	//	Determine a list of all quasihole positions without repeated values

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

	m_thetaFuncsQh=new (std::nothrow) utilities::ThetaLookUp*[m_nbrQhPositions];

	for(i=0;i<m_nbrQhPositions;i++)
	{
		m_thetaFuncsQh[i] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,0.5,-0.5,m_uniqueQhPositions[i],true,1);
	}

	//////////////////////////////////////////////////////////////////////////////////
	//  Tabulate COM theta functions (using the generalized theta function notation)
	#if _WITH_GENERALISED_THETA_FUNCS_
	
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
    std::cout<<"\reduced number of (w1+w2+...) is "<<m_nbrCom<<std::endl;
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
	    
		    m_uniqueCom[counter]=tot1;
		    
		    //////////////////////////////////////////////////////////////////////////////////
			#if _TEST_MODE_
		    std::cout<<"\t filtered (w1+w2+...) "<<counter<<" is "<<m_uniqueCom[counter]<<std::endl;
			#endif
		    //////////////////////////////////////////////////////////////////////////////////
		    
		    counter++;
	    }
    }
    
    delete[] keepIndexes;
    
    //  Generate theta function objects with the corresponding offsets
    //	Later on we need to evaluate things with arguments up to Lx*jastrowExponent and Ly*jastrowExponent
    //	Hence the last argument is set to m_wfData->jastrowExponent
    //  The first m_nbrCom are a=0, the next m_nbrCom are a=1 etc...

    m_thetaFuncsCom = new (std::nothrow) utilities::ThetaLookUp*[m_nbrCom*m_wfData->torusDegeneracy];

    double thetaArg=(m_wfData->flux+nbrQh-m_wfData->jastrowExponent)/2.0;

    counter=0;

    for(int a=0;a<m_wfData->torusDegeneracy;a++)
    {

        for(i=0;i<m_nbrCom;i++)
        {
            m_thetaFuncsCom[counter] = new utilities::ThetaLookUp(m_wfData->Lx,m_wfData->Ly,thetaArg+(double)a/m_wfData->jastrowExponent,thetaArg,-m_uniqueCom[i],true,m_wfData->jastrowExponent);
            counter++;
        }
    }

	#endif
	//////////////////////////////////////////////////////////////////////////////////

	return;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_TORUS_GEOMETRY_

//////////////////////////////////////////////////////////////////////////////////
//! \brief 	Laughlin wave function algorithm for a lattice system, in the torus
//!	geometry. This is just the ground state. 
//! 
//!  The look-up table (thetaLookUpTable) should contain a 1D array such that
//!  the log of the theta function which would be evaluated for a given
//!  co-ordinate separation dx and dy lattice spacings 
//!	
//!	 See PRB 54, 16864 (1996) for more details about these wave functions.
//!  See also PRB 77, 125321 (2008).
//! 
//////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::Laughlin::EvaluateWfTorus(
	const int n,					    //!<	Number of particles to perform 
									    //!		the evaluation for
	std::complex<int>* latticeConfig,	//!<	An array of complex integers representing 
										//!		the electron configuration on the lattice
	const int whichState)			    //!<	An integer denoting which of the degenerate 
										//!		torus states to evaluate for
	const
{

    //  Local variables

    int i,j;
    double fluxFactor;
    dcmplx tmpDcmplx,zcm;
    std::complex<int> tmpIntCmplx,diff;
    int Lx=m_wfData->Lx,Ly=m_wfData->Ly;
    int diffx,diffy;
	dcmplx wfValue;

    //////////////////////////////////////////////////////////////////////////////////
    //	Print co-ordinates such that they can be imported into mathematica
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

    if(whichState>=m_wfData->jastrowExponent)
    {
        std::cout<<"ERROR: Laughlin::evaluate_wf(int n, dcmplx *z,int whichState) whichState argument must be between 0 and jastrowExponent-1"<<std::endl;
        exit(-1);
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
			
            diffx=real(diff);
            diffy=imag(diff);

            if(diffx>=0 && diffy>=0)
            {
            	wfValue+=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly+diffy);
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly+diffy);
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
            }
            else if(diffx<0 && diffy>=0)
            {
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable-diffx*Ly+diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
            }
            else if(diffx>=0 && diffy<0)
			{
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly-diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
			}
            else if(diffx<0 && diffy<0)
			{
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable-diffx*Ly-diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
			}
			
			//////////////////////////////////////////////////////////////////////////////////
	        //  Check all values 
	        #if _TEST_MODE_ == 2
	        
	        std::cout<<"\n\tCheck z"<<i<<"-z"<<j<<" table value (should match)"<<std::endl;
	        std::cout<<"\n\t"<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
	        (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j])))/(double)Lx,(double)Ly/Lx)<<std::endl;
	        
	        #endif
	        //////////////////////////////////////////////////////////////////////////////////
		}
	}
	
	wfValue*=(double)m_wfData->jastrowExponent;
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Jastrow term: "<<wfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating Gaussian term..."<<std::endl;

	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
    //	Implement Gaussian exponential part

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_GENERALISED_THETA_FUNCS_

	fluxFactor=(PI)*(double)(m_wfData->flux)/(Lx*Ly);
	 
	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];
		
		wfValue-=fluxFactor*imag(tmpDcmplx)*imag(tmpDcmplx);
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_COM_ZEROES

	fluxFactor=(PI/2.0)*(double)(m_wfData->flux)/(Lx*Ly);

	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];

		wfValue-=fluxFactor*abs(tmpDcmplx)*abs(tmpDcmplx);
		wfValue+=fluxFactor*tmpDcmplx*tmpDcmplx;
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////
    
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
	
	//////////////////////////////////////////////////////////////////////////////////
    #if _WITH_COM_ZEROES

	//	Define COM zeros:
	
	dcmplx laughlinComZeros[m_wfData->jastrowExponent];
	
	for(i=0;i<m_wfData->jastrowExponent;i++)
	{
		laughlinComZeros[i]=sqrt(2.0)*exp(2.0*(double)PI*(double)i/m_wfData->jastrowExponent*dcmplx(0,1))+(double)whichState;
	}
	
	//	NOTE: COM ZEROS HAVE ONLY BEEN FIXED CORRECTLY FOR NU=1/2

	laughlinComZeros[0]=1.0/sqrt(2.0);
	laughlinComZeros[1]=-sqrt(2.0);

    //////////////////////////////////////////////////////////////////////////////////
    #if _TEST_MODE_ == 2

        std::cout<<"\n\tCenter-of-mass zeros are: "<<laughlinComZeros[0]<<" "<<laughlinComZeros[1]<<std::endl;

    #endif
    //////////////////////////////////////////////////////////////////////////////////
	
	for(i=0;i<m_wfData->jastrowExponent;i++)
	{
		wfValue+=utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(zcm-laughlinComZeros[i])/(double)Lx,(double)Lx/Ly);
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_GENERALISED_THETA_FUNCS_
	
	double thetaArg=(m_wfData->flux-m_wfData->jastrowExponent)/2.0;
	
	wfValue+=utilities::thetaFunction::GeneralisedJacobi(thetaArg+(double)whichState/m_wfData->jastrowExponent,thetaArg,((double)m_wfData->jastrowExponent/Lx)*zcm,(double)m_wfData->jastrowExponent*Ly/Lx);
	
    #endif	
    //////////////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////////////
//! \brief 	Laughlin wave function algorithm for a lattice system, in the torus
//!	geometry and with some number of quasi-holes.
//!
//!	This function consults various look up tables for efficient evaluation
//!	of theta functions. In order to use the look up tables you must call
//!	LatticeDataPrecalculation() first. This function will still be evaluated 
//!	even if the look up tables are not provided. 
//! 
//////////////////////////////////////////////////////////////////////////////////
dcmplx FQHE::Laughlin::EvaluateQuasiholeWfTorus(
	const int n,						//!<	Number of particles to perform 
									    //!		the evaluation for
	std::complex<int>* latticeConfig,   //!<	An array of complex integers representing
									    //!	    the electron configuration on the lattice
	const int whichState,				//!<	An integer to specify which degenerate 
									    //!		torus state we're in.
	const int nbrQh,					//!<	Number of quasi holes (NOTE the quasi holes are not 
									    //!		necessarily confined to the lattice
	dcmplx* zQh)				        //!<	An array of quasi hole coordinates
	const
{
    ////////////////////////////////////////////////////////////////////////////////
    #if	_BENCHMARK_MODE_

	m_timer=clock();

    #endif
    ////////////////////////////////////////////////////////////////////////////////

	int i,j;
	int flux=m_wfData->flux+nbrQh;
	double fluxFactor;
	dcmplx tmpDcmplx,zcm;
	std::complex<int> tmpIntCmplx,diff;
	int Lx=m_wfData->Lx,Ly=m_wfData->Ly;
	int diffx,diffy;
	dcmplx wfValue;
   
	////////////////////////////////////////////////////////////////////////////////
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
	//  Print error messages
	
	if(whichState>=m_wfData->jastrowExponent)
	{
		std::cout<<"ERROR: Laughlin::evaluate_wf(int n, dcmplx *z,int whichState) whichState argument must be between 0 and jastrowExponent-1"<<std::endl;
		exit(-1);
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
			
            diffx=real(diff);
            diffy=imag(diff);

            if(diffx>=0 && diffy>=0)
            {
            	wfValue+=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly+diffy);
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly+diffy);
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
            }
            else if(diffx<0 && diffy>=0)
            {
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable-diffx*Ly+diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),PI-imag(tmpDcmplx));
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
            }
            else if(diffx>=0 && diffy<0)
			{
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable+diffx*Ly-diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),-imag(tmpDcmplx));
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
			}
            else if(diffx<0 && diffy<0)
			{
            	tmpDcmplx=*(m_thetaFuncs->thetaLookUpTable-diffx*Ly-diffy);
            	wfValue+=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
            	
            	//////////////////////////////////////////////////////////////////////////////////
            	#if _TEST_MODE_ == 2
            	
            	    m_testValA=dcmplx(real(tmpDcmplx),imag(tmpDcmplx)-PI);
            	
            	#endif
            	//////////////////////////////////////////////////////////////////////////////////
			}
			
			//////////////////////////////////////////////////////////////////////////////////
	        //  Check all values 
	        #if _TEST_MODE_ == 2
	        
	        std::cout<<"\tCheck z"<<i<<"-z"<<j<<" table value (should match)"<<std::endl;
	        std::cout<<"\t"<<m_testValA<<"\t"<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,
	        (dcmplx(real(latticeConfig[i]-latticeConfig[j]),imag(latticeConfig[i]-latticeConfig[j])))/(double)Lx,(double)Ly/Lx)<<std::endl;
	        
	        #endif
	        //////////////////////////////////////////////////////////////////////////////////
		}
	}
	
	wfValue*=(double)m_wfData->jastrowExponent;
	
	//////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Jastrow term: "<<wfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating Gaussian term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_    

	m_timeJastrow+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

	m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	//	Implement Gaussian exponential part

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_GENERALISED_THETA_FUNCS_

	fluxFactor=(PI)*(double)(flux)/(Lx*Ly);
	 
	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];
		
		wfValue-=fluxFactor*imag(tmpDcmplx)*imag(tmpDcmplx);
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_COM_ZEROES

	fluxFactor=(PI/2.0)*(double)(flux)/(Lx*Ly);

	for(i=0;i<n;i++)
	{
		tmpDcmplx=latticeConfig[i];

		wfValue-=fluxFactor*abs(tmpDcmplx)*abs(tmpDcmplx);
		wfValue+=fluxFactor*tmpDcmplx*tmpDcmplx;
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////
    
    //////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after Gaussian part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating Quasi-hole term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

	m_timeGauss+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

	m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
    //  Implement quasi-hole part
	
	for(j=0;j<nbrQh;j++)
	{
		int currentPosition=0;

		if(m_uniqueQhPositions!=0)		
		{
			//	determine which quasihole position corresponds to zQh[j]

			for(int k=0;k<m_nbrQhPositions;k++)
			{
				if(m_uniqueQhPositions[k]==zQh[j])
				{
					currentPosition=k;
				}
			}
		}

	    //////////////////////////////////////////////////////////////////////////////////
	    #if _TEST_MODE_ == 2
	
	        std::cout<<"\n\tLocation of quasi-hole in list "<<currentPosition<<std::endl;
	
	    #endif
	    //////////////////////////////////////////////////////////////////////////////////

		if(m_uniqueQhPositions==0)        //	standard evaluation if a pre-calcualted list of positions
									    //	has not been assigned
		{
			for(i=0;i<n;i++)
			{
				wfValue+=utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[j])/(double)Lx,(double)Ly/Lx);
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{

				wfValue+=*(m_thetaFuncsQh[currentPosition]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]));
			
				//////////////////////////////////////////////////////////////////////////////////
				#if _TEST_MODE_ == 2
			
					std::cout<<"\t\tTable lookup value vs. standard evaluation:"<<std::endl;
					std::cout<<"\tTable value: "<<*(m_thetaFuncsQh[currentPosition]->thetaLookUpTable+real(latticeConfig[i])*Ly+imag(latticeConfig[i]))<<std::endl;
					std::cout<<"\tShould be: "<<utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(dcmplx(real(latticeConfig[i]),imag(latticeConfig[i]))-zQh[j])/(double)Lx,(double)Ly/Lx)<<std::endl;
			
				#endif
				//////////////////////////////////////////////////////////////////////////////////

			}
		}

	}

    //////////////////////////////////////////////////////////////////////////////////
	#if _TEST_MODE_
	
	    std::cout<<"\n\tWave function value after quasi-hole part: "<<wfValue<<std::endl;
	    std::cout<<"\tIncreased by: "<<wfValue-m_prevWfValue<<std::endl;
	    
	    m_prevWfValue=wfValue;
	    
	    std::cout<<"\n\tCalculating COM term..."<<std::endl;
	
	#endif
	//////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

	m_timeQuasihole+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

	m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    //  Implement centre of mass part

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_COM_ZEROES

	//	Define COM zeros:
	
	dcmplx laughlinComZeros[m_wfData->jastrowExponent];
	
	for(i=0;i<m_wfData->jastrowExponent;i++)
	{
		laughlinComZeros[i]=sqrt(2.0)*exp(2.0*(double)PI*(double)i/m_wfData->jastrowExponent*dcmplx(0,1))+(double)whichState;
	}

	laughlinComZeros[0]=sqrt(2.0)*pow(dcmplx(0,1),(double)whichState);
	laughlinComZeros[1]=-sqrt(2.0)*pow(dcmplx(0,1),(double)whichState);
	
	//	NOTE: COM ZEROS HAVE ONLY BEEN FIXED CORRECTLY FOR NU=1/2

    //////////////////////////////////////////////////////////////////////////////////
    #if _TEST_MODE_ == 2

        std::cout<<"\n\tCenter-of-mass zeros are: "<<laughlinComZeros[0]<<" "<<laughlinComZeros[1]<<std::endl;

    #endif
    //////////////////////////////////////////////////////////////////////////////////
	
	zcm=0.0;
				
	for(i=0;i<n;i++)
	{
		zcm+=latticeConfig[i];
	}
	
	for(j=0;j<nbrQh;j++)
	{
		zcm+=zQh[j]/(double)m_wfData->jastrowExponent;
	}    

	for(i=0;i<m_wfData->jastrowExponent;i++)
	{
		wfValue+=utilities::thetaFunction::GeneralisedJacobi(0.5,-0.5,(zcm-laughlinComZeros[i])/(double)Lx,(double)Lx/Ly);
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    #if _WITH_GENERALISED_THETA_FUNCS_
	
	if(m_thetaFuncsCom!=0)        //  if a table of value is stored, use that
	{
	    //  determine current (w1+w2+...) value
	    
	    int currCom=0;
	    dcmplx tot1;
	    
	    for(i=0;i<m_nbrCom;i++)
		{
		    tot1=0.0;
            
            for(j=0;j<nbrQh;j++)
            {
                tot1+=zQh[j];
            }

			if( abs(m_uniqueCom[i]-tot1) < sameTol)
			{
				currCom=i;
			}
		}
		
	    //////////////////////////////////////////////////////////////////////////////////
		#if _TEST_MODE_ == 2
	
	        std::cout<<"\n\tLocation of w1+w2+... in list "<<currCom<<std::endl;
	        
	        std::cout<<"\tcheck: "<<m_uniqueCom[currCom]<<" compare with "<<zQh[0]+zQh[1]<<std::endl;
	
	    #endif
	    //////////////////////////////////////////////////////////////////////////////////

        dcmplx thetaArgWithQuasiHole;
        dcmplx thetaArgNoQuasiHole;

        zcm=0.0;
                		
		for(i=0;i<n;i++)
		{
			zcm+=latticeConfig[i];
		}

		thetaArgNoQuasiHole=((double)m_wfData->jastrowExponent/Lx)*zcm;
		
		for(j=0;j<nbrQh;j++)
		{
			zcm+=zQh[j]/(double)m_wfData->jastrowExponent;
		}
        
		thetaArgWithQuasiHole=((double)m_wfData->jastrowExponent/Lx)*zcm;
		
		//  Shift sum of z_i such that the theta function has an argument 
		//	within the range 0<sum_i z_i /Lx <1 (the allows us to use a 
		//	look-up table
		
        dcmplx xShiftFactor=0;
        dcmplx yShiftFactor=0;
        double t = (double)m_wfData->jastrowExponent*Ly/Lx;
        double arg1=(flux-m_wfData->jastrowExponent)/2.0+(double)whichState/m_wfData->jastrowExponent;
        double arg2=(flux-m_wfData->jastrowExponent)/2.0;
	    
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

	    wfValue+=*(m_thetaFuncsCom[whichState*m_nbrCom+currCom]->thetaLookUpTable+(int)real(thetaArgNoQuasiHole*(double)Lx)*Ly*m_wfData->jastrowExponent+(int)imag(thetaArgNoQuasiHole*(double)Lx));

	    //  correct for the shift factors
	    
	    wfValue+=xShiftFactor+yShiftFactor;
	    
	    //////////////////////////////////////////////////////////////////////////////////
	    #if _TEST_MODE_ == 2

        std::cout<<"\n\tCheck COM term vs table value. "<<std::endl;
        std::cout<<"\tTable value: "<<wfValue-m_prevWfValue<<std::endl;
        std::cout<<"\tShould be: "<<utilities::thetaFunction::GeneralisedJacobi(
	    arg1,arg2,((double)m_wfData->jastrowExponent/Lx)*zcm,t)<<std::endl;

        #endif
	    //////////////////////////////////////////////////////////////////////////////////
	}
	else                        //  otherwise calculate from scratch
	{
	    wfValue+=utilities::thetaFunction::GeneralisedJacobi(
	    (flux-m_wfData->jastrowExponent)/2.0+(double)whichState/m_wfData->jastrowExponent,
	    (flux-m_wfData->jastrowExponent)/2.0,
	    ((double)m_wfData->jastrowExponent/Lx)*zcm,
	    (double)m_wfData->jastrowExponent*Ly/Lx);
	}

    #endif
    //////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
	#if	_BENCHMARK_MODE_

	m_timeCM+=1.0/3600.0*( double )( clock() - m_timer ) / CLOCKS_PER_SEC;

	m_timer=clock();

	#endif
	//////////////////////////////////////////////////////////////////////////////////

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
