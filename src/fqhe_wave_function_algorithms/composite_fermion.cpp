////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 22/12/2014
//!
//!  \file 
//!		This cpp file contains functions to evaluate composite fermion and
//!		Bonderson-Slingerland wave functions. The algorithms themselves are
//!		written as templates in composite_fermion.h.
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

#include "composite_fermion.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Default constructor for the composite fermion data struct
//!
//! Sets default parameter values
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::CompositeFermionData::CompositeFermionData()
    :
    halfFluxAttatched(1),
    effectiveMonopole(1),
    LLup(2),
    LLdown(0),
    NEF(false),
    spinPolarised(true),
    altProjection(false),
    bsMode(false),
    unprojMode(false)
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Sets parameter values using command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::CompositeFermionData::InitFromCommandLine(
    #if _ENABLE_MPI_
    boost::program_options::variables_map* options,  //!<    Command line argument list
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
    #else
    boost::program_options::variables_map* options)  //!<    Command line argument list
    #endif
{
    #if _ENABLE_MPI_
    if(0 == mpi.m_id)    //  For the master node
    #endif
    {
        LLup          = (*options)["llup"].as<int>();
        LLdown        = (*options)["lldown"].as<int>();    
        altProjection = (*options)["altProjection"].as<bool>();
        unprojMode    = (*options)["un-proj"].as<bool>();        
    }
    
    #if _ENABLE_MPI_
    this->MpiSync(0,mpi);
    #endif
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#if _ENABLE_MPI_ 

////////////////////////////////////////////////////////////////////////////////
//!	\brief	MPI - communicate values required on parallel processes to the other
//! nodes
//!
////////////////////////////////////////////////////////////////////////////////
                                     
void FQHE::CompositeFermionData::MpiSync(
    const int syncId,   //!<    Node to sync with
    const utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    mpi.Sync(&LLup,1,syncId);
	mpi.Sync(&LLdown,1,syncId);
	mpi.Sync(&altProjection,1,syncId);
	mpi.Sync(&unprojMode,1,syncId);
}	
	
#endif

////////////////////////////////////////////////////////////////////////////////
//!	\brief The constructor for a CompositeFermion class.
//!
//! The constructor prepares the object for evaluation of a selected
//! composite fermion wave function/ Bonderson-Slingerland wave function.
//!     - Interprets argument data to decide what type of CF wave function 
//!     to prepare for
//!     - Pre calculates binomial coefficients and factorials required	
//!     - Allocates memory to store the Slater determinant and any required
//!     intermediate data for the calculation
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::CompositeFermion::CompositeFermion(
	FQHE::WaveFunctionData *wavefunctionData,		 	 //!<	Address of a struct of wave function data
	FQHE::CompositeFermionData *compositeFermionData)	 //!<	Address of a struct of composite fermion data
	: FQHE::WaveFunction(wavefunctionData),		         //  	Call base class constructor
	m_cfData(compositeFermionData)
{
	//////  1.  INTERPRET ARGUMENT DATA TO FIND WHAT TYPE OF 	////////////////////
	//////	CF WAVE FUNCTION TO EVALUATE

	//	Declare instance of pre-data struct to contain pre-calculated data for the algorithm

	m_preData = new PreData;

	//	NB in the case of CF wave functions the variable "Jastrow exponent" actually contains 
	//	half the "flux attached" (so the wave function contains (z_i-z_j)^(2*jastrowExponent)
	
	m_cfData->halfFluxAttatched=wavefunctionData->jastrowExponent;

	//	Temporary working values of LLup and LLdown for the constructor

	int LLup = m_cfData->LLup,LLdown=m_cfData->LLdown;

	// test if the effective Landau levels are fully occupied with electrons
	
	bool fullyOccupied=true;
	
	if(m_wfData->type==_BONDERSON_SLINGERLAND_)
	{
		m_cfData->bsMode=true;
	}
	else
	{
		m_cfData->bsMode=false;
	}

	//	Resolve conflicts
	
	if(LLup==0 &&LLdown==0)
	{
		std::cerr<<"\n\tERROR: CF LLs cannot both be set to zero."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}
	
	if((LLup<0&&LLdown>0) || (LLup>0 &&LLdown<0))
	{
		std::cerr<<"\n\tERROR: CF LLs should have the same sign."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}
	
	if(abs(LLdown)>abs(LLup))
	{
		//	Set up so that the order of specified LLs does not matter
		m_cfData->LLup=LLdown;
		m_cfData->LLdown=LLup;
		LLup=m_cfData->LLup;
		LLdown=m_cfData->LLdown;
	}
	
	if(LLdown==0)
	{
		m_cfData->spinPolarised=true;
	}
	else
	{
		m_cfData->spinPolarised=false;
	}
	
	if(m_wfData->jastrowExponent<=0)
	{
		std::cerr<<"\n\tWARNING: p<=0. Must have at least 1 flux attached. Setting fluxAttatched=1."<<std::endl<<std::endl;
		m_wfData->jastrowExponent=1;
	}
	
	if(m_cfData->unprojMode && !m_cfData->spinPolarised)
	{
		std::cerr<<"\n\tERROR: un-projected and spin un-polarised not available."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}
	
	//	Detect +ve or -ve effective field based on sign of LLup variable
	
	if(LLup<0)
	{
	   //	Negative effective field
	   
	   m_cfData->NEF=true;
	   m_wfData->fillNumerator=-LLup-LLdown;
	   m_wfData->fillDenominator=-((int)2*m_wfData->jastrowExponent*(LLup+LLdown)+1);
	}
	else
	{
	   //	Positive effective field
	   
	   m_cfData->NEF=false;
	   m_wfData->fillNumerator=LLup+LLdown;
	   m_wfData->fillDenominator=(int)2*m_wfData->jastrowExponent*(LLup+LLdown)+1;
	}  
	
	if(m_cfData->NEF && m_wfData->jastrowExponent!=1)
	{
		std::cerr<<"\n\tERROR: For negative effective field, no coded algorithm for flux attatched >=2."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}

	//	Check for partially filled effective Landau levels and set the monopole
	//	strength to be that of a fully filled set of lls
	
	//	Determine monopole strength m_preData->effectiveMonopole
	// See Appendix A of "Spinful CFs in a negative effective field"

	m_cfData->effectiveMonopole=(double)((int)m_wfData->nbr-LLup*LLup-LLdown*LLdown)/(2*LLup+2*LLdown);
	m_preData->effectiveMonopole = m_cfData->effectiveMonopole;
	m_wfData->monopoleStrength = (m_wfData->nbr*m_wfData->fillDenominator/m_wfData->fillNumerator-m_wfData->shift)/2.0;
    
	//	check if 2*effectiveMonopole is an integer

	double intPart;
			
	if(modf(2*m_preData->effectiveMonopole,&intPart)!=0)
	{
		//	if it's not an integer, then we have partially filled LLs
		//	Round effectiveMonopole up to the nearest filled LL value
		
		if(m_preData->effectiveMonopole>=0)
		{
			m_preData->effectiveMonopole=ceil(2*m_preData->effectiveMonopole)/2.0;
		}
		else
		{
			m_preData->effectiveMonopole=-ceil(2*fabs(m_preData->effectiveMonopole))/2.0;
		}

		fullyOccupied=false;
		
		std::cerr<<"\tWARNING PARTIALLY FILLED LL."<<std::endl<<std::endl;
	}

	
	if(!m_cfData->spinPolarised && !fullyOccupied)
	{
		std::cerr<<"\n\tERROR: Cannot work with partially occupied spinful CF wave functions."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}
	
	if(!m_cfData->NEF && m_cfData->LLdown>1 &&  m_cfData->LLup!=m_cfData->LLdown)
	{
		std::cerr<<"\n\tERROR: Cannot work with such cases yet - a bug remains in NonPolarisedPEF algorithm in this case."<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}

	m_wfData->nbrUp=(int)(2.0*m_preData->effectiveMonopole*LLup+LLup*LLup);
	m_wfData->nbrDown=(int)(2.0*m_preData->effectiveMonopole*LLdown+LLdown*LLdown);
	
	m_wfData->shift=m_wfData->nbr*(double)m_wfData->fillDenominator/m_wfData->fillNumerator-(2*m_preData->effectiveMonopole+2*m_wfData->jastrowExponent*(m_wfData->nbr-1));
	
	if(m_wfData->statistics==_BOSONS_)
	{
		m_wfData->fillDenominator-=m_wfData->fillNumerator;
		m_wfData->shift--;
	}
	
	if(m_cfData->bsMode)
	{
		if(!m_cfData->spinPolarised)
		{
			std::cerr<<"\n\tERROR: Bonderson--Slingerland wavefunctions ";
			std::cerr<<"are only available for the spin polarised case."<<std::endl<<std::endl;
			exit(EXIT_FAILURE);
		}
		m_wfData->shift++;
		m_wfData->fillDenominator+=m_wfData->fillNumerator;
			
	}
	
	if(m_cfData->spinPolarised)
	{
		m_cfData->altProjection=false;
	}
	   
	m_preData->effectiveMonopole=fabs(m_preData->effectiveMonopole);  // Take |Q| for use in subsequent calculations
	LLup=abs(LLup);
	LLdown=abs(LLdown);
	m_cfData->LLup=abs(LLup);
	m_cfData->LLdown=abs(LLdown);
	
	m_wfData->fillingFactor = (double)m_wfData->fillNumerator/m_wfData->fillDenominator;
	m_wfData->flux = (double)m_wfData->nbr/(m_wfData->fillingFactor)-m_wfData->shift;

	//	Determine default occupation numbers and Lz values
	//	These can be changed by calling the changeLLdata(llOccupationNumbers,cfLz) function
	
	m_preData->llOccupationNumbers = new (std::nothrow) double[LLup];
	m_preData->cfLz = new (std::nothrow) double[m_wfData->nbr];
    
    {
	    double m;
	    int d,ll,degen,totdegen;
	
	    ll=0;
	    m=-m_preData->effectiveMonopole;
	    degen=(int)(2.0*m_preData->effectiveMonopole)+1;
	    totdegen=degen;

	    m_preData->llOccupationNumbers[ll]=degen;
	
	    for(d=0;d<degen;d++,m++)
	    {
		    m_preData->cfLz[d]=m;
		    utilities::cout.DebuggingInfo()<<m_preData->cfLz[d]<<std::endl;
		
		    if(!m_cfData->spinPolarised)
		    {
			    //	set the values for the spin-down LL also
			    m_preData->cfLz[d+m_wfData->nbrUp]=m;
		    }
	    }

	    utilities::cout.DebuggingInfo()<<"occ "<<m_preData->llOccupationNumbers[ll]<<std::endl;
	
	    while(ll+1<LLup)
	    {
		    ll++;
		    m=-m_preData->effectiveMonopole-ll;
		    degen=(int)(2.0*(m_preData->effectiveMonopole+ll))+1;
		    m_preData->llOccupationNumbers[ll]=degen;

		    for(d=totdegen;d<(totdegen+degen);d++,m++)
		    {
			    if(d>=(int)m_wfData->nbrUp)	
			    {		
				    m_preData->llOccupationNumbers[ll]=d-totdegen;
				    break;
			    }

			    m_preData->cfLz[d]=m;
			    utilities::cout.DebuggingInfo()<<m_preData->cfLz[d]<<std::endl;
		    }
		    utilities::cout.DebuggingInfo()<<"occ "<<m_preData->llOccupationNumbers[ll]<<std::endl;
		    totdegen+=degen;
	    }
	
	    //	Do the same for spin-down LLs
	
	    ll=0;
	
	    degen=(int)(2.0*m_preData->effectiveMonopole)+1;
	    totdegen=degen;
	
	    while(ll+1<LLdown)
	    {
		    ll++;
		    m=-m_preData->effectiveMonopole-ll;
		    degen=(int)(2.0*(m_preData->effectiveMonopole+ll))+1;

		    for(d=totdegen;d<(totdegen+degen);d++,m++)
		    {
			    m_preData->cfLz[d+m_wfData->nbrUp]=m;
			    utilities::cout.DebuggingInfo()<<m_preData->cfLz[d]<<std::endl;
		    }
		    utilities::cout.DebuggingInfo()<<"occ "<<m_preData->llOccupationNumbers[ll]<<std::endl;
		    totdegen+=degen;
	    }
	
	    //	Check to make sure that all the lls contain some number of electrons,
	    //	Otherwise quit
	
	    for(ll=0;ll<LLup;ll++)
	    {
		    if(m_preData->llOccupationNumbers[ll]==0)
		    {
			    std::cerr<<"\n\tERROR: CF LL occupation too low - use fewer LLs."<<std::endl<<std::endl;
			    exit(EXIT_FAILURE);
		    }
	    }
	}
	
	//////      Generate output file name       ////////////////////////////////////
	
	{
	    std::stringstream name;

	    name.str("");
	
	    if(m_cfData->bsMode)
	    {
		    if(m_cfData->NEF)
		    {
			    name<<"Pfaffian*2CF_(-"<<LLup<<",-"<<LLdown<<") Bondersen--Slingerland (negative effective field)";
		    }
		    else
		    {
			    name<<"Pfaffian*2CF_("<<LLup<<","<<LLdown<<") Bondersen--Slingerland (positive effective field)";
		    }
	    }
	    else if(m_cfData->NEF)
	    {
		    name<<"2CF_(-"<<LLup<<",-"<<LLdown<<") Composite Fermion (negative effective field)";
	    }
	    else
	    {
		    name<<"2CF_("<<LLup<<","<<LLdown<<") Composite Fermion (positive effective field)";
	    }

	    //	Generate file name identifier
	    name.str("");
	    m_wfData->GenerateFileName();
	    name<<m_wfData->wfFileName;
	
	    if(m_cfData->NEF)
	    {
		    name<<"_LL_-"<<LLup;
		    if(LLdown==0)
		    {
			    name<<LLdown;
		    }
		    else
		    {
			    name<<"-"<<LLdown;
		    }
	    }
	    else
	    {
		    name<<"_LL_"<<LLup<<LLdown;
	    }
	
	    m_wfData->wfFileName = name.str();
	}

    //  Generate sphere/disc radius value
    m_wfData->InitRadius();

	//  Print a summary
	m_wfData->CheckAndPrint();	
		
	utilities::cout.MainOutput()<<"\tNo. 'flux attached'\t\t"<<2*m_wfData->jastrowExponent;
	utilities::cout.MainOutput()<<"\n\tEffective monopole strength\t"<<2*m_preData->effectiveMonopole;
	utilities::cout.MainOutput()<<"\n\tn_up is\t\t\t\t"<<	m_wfData->nbrUp;
	utilities::cout.MainOutput()<<"\n\tn_down is\t\t\t"<<	m_wfData->nbrDown;
	utilities::cout.MainOutput()<<"\n\tLLs "<<((fullyOccupied) ? ("fully") : ("partially"))<<" occupied";
	utilities::cout.MainOutput()<<"\n\tAlternative projection\t\t"<<((m_cfData->altProjection) ? ("on") : ("off"));
	utilities::cout.MainOutput()<<"\n\tNo projection at all\t\t"<<((m_cfData->unprojMode) ? ("on") : ("off"));
	utilities::cout.MainOutput()<<std::endl<<std::endl;

//////  2 and 3.   PRE-CALCULATE BINOMIALS AND ALLOCATE MEMORY REQUIRED 	////
//////	FOR INTERNAL OPERATIONS    

	//	set the maximum size for binomials stored

	m_preData->maxBinom=(int)(2*m_preData->effectiveMonopole+std::max(abs(m_cfData->LLup),abs(m_cfData->LLdown)))+1;

	//	To save processing time, store a list of binomial coefficients

	m_preData->binomList=new (std::nothrow) long int[m_preData->maxBinom*m_preData->maxBinom];

	for(int i=0;i<m_preData->maxBinom;i++)
	{
	   for(int j=0;j<m_preData->maxBinom;j++)
	   {
		   m_preData->binomList[i*m_preData->maxBinom+j]=utilities::Binomial<long int>(i,j);
		   utilities::cout.DebuggingInfo()<<"\t"<<i<<" C "<<j<<" = "<<m_preData->binomList[i*m_preData->maxBinom+j]<<std::endl;
	   }
	}
	
	//	allocate array to store number of terms in each LL (negative effective field case only)
	m_preData->nbrEachLL=0;
	if(m_cfData->NEF)
	{
		m_preData->nbrEachLL= new int[std::max(abs(m_cfData->LLup),abs(m_cfData->LLdown))];
		this->GeneratePreFactorNumber(m_wfData->nbr,std::max(abs(m_cfData->LLup),abs(m_cfData->LLdown)),m_preData->nbrEachLL);
	}
	
	////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_SPHERE_GEOMETRY_
	
	//	declare instances of template classes for each type of cf wave function evaluation
	//	in the sphere geometry case

	//	m_lowPrec uses standard double and complex<double> variables
	m_lowPrec = new ArbitraryPrecisionCfWavefunctions<dcmplx,double,64>(wavefunctionData,compositeFermionData,m_preData);

	//	high precision using the mpf type form the gmp library
	m_precLevel1 = new ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL1_>(wavefunctionData,compositeFermionData,m_preData);
	m_precLevel2 = new ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL2_>(wavefunctionData,compositeFermionData,m_preData);
	m_precLevel3 = new ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL3_>(wavefunctionData,compositeFermionData,m_preData);
	m_precLevel4 = new ArbitraryPrecisionCfWavefunctions<mpfCmplx,mpf_t,_PREC_LEVEL4_>(wavefunctionData,compositeFermionData,m_preData);

	//	NOTE the maximum precision that can be used at the moment is _PREC_LEVEL4_ (see header file for definitions)
	
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Destructor for CompositeFermion class.
//!
//!		Used to deallocate memory required for internal calculations 
//!		Called automatically when the object is deallocated	
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::CompositeFermion::~CompositeFermion()
{	
	//	Call destructors for low and high precision variables

	////////////////////////////////////////////////////////////////////////////////
	#if _ENABLE_SPHERE_GEOMETRY_
	
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		utilities::cout.SecondaryOutput()<<"m_lowPrec : ";
		#endif
		////////////////////////////////////////////////////////////////////////////////
		
		delete m_lowPrec;
		
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		utilities::cout.SecondaryOutput()<<"m_precLevel1 : ";
		#endif
		////////////////////////////////////////////////////////////////////////////////
		
		delete m_precLevel1;
			
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		utilities::cout.SecondaryOutput()<<"m_precLevel2 : ";
		#endif
		////////////////////////////////////////////////////////////////////////////////
		
		delete m_precLevel2;
			
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		utilities::cout.SecondaryOutput()<<"m_precLevel3 : ";
		#endif
		////////////////////////////////////////////////////////////////////////////////
		
		delete m_precLevel3;
		
		////////////////////////////////////////////////////////////////////////////////
		#if _CF_BENCHMARK_MODE_
		utilities::cout.SecondaryOutput()<<"m_precLevel4 : ";
		#endif
		////////////////////////////////////////////////////////////////////////////////
		
		delete m_precLevel4;
	
	#endif
	////////////////////////////////////////////////////////////////////////////////

	//	Delete some other memory allocations

	utilities::cout.DebuggingInfo()<<"delete m_cfData->m_preData->binomList"<<std::endl;
	delete[] m_preData->binomList;
	utilities::cout.DebuggingInfo()<<"delete m_cfData->m_preData->llOccupationNumbers"<<std::endl;
	if(m_preData->llOccupationNumbers!=0)	delete[] m_preData->llOccupationNumbers;
	utilities::cout.DebuggingInfo()<<"delete m_preData->cfLz"<<std::endl;
	if(m_preData->cfLz!=0)					delete[] m_preData->cfLz;
	utilities::cout.DebuggingInfo()<<"did dtor"<<std::endl;
	if(m_preData->nbrEachLL!=0)				delete[] m_preData->nbrEachLL;
	
	delete m_preData;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function evaluates a composite fermion wave function of the
//!	type specified by the class parameters, and for the sphere geometry.
//!
//! FIRST: make an object in the CompositeFermion class	
//! The data required to specify the wave function is contained in the
//! two data structs (see the description of the class constructor).
//!
//! Example of data structure set up:
//!     wfData.geometry=sphere 		Sets the sphere geometry
//!     wfData.type=CompositeFermion	Sets a composite fermion wave function
//!     wfData.nbr=36				Sets the number of particles for evaluation
//!     wfData.jastrowExponent=1	This will be used as half the number of 
//!										flux attached
//!     cfData.statistics=_FERMIONS_	evaluate for fermions
//!
//!     cfData.LLup=2			Set two spin-up effective LLs
//!									(- sign here sets negative effective field)
//!     cfData.LLdown=0			Set no spin down effective LLs (polarized)
//!
//!     cfData.altProjection=false	Turn off separate projection of spin-up
//!										and spin-down wave functions (do combined 
//!										projection).
//!     cfData.unprojMode=false	Set to true to use the unprojected wave
//!									function i.e.  (IQHE*jastrow factor)
//!
//! The the EvaluateWfSphere method can be called on the object to return
//! the wave function value	for a given set of co-ordinates	
//!	
//! \return Returns log of wave function value.
//!
////////////////////////////////////////////////////////////////////////////////
	
dcmplx FQHE::CompositeFermion::EvaluateWfSphere(
	const int n,    //!<	The number of particles to perform the evaluation for
	dcmplx *u,	    //!<	The memory address of the start of an array of spinor u coordinates
	dcmplx *v)	    //!<	The memory address of the start of an array of spinor v coordinates
	const
{
    int maxPrecLevel = 0;       //  Track the maximum precision level used
	dcmplx testVal,testVal1;	//	Test evaluations of the wave function as different
								//	levels of precision	
	dcmplx wfValue = 0.0;	    //	Value of wave function

	if(m_wfData->nbr!=n)
	{
		std::cerr<<"\tevaluate_wf_sphere WARNING: argument n does not match object WaveFunctionData nbr paramater."<<std::endl;
		getchar();
	}

	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_

		utilities::cout.DebuggingInfo().precision(15);

		utilities::cout.DebuggingInfo()<<"Here are the u co-ordinates:"<<std::endl;
		for (int j=0;j<n;j++){utilities::cout.DebuggingInfo()<<real(*(u+j))<<"+"<<imag(*(u+j))<<"I,"<<std::endl;}
		utilities::cout.DebuggingInfo()<<"Here are the v co-ordinates:"<<std::endl;
		for (int j=0;j<n;j++){utilities::cout.DebuggingInfo()<<real(*(v+j))<<"+"<<imag(*(v+j))<<"I,"<<std::endl;}

		utilities::cout.DebuggingInfo()<<"LL occupation numbers:"<<std::endl;

		for(int i=0;i<m_cfData->LLup;i++)
		{
			utilities::cout.DebuggingInfo()<<m_preData->llOccupationNumbers[i]<<std::endl;
		}

		utilities::cout.DebuggingInfo()<<"Lz values for each electron:"<<std::endl;

		for(int i=0;i<n;i++)
		{
			utilities::cout.DebuggingInfo()<<m_preData->cfLz[i]<<std::endl;
		}

		getchar();
	#endif
	////////////////////////////////////////////////////////////////////////////////

    //  Test for any particle being close to the north pole - if so , perform a 
    //  global roation to move it away from that position. This will not affect
    //  the value of the wave function, but will improve the stability of the
    //  algorithm
    
    //dcmplx* uTemp = 0;
    //dcmplx* vTemp = 0;
    
    //while(this->IsCloseToPole(n,u,v))
    //{
    //    //  Keep track of the existing co-ordinates
    //    
    //    uTemp = new (std::nothrow) dcmplx[n];
    //    vTemp = new (std::nothrow) dcmplx[n];
    //    
    //    memcpy(uTemp,u,n*sizeof(dcmplx));
    //    memcpy(vTemp,v,n*sizeof(dcmplx));
    //
    //    this->RotateAwayFromPole(n,u,v);
    //}

	//	The following if statements will call the appropriate wave function algorithm

    double conditionNumber;

	if(m_cfData->NEF)       //	NEGATIVE EFFECTIVE FIELD CF WAVE FUNCTIONS
	{
		//	Calculate the condition number
		//	This will be used to determine whether to go to high precision or not

		conditionNumber = m_lowPrec->CalculateConditionNumber(n,u,v);

		if(m_cfData->spinPolarised)     //	SPIN POLARIZED
		{
			testVal=m_lowPrec->PolarisedNEF(n,m_cfData->LLup,u,v);

            if(std::isinf(real(testVal)) || std::isnan(real(testVal)))
            {
                std::cerr<<"RECORDED OVERFLOW ERROR"<<std::endl;
            }

			if(conditionNumber>FQHE::condTolL2 || std::isinf(real(testVal)) || std::isnan(real(testVal)))
			{
			    maxPrecLevel++;
			    testVal1 = m_precLevel1->PolarisedNEF(n,m_cfData->LLup,u,v);
			
			    ////////////////////////////////////////////////////////////////////////////////
			    #if _DEBUG_
			    utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
			    #endif
			    ////////////////////////////////////////////////////////////////////////////////

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->PolarisedNEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->PolarisedNEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////

				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->PolarisedNEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE.SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}
			}
			else
			{
			    testVal1 = testVal;
			}

			wfValue = testVal1;
		}
		else if(m_cfData->altProjection)		//	SPIN NON-POLARIZED SEPARATE UP AND DOWN PROJECTION
		{
			testVal = m_lowPrec->PolarisedNEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
			testVal += m_lowPrec->PolarisedNEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

            if(std::isinf(real(testVal)) || std::isnan(real(testVal)))
            {
                std::cerr<<"RECORDED OVERFLOW ERROR"<<std::endl;
            }
    
			if(conditionNumber>FQHE::condTolL2 || std::isinf(real(testVal)) || std::isnan(real(testVal)))
			{
			    maxPrecLevel++;
			    testVal1 = m_precLevel1->PolarisedNEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
			    testVal1 += m_precLevel1->PolarisedNEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);
			
			    ////////////////////////////////////////////////////////////////////////////////
			    #if _DEBUG_
			    utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
			    #endif
			    ////////////////////////////////////////////////////////////////////////////////
			
				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->PolarisedNEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel2->PolarisedNEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->PolarisedNEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel3->PolarisedNEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->PolarisedNEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel4->PolarisedNEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE.SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}
			}
			else
			{
			    testVal1 = testVal;
			}

			wfValue = testVal1;
		}
		else    //	SPIN NON-POLARIZED COMBINED UP-DOWN PROJECTION
		{

			testVal = m_lowPrec->NonPolarisedNEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

            if(std::isinf(real(testVal)) || std::isnan(real(testVal)))
            {
                std::cerr<<"RECORDED OVERFLOW ERROR"<<std::endl;
            }

			if(conditionNumber>FQHE::condTolL2 || std::isinf(real(testVal)) || std::isnan(real(testVal)))
			{
			    maxPrecLevel++;
                testVal1 = m_precLevel1->NonPolarisedNEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);
                
                ////////////////////////////////////////////////////////////////////////////////
			    #if _DEBUG_
			    utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
			    #endif
			    ////////////////////////////////////////////////////////////////////////////////
                
				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->NonPolarisedNEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->NonPolarisedNEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->NonPolarisedNEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE.SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}
			}
			else
			{
			    testVal1 = testVal;
			}

			wfValue = testVal1;
		}
	}
	else								//	POSITIVE EFFECTIVE FIELD CF WAVE FUNCTIONS
	{
		if(m_cfData->spinPolarised)     //	SPIN POLARIZED
		{
			if(n>FQHE::maxSwitch)			//	above a certain number, use high precision
			{
				testVal = m_lowPrec->PolarisedPEF(n,m_cfData->LLup,u,v);
				testVal1 = m_precLevel1->PolarisedPEF(n,m_cfData->LLup,u,v);
				maxPrecLevel++;

				////////////////////////////////////////////////////////////////////////////////
				#if _DEBUG_
				utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
				#endif
				////////////////////////////////////////////////////////////////////////////////

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->PolarisedPEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->PolarisedPEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////

				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->PolarisedPEF(n,m_cfData->LLup,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////

				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE. SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}

				wfValue = testVal1;
			}
			else
			{
				wfValue = m_lowPrec->PolarisedPEF(n,m_cfData->LLup,u,v);
			}

		}
		else if(m_cfData->altProjection)		//	SPIN NON-POLARIZED SEPARATE UP AND DOWN PROJECTION
		{
			if(n>FQHE::maxSwitch)			//	above a certain number, use high precision
			{
				testVal = m_lowPrec->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
				testVal += m_lowPrec->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

				testVal1 = m_precLevel1->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
				testVal1 += m_precLevel1->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

                maxPrecLevel++;

				////////////////////////////////////////////////////////////////////////////////
				#if _DEBUG_
				utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
				#endif
				////////////////////////////////////////////////////////////////////////////////

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel2->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel3->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
					testVal1 += m_precLevel4->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE. SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}

				wfValue = testVal1;
			}
			else
			{
				wfValue = m_lowPrec->PolarisedPEF(m_wfData->nbrUp,m_cfData->LLup,u,v);
				wfValue += m_lowPrec->PolarisedPEF(m_wfData->nbrDown,m_cfData->LLdown,u+m_wfData->nbrUp,v+m_wfData->nbrUp);
			}

		}
		else    //	SPIN NON-POLARIZED COMBINED UP-DOWN PROJECTION
		{
			if(n>FQHE::maxSwitch)			//	Above a certain number, use high precision
			{
				testVal = m_lowPrec->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);
				testVal1 = m_precLevel1->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);
    
                maxPrecLevel++;

				////////////////////////////////////////////////////////////////////////////////
				#if _DEBUG_
				utilities::cout.DebuggingInfo()<<"\n\tdcmplx value value: "<<testVal<<"\tmpf_t Level1-precision value: "<<testVal1<<std::endl;
				#endif
				////////////////////////////////////////////////////////////////////////////////

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel2->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level2-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel3->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level3-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					testVal = testVal1;

					testVal1 = m_precLevel4->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);

					////////////////////////////////////////////////////////////////////////////////
					#if _DEBUG_
					utilities::cout.DebuggingInfo()<<"\n\tmpf_t Level4-precision value: "<<testVal1<<std::endl;
					#endif
					////////////////////////////////////////////////////////////////////////////////
				}

				if(this->IncreasePrecCondition(testVal,testVal1,maxPrecLevel))
				{
					std::cerr<<"\n\tERROR: MAXIMUM PRECISION REACHED WITHOUT CONVERGENCE. SINGULAR WAVE FUNCTION OR OVERFLOW ERROR. "<<std::endl;
					std::cerr<<"\n\tPRESS ANY KEY TO CONTINUE..."<<std::endl;
					getchar();
				}

				wfValue = testVal1;
			}
			else
			{
				wfValue = m_lowPrec->NonPolarisedPEF(m_wfData->nbrUp,m_wfData->nbrDown,m_cfData->LLup,m_cfData->LLdown,u,v);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	//	Put in up-down Jastrow factor that was not projected in alternative projection method
	if(m_cfData->altProjection)
	{
		//multiply by Jastrow of (z_up - z_down)^2p
		for (int i=0;i<m_wfData->nbrUp;i++)
		{
			for (int j=m_wfData->nbrUp;j<n;j++)
			{
				wfValue+=4.0*m_wfData->jastrowExponent*log(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	//	For the bosons case, divide by a Jastrow factor
	if(m_wfData->statistics==_BOSONS_)
	{
		if(m_cfData->spinPolarised)
		{
			//divide by Jastrow of (z_i - z_j)

			for (int i=0;i<n;i++)
			{
				for (int j=i+1;j<n;j++)
				{
					wfValue-=log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
		}
		else
		{
			//divide by Jastrow of (z_up - z_down)
			for (int i=0;i<m_wfData->nbrUp;i++)
			{
				for (int j=m_wfData->nbrUp;j<n;j++)
				{
					wfValue-=log(*(u+i)**(v+j)-*(u+j)**(v+i));
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	//	For Bonderson-Slingerland case, multiply by a Pfaffian
	if(m_cfData->bsMode)
	{
		dcmplx *m_pfaff = new dcmplx[n*n];
		dcmplx *p_row;
		//	Calculate Pfaffian part

		//	Set upper triangular parts to be 1/(u[i]v[j]-u[j]v[i])

        p_row=m_pfaff;    
        
		for(int i=0;i<n;i++,p_row+=n)
		{
			for(int j=i+1;j<n;j++)
			{
				*(p_row+j)=dcmplx(1,0)/(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}

		//PrintPfaffian(n);

		wfValue+=utilities::linearAlgebra::LogPfaffian(m_pfaff,n);

		//PrintPfaffian(n);

		//	Multiply by a Jastrow factor

		for (int i=0;i<n;i++)
		{
			for (int j=i+1;j<n;j++)
			{
				wfValue+=log(*(u+i)**(v+j)-*(u+j)**(v+i));
			}
		}

		delete[] m_pfaff;
	}

	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_
    utilities::cout.DebuggingInfo()<<"\tCondition Number: "<<conditionNumber<<std::endl;
    utilities::cout.DebuggingInfo()<<"\tFinal wave function value : "<<wfValue<<std::endl;
    utilities::cout.DebuggingInfo()<<"\tPREC LEVEL:\t"<<maxPrecLevel<<std::endl;
    getchar();
	#endif
	////////////////////////////////////////////////////////////////////////////////

	if(std::isnan(real(wfValue)) || std::isinf(real(wfValue)))
	{
		std::cerr<<"\tWARNING: nan/inf wave function value detected "<<std::endl;
		getchar();
	}

    //  If we rotated the co-ordinates, return the values in the u and v arrays 
    //  to their original values

    //if(uTemp != 0)
    //{
    //    memcpy(u,uTemp,n*sizeof(dcmplx));
    //    memcpy(v,vTemp,n*sizeof(dcmplx));
    //    
    //    delete[] uTemp;
    //    delete[] vTemp;
    //}

	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////
#if _ENABLE_DISC_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function is intended for evaluation of CF wave functions in the 
//!	disc geometry. THIS FUNCITON IS NOT COMPLETED.
//!
//!	\return	The log of the wave function. N.B. the result is not normalized
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::CompositeFermion::EvaluateWfDisc(
	const int n,		//!<	The number of particles to perform the evaluation for
	dcmplx* z)          //!<	The memory address of the start of an array of disc coordinates
	const
{
	//    !!!!!	  THIS FUNCITON IS NOT COMPLETED    !!!!!

	// 	Composite fermion wave function algorithm for disc geometry
	dcmplx wfValue=0.0;

	for(int i=0;i<n;i++)
	{
		//	Implement Gaussian exponential part
		double tmp=abs(*(z+i));
		wfValue-=tmp*tmp/4;
	}

	//    !!!!!	  THIS FUNCITON IS NOT COMPLETED    !!!!!

	return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function allows the user to reset the internal cfLz and
//!	llOccupationNumbers variables. This would apply, for example, where
//!	one wants to calculate a partially filled CF Slater determinant
//!	
//!	Sets internal PreData data values in the CompositeFermion class.
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::CompositeFermion::SetLLdata(
	const int n,	        //!<	Total number of particles in the wave function
							//!		(must equal the length of the cfLz array)
	const int nbrLL,		//!<	Total number of effective LLs in the wave function
							//!		(must equal the length of the llOccupationNumbers array)
	const int* cfLz,		//!<	Pointer to a list of all Lz values of the electrons
							//!		including all effective LLs e.g. (for LL=2 and n=10)
							//!		[-1.5,-0.5,0.5,1.5,-2.5,-1.5,-0.5,0.5,1.5,2.5]
	const double* llOccupationNumbers)
	                        //!<	Pointer to a list of occupation numbers	
							//!	    for each effective LL. e.g. for	LL=2 and n=10 is [4,6]	
{
	//	Check that n is the same as the number of particles:

	if(m_wfData->nbr!=n)
	{
		std::cerr<<"\tERROR: SetLLdata() first argument must be same as particle number."<<std::endl;
		exit(EXIT_FAILURE);
	}

	if(m_cfData->LLup!=nbrLL)
	{
		std::cerr<<"\tERROR: SetLLdata() second argument must be same as Landau level number."<<std::endl;
		exit(EXIT_FAILURE);
	}

	for(int i=0;i<n;i++)
	{
		m_preData->cfLz[i]=cfLz[i];
	}

	for(int j=0;j<nbrLL;j++)
	{
		m_preData->llOccupationNumbers[j]=llOccupationNumbers[j];
	}

	return;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Pre calculation of factorial and binomial factors for each term	
//!	occurring the the CF wave function. Determines the size of m_facList
//!	to be allocated.
//!
//!	This calculation only needs to be done in the negative effective
//!	case, which requires additional factorial factors.
//!
//!	\return Size of m_facList array to be allocated	
//!
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

void FQHE::CompositeFermion::GeneratePreFactorNumber(
	const int n,	    //!<	Number of particles
	const int nbrLL,    //!<	Number of Landau levels
	int* nbrEachLL)	    //!<	An array containing the cumulative number 
					    //!		of states in each LL
{
	//	Declare local variables	

	nbrEachLL[0] = 0;

	for(int d=0;d<m_preData->llOccupationNumbers[0];d++)
	{
		double m = m_preData->cfLz[d];

		for(int t=(int)(m_preData->effectiveMonopole+m);t<=(int)(n-1-m_preData->effectiveMonopole+m);t++)
		{
			nbrEachLL[0]++;
		}
	}

	utilities::cout.DebuggingInfo()<<"\tNumber in ll "<<0<<"\t"<<nbrEachLL[0]<<std::endl;

    for(int ll=1,degen=0,totdegen=0;ll<nbrLL;ll++,totdegen+=degen)
    {
        nbrEachLL[ll] = nbrEachLL[ll-1];
		degen = m_preData->llOccupationNumbers[ll];
		
		for(int d=(int)totdegen;d<(int)(totdegen+degen);d++)
		{
			double m = m_preData->cfLz[d];

			for(int s=0;s<=ll;s++)
			{
				//	skip if we are multiplying by a 0 binomial coefficient
				if(m_preData->effectiveMonopole+m+s>=0 && m_preData->effectiveMonopole+m+s<=2*m_preData->effectiveMonopole+ll)
				{
					for(int t=(int)(m_preData->effectiveMonopole+m+s);t<=(int)(n-1-m_preData->effectiveMonopole+m-ll+s);t++)
					{
						if(t>=0 && t<n)
						{
							nbrEachLL[ll]++;
						}
					}
				}
			}
		}
		utilities::cout.DebuggingInfo()<<"\tNumber in ll "<<ll<<"\t"<<nbrEachLL[ll]<<std::endl;
	}

	return;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function tests CF wave functions evaluated at different precision
//!	levels, and tries to determine if it is necessary to go to even higher
//!	precision level to obtain an accurate result.
//!
//!	\return	Returns true if even higher precision is required
//!
////////////////////////////////////////////////////////////////////////////////

bool FQHE::CompositeFermion::IncreasePrecCondition(
	const dcmplx lowPrecResult,		//!<	A low precision result for the wave function value
	const dcmplx highPrecResult,	//!<	A higher precision result for the wave function value
	int& maxPrecLevel)              //!<    Set the precision level
	const
{
	//	Set the condition for increasing the precision of a calculation
	bool decision;

	decision = fabs(real(highPrecResult-lowPrecResult))>FQHE::tol;
	decision = decision || std::isinf(real(lowPrecResult)) || std::isnan(real(lowPrecResult));

	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_
	utilities::cout.DebuggingInfo()<<fabs(real(highPrecResult-lowPrecResult))<<std::endl;
	if(decision)    maxPrecLevel++;
	#endif
	////////////////////////////////////////////////////////////////////////////////

	return decision;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//! \brief Calculate the normalization factor of the monopole harmonic functions
//!
//!	\return Value of the normalization coefficient.
//!
////////////////////////////////////////////////////////////////////////////////

double FQHE::MonopoleHarmonicNorm(
    const double q,   //!<    monopole strength
    const int l,      //!<    Landau level index
    const double m)   //!<    z-component of angular momentum: from -q-l to q+l   
{
    double norm = std::pow(2.0,(-q-l))*sqrt((2*(q+l)+1)/(4*PI));

    if((q+l-m)>(q+l+m))
    {
        norm *= sqrt(utilities::DivFactorial<double>((q+l-m),(2.0*q+l))*utilities::DivFactorial<double>((q+l+m),l));
    }
    else
    {
        norm *= sqrt(utilities::DivFactorial<double>((q+l-m),l)*utilities::DivFactorial<double>((q+l+m),(2.0*q+l)));
    }

    return norm;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Test to see if any co-ordinate lies close to the north pole.
//!
////////////////////////////////////////////////////////////////////////////////

bool FQHE::CompositeFermion::IsCloseToPole(const int n,dcmplx* u,dcmplx* v) const
{
    //  Test if any of the u co-ordiantes has a very small value
    
    double tol = 0.0001;
    
    for(int i=0;i<n;i++)
    {
        if(abs(u[i])<tol)   return true;    
    }
    
    return false;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Attempt to rotate all co-ordinates to a frame where the wave 
//! function is invariant, but without particles close to the north pole
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::CompositeFermion::RotateAwayFromPole(const int n,dcmplx* u,dcmplx* v) const
{
	for(int k=0;k<n;k++)
	{
        //  Try reflecting all co-ordinates about the equator

		//double cosTheta = std::real(sqrt(1.0-4.0*u[k]*u[k]*v[k]*v[k]));

		//u[k] *= sqrt((1-cosTheta)/(1+cosTheta));
		//v[k] *= sqrt((1+cosTheta)/(1-cosTheta));
		
		u[k] = v[k];
		v[k] = u[k];
	}
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
