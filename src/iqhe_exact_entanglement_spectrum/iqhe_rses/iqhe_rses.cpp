////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 10/07/2014
//!
//!  \file 
//!		This file defines a class to calculate the numerically exact real
//!     space entanglement spectrum for integer quantum Hall states
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

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "iqhe_rses.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Specify the integrand function for overlap calculations, using
//! a syntax compatible with the gsl_integration library
//!
////////////////////////////////////////////////////////////////////////////////
double Integrand(
    double x,       //!<    Integration variable
    void* params)   //!<    A list of function parameters:
                    //!<    params[0] is the monopole strength
                    //!<    params[1] is the common z-component angular momentum
                    //!<    params[2] is the Landau level index of the first function
                    //!<    params[3] is the Landau level index of the second function                  
{
    const double q = ((double*) params)[0];
    const double m = ((double*) params)[1];
    const int ll1  = (int) ((double*) params)[2];
    const int ll2  = (int) ((double*) params)[3];
    
    IqheRses obj;

    double val = obj.EvaluateMonopoleHarmonic(q,ll1,m,x)*obj.EvaluateMonopoleHarmonic(q,ll2,m,x);
    
    //utilities::cout.DebuggingInfo()<<"current integrand value:"<<val<<std::endl;

    return val;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Class constructor
//!
////////////////////////////////////////////////////////////////////////////////

IqheRses::IqheRses()
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Class destructor
//!
////////////////////////////////////////////////////////////////////////////////

IqheRses::~IqheRses()
{}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Calculate the monopole harmonic function, which describes the single
//! particle orbitals in the quantum Hall problem in a spherical geometry
//!
////////////////////////////////////////////////////////////////////////////////

double IqheRses::EvaluateMonopoleHarmonic(
    const double q,   //!<    monopole strength
    const int l,      //!<    Landau level index
    const double m,   //!<    z-component of angular momentum: from -q-l to q+l 
    const double x)   //!<    position to evaluate the monopole harmonic at 
    const
{

    double sum = 0;
    const int max = (int)(q+l-m);
    
    for(int s=0;s<=max;s++)
    {
        //  Calculate -1^(q+l-m-s)
    
        const int sign = (int)(q+l-m-s) & 1 ? -1 : 1;
    
        sum +=  sign*utilities::BinomialFromTable(l,s)*utilities::BinomialFromTable((int)(2*q+l),(int)(q+l-m-s))*std::pow((1-x),(l-s))*std::pow((1+x),s);
    }

    sum *= std::pow((1-x),((q-m)/2.0))*std::pow((1+x),((q+m)/2.0));
    
    return sum;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Calculate the normalization factor of the monopole harmonic functions
//!
////////////////////////////////////////////////////////////////////////////////

double IqheRses::MonopoleHarmonicNorm(
    const double q,   //!<    monopole strength
    const int l,      //!<    Landau level index
    const double m)   //!<    z-component of angular momentum: from -q-l to q+l   
    const 
{
    double norm = std::pow(2.0,(-q-l))*sqrt((2*(q+l)+1)/(4*PI));

    if((q+l-m)>(q+l+m))
    {
        norm*=sqrt(utilities::DivFactorial<double>((q+l-m),(2.0*q+l))*utilities::DivFactorial<double>((q+l+m),l));
    }
    else
    {
        norm*=sqrt(utilities::DivFactorial<double>((q+l-m),l)*utilities::DivFactorial<double>((q+l+m),(2.0*q+l)));
    }

    return norm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Public interface to sequentially call all internal calculations
//!
////////////////////////////////////////////////////////////////////////////////

void IqheRses::GenerateIqheSpectrum(
    boost::program_options::variables_map* optionList,  //!<    Command line option list
    const int lzaMin,                                   //!<    Minimum sector calculated
    const int lzaMax)                                   //!<    Maximum sector calculated
{
    utilities::cout.MainOutput()<<"\n\t=========================    CALCUALTING IQHE RSES    ========================="<<std::endl;

    //  Import command line variables
    
    int nbrParticles = (*optionList)["nbr-n"].as<int>();
    int nbrA = (*optionList)["nbr-a"].as<int>();
    
    if(nbrParticles>_MAX_NUMBER_)
    {
        std::cerr<<"ERROR: IQHE spectrum not calculable for more than "<<_MAX_NUMBER_<<" particles - the code needs to be modified."<<std::endl;
        return;
    }
    
    int nbrLevels = (*optionList)["nbr-ll"].as<int>();
    double monopoleStrength = (double)(nbrParticles-nbrLevels*nbrLevels)/(2*nbrLevels);
    //  Dimension of the list of overlaps/pseudoenergies from -q-l to q+l with l=nbrLevels-1
    int nbrStates = (int)(1+2*(monopoleStrength+nbrLevels-1));

    utilities::cout.MainOutput()<<"\n\t\tMONOPOLE STRENGTH (=FLUX/2)\t\t"<<monopoleStrength<<std::endl;

    //////      Calculate list of overlap matrices
    
    //  Allocate memory to store overlap matrix list
     
    double* overlapMatrixBuffer = new double[nbrStates*nbrLevels*nbrLevels];

    //  Populate the list
    
    this->GenerateOverlapMatrix(overlapMatrixBuffer,monopoleStrength,nbrLevels,nbrStates);

    //  Allocate space to store new set of energy levels
    
    double* energyLevels = new double[nbrStates*nbrLevels];
    double normalization = 0.0;
    
    //  Generate new energy levels and normalization
    
    this->GenerateNewLevels(overlapMatrixBuffer,energyLevels,normalization,nbrLevels,nbrStates);
    
    //  Free memory associated with overlap matrices - these are no longer required
    
    delete[] overlapMatrixBuffer;
    
    for(int i=0;i<nbrStates*nbrLevels;i++)
    {
        //utilities::cout.DebuggingInfo()<<"energyLevel["<<i<<"] = "<<energyLevels[i]<<std::endl;
    }
    
    //  Represent occupation patterns of single particle orbital populations by 
    //  a fixed vector of 128-bit bitsets. Each bitset represents the occupation
    //  of single particle orbitals within one Landau level. The vector of these
    //  containers is a list of all possible occupations
    
    for(int s=lzaMin;s<=lzaMax;s++)
    {
        std::vector<std::vector<std::bitset<_MAX_NUMBER_>> > occupationsList;
        std::vector<int> spectrumLabels;
        
        occupationBasis::GenerateLevelOccupations(occupationsList,spectrumLabels,
                                       nbrParticles,nbrA,nbrLevels,s);
        
        utilities::cout.DebuggingInfo()<<"\n\t\tOCCUPATION PATTERNS FOR SECTOR "<<s<<std::endl;
        
        for(unsigned int i=0;i<occupationsList.size();i++)
        {
            utilities::cout.DebuggingInfo()<<utilities::cout.HyphenLine()<<std::endl;
            
            for(int j=0;j<nbrLevels;j++)
            {
                utilities::cout.DebuggingInfo()<<"\t\t"<<occupationsList[i][j]<<std::endl;
            }
        }
        
        #if _DEBUG_
        this->CheckAngularMomentum(occupationsList,monopoleStrength,nbrA,nbrLevels,s);
        #endif

        //  Allocate memory to store the RSES

        double* spectrum = new double[occupationsList.size()];
        
        occupationBasis::AssignEnergies(occupationsList,energyLevels,spectrum,normalization,nbrStates,nbrLevels);
        
        this->SpectrumToFile(optionList,spectrum,//spectrumLabels,
                             occupationsList.size(),s);
        
        delete[] spectrum;
    }

    delete[] energyLevels;

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the overlap matrix of monopole harmonic functions in the upper
//! half-sphere (assuming an equatorial real space cut) for functions with the 
//! same z-component of angular momentum but occurring in different Landau levels. 
//! 
//! This calculation is necessary because such functions are not orthonormal when
//! integrated over only half of the sphere. The overlap matrix generated here
//! is diagonalized in order to generate a new set of single particle entanglement
//! energies for construction of the RSES
//!
////////////////////////////////////////////////////////////////////////////////

void IqheRses::GenerateOverlapMatrix(
    double* overlapMatrixBuffer,    //!<    Buffer to store list of overlap matrices
    const double q,                 //!<    Monopole strength
    const int nbrLevels,            //!<    Number of Landau levels
    const int nbrStates)            //!<    Number of single particle orbitals
    const
{
    //  Numerical integration parameters
    
    gsl_function integrand;         //  gsl_function object sets integrand function
    double integrandParams[4];      //  List of integration parameters stored in double format   
    double lowerLimit = 0.0;        //  Lower integration limit (orthogonal for x between -1 and 1)
    double upperLimit = 1.0;        //  Upper integration limit
    double epsabs = 1e-10;          //  Set relative error tolerance
    double epsrel = 1e-20;          //  Set absolute error tolerance
    double result;                  //  Container for integration result
    double error;                   //  Container for integration error
    size_t nevals;                  //  Container for number of function evaluations
    size_t maxIntervals=1000;       //  Maximum number of intervals in the integration routine

    //  Allocate work space for numerical integration routine

    gsl_integration_cquad_workspace* gsl_work=gsl_integration_cquad_workspace_alloc(maxIntervals);
    
    //  Populate the list of overlap matrices

    double* p_buffer = overlapMatrixBuffer;
    double zComponentAngularMomentum = -q-(nbrLevels-1);
    
    integrandParams[0] = q;
    integrand.function = &Integrand;
    
    for(int m=0;m<nbrStates;m++,zComponentAngularMomentum+=1.0)
    {
        integrandParams[1] = zComponentAngularMomentum;
    
        for(int i=0;i<nbrLevels;i++)
        {
            integrandParams[2] = i;
            const double iNorm =  this->MonopoleHarmonicNorm(q,i,zComponentAngularMomentum);
        
            for(int j=0;j<nbrLevels;j++)
            {
                integrandParams[3] = j;
                const double jNorm = this->MonopoleHarmonicNorm(q,j,zComponentAngularMomentum);
            
                //  Update integration parameters
                integrand.params = integrandParams;
            
                #if _DEBUG_
                std::cout<<"Integration parameters:\tq="<<q<<"\tl1="<<i<<"\tl2="<<j<<"\tm="<<zComponentAngularMomentum<<std::endl;
                #endif
            
                //  Calculate integral
                gsl_integration_cquad(&integrand,lowerLimit,upperLimit,epsabs,epsrel,gsl_work,&result,&error,&nevals);

                #if _DEBUG_
                std::cout<<"Integration result = "<<result<<" +/- "<<error<<std::endl;
                std::cout<<"norm = "<<2*PI*iNorm*jNorm<<std::endl;
                std::cout<<"total = "<<2*PI*iNorm*jNorm*result<<std::endl;
                std::cout<<"number of function evals = "<<nevals<<std::endl;
                #endif

                *(p_buffer) = 2*PI*iNorm*jNorm*result;
                
                p_buffer++;
            }
        }
    }
    #if _DEBUG_
    PAUSE();
    #endif
    
    //  Free integration workspace
    gsl_integration_cquad_workspace_free(gsl_work);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief This function calls a LAPACK routine to diagonalize the overlap matrix
//! generated by the GenerateOverlapMatrix() function. This defines a new set of 
//! entanglement energy levels
//!
////////////////////////////////////////////////////////////////////////////////

void IqheRses::GenerateNewLevels(
    double* overlapMatrixBuffer,    //!<    Buffer to store list of overlap matrices
    double* energyLevels,           //!<    On input, contain alpha_m
                                    //!     On output contain epsilon_m
    double& normalization,          //!<    A value for the normalization, to be set
    const int nbrLevels,            //!<    Number of Landau levels
    const int nbrStates)            //!<    Number of single particle orbitals   
    const
{
    double* p_buffer = overlapMatrixBuffer;
    double* p_levels = energyLevels;

    for(int m=0;m<nbrStates;m++)
    {
        utilities::linearAlgebra::DiagonalizeSymmetricMatrix<double>(p_buffer,p_levels,nbrLevels,'U');
    
        //  shift pointers to next matrix
        p_buffer+=nbrLevels*nbrLevels;
        p_levels+=nbrLevels;
    }
    
    //  Here we can exploit the fact that the spectrum is symmetrical 
    //  Otherwise numerical errors can build up from the calculation of the 
    //  overlap matrices
   
    //  Generate entanglement energy levels defined by the function log(e_n/(1-e_n))
    //  for a half-cut in real space (e_n are the eigenvalues of the overlap matrix)
    
    //  The normalization of the energy levels is defined by sum log(1-e_n)
    //  for all states.

    normalization = 0.0;
    
    utilities::cout.DebuggingInfo()<<std::endl;
    
    //  Temporary store for the new energy levels
    
    double* newLevels = new double[nbrStates*nbrLevels];
    
    p_levels = newLevels;
    
    for(int i=0;i<nbrStates*nbrLevels;i++,p_levels++)
    {
        //utilities::cout.DebuggingInfo()<<"energyLevel["<<i<<"] = "<<*(p_levels)<<std::endl;

        const double alpha = energyLevels[i];
        const double oneMinusAlpha = 1.0-alpha;
        
        double currentNormalization = log(oneMinusAlpha);

        *(p_levels) = log(alpha/oneMinusAlpha);
        
        //  inf or nan results can occur in the case where we calculated
        //  the overlap matrix at the end of the spectrum - these contain dummy states.
        //  OR they can occur due to numerical integration overflow
        
        if(std::isinf(*(p_levels)) || std::isnan(*(p_levels)))
        {
            //  in case of overflow error, set equal to - the value at the other 
            //  end of the spectrum (assuming symmetry)
        
            *(p_levels) = -newLevels[nbrStates*nbrLevels-i-1];
        
            //  if both ends are inf/nan then it corresponds to a dummy state 
            //  e.g. in a lower Landau level at the end of the spectrum
            //  We can safely set the energy value to 0, and this will avoid propagation
            //  of inf and nan results (the corresponding state is never occupied)
        
            if(std::isinf(*(p_levels)) || std::isnan(*(p_levels)))
            {
                *(p_levels) = 0;
            }
        }

        if(std::isinf(currentNormalization) || std::isnan(currentNormalization))
        {
            //std::cerr<<"\n\tWARNING: log("<<oneMinusAlpha<<") is inf or nan for state "<<i<<std::endl;
        
            //  in case of overflow error, set equal to the value at the other 
            //  end of the spectrum (assuming symmetry)
        
            currentNormalization = log(energyLevels[nbrStates*nbrLevels-i-1]);
        
            if(std::isinf(currentNormalization) || std::isnan(currentNormalization))
            {
                //  If it's still inf or nan then we set the normalization contribution 
                //  to be zero to avoid propagation of errors
                
                currentNormalization = 0.0;
            }
        }
        
        //  Accumulate the normalization value
        
        normalization += currentNormalization;

        utilities::cout.DebuggingInfo()<<"\tSingle particle pseudoenergy ["<<i<<"] = "<<*(p_levels)<<std::endl;
    }
    
    //  copy the new levels to overwrite the energy levels buffer
    
    memcpy(energyLevels,newLevels,nbrStates*nbrLevels*sizeof(double));
    
    //  clear up temporary memory
    
    delete[] newLevels;

    utilities::cout.DebuggingInfo()<<"\n\tNormalization = "<<normalization<<std::endl;
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if _DEBUG_

////////////////////////////////////////////////////////////////////////////////
//! \brief A debugging function that can be used to check the angular momentum 
//! and counting of the generated configurations
//!
////////////////////////////////////////////////////////////////////////////////

void IqheRses::CheckAngularMomentum(
    std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,   //!<    List of LL occupations
    const double monopoleStrength,//!<   Monopole strength
    const int nbrA,           //!<    Number in the A region
    const int nbrLevels,      //!<    Number of Landau levels
    const int sector)         //!<    Given angular momentum sector
    const
{
    //  Generate a full list of angular momentum labels for the highest Landau level
    
    //PRINT("sector",sector);
    
    int mSize = 1+(monopoleStrength+nbrLevels-1)*2;
    
    std::vector<double> mList(mSize);
    
    double m = -monopoleStrength-(nbrLevels-1);
    
    for(int i=0;i<mSize;i++,m+=1.0)
    {
        mList[i] = m;
    }
    
    PRINT("mList",&mList[0],mSize);

    //  Check the angular momentum eigenstates represented by each occupationsList

    for(unsigned int i=0;i<occupationsList.size();i++)
    {
        double mCheck = 0;
    
        for(int j=0;j<nbrLevels;j++)
        {
            for(int k=0;k<mSize;k++)
            {
                //PRINT("occupationsList[i][j][k]",occupationsList[i][j][k]);
                
                //PRINT("mList[k]",mList[k]);
            
                mCheck += occupationsList[i][j][k]*mList[k];
                
                //getchar();
            }
        }
        
        PRINT("Angular momentum check: ",mCheck);
    }

    return;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function generates a file name string to identify output files
//!
////////////////////////////////////////////////////////////////////////////////

std::string IqheRses::GenFileName(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int sector)                                   //!<    Angular momentum
    const
{
	std::stringstream fileName;	//	string stream to build file name

	fileName.str("");

	fileName<<(*optionList)["out-path"].as<std::string>()<<"/Exact_IQHE_Eigenvalues_n_"
	<<(*optionList)["nbr-n"].as<int>()<<"_na_"<<(*optionList)["nbr-a"].as<int>()
	<<"_LL_"<<(*optionList)["nbr-ll"].as<int>()<<".dat";

    return fileName.str();
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Assign the energy levels calculated by the function GenerateNewLevels()
//! to the set of orbital occupations generated by the function GenerateLevelOccupations
//!
////////////////////////////////////////////////////////////////////////////////

void IqheRses::SpectrumToFile(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    double* spectrum,           //!<    Sorted list of entanglement eigenvalues
    //std::vector<int>& spectrumLabels,//!<    List of lables identifying the RSES branches
    const int spectrumDim,      //!<    Length of the spectrum array
    const int sector)           //!<    Angular momentum sector labelling the current set of data
    const
{
    std::string fileName = this->GenFileName(optionList,sector);
        
    std::ofstream f_out;
    
    if(sector == (*optionList)["min-sector"].as<int>())
    {
        f_out.open(fileName.c_str());
    }
    else
    {
        //  Append to existing file
        
        f_out.open(fileName.c_str(),std::ios::app);
    }

    if(!f_out.is_open())
    {
	    std::cerr<<"\tFile: "<<fileName<<" not found !!! "<<std::endl<<std::endl;
	    
	    exit(EXIT_FAILURE);
    }
    
    //  Sort the spectrum in descending order
    
    std::sort(spectrum,spectrum+spectrumDim);

    for(int i=0;i<spectrumDim;i++)
    {
        f_out<<sector<<"\t";
        f_out<<spectrum[i]<<"\t";
	    f_out<<exp(-spectrum[i])<<"\n";
    }

    f_out.close();

    utilities::cout.MainOutput()<<"\n\t- SECTOR "<<sector<<" DATA WRITTEN TO FILE: "<<fileName<<std::endl<<std::endl;

    utilities::cout.SecondaryOutput()<<"\n\tPseudoenergy"<<"\t\tSector"<<std::endl;

    for(int i=0;i<spectrumDim;i++)
    {
        utilities::cout.SecondaryOutput()<<"\t"<<spectrum[i]<<"\t\t\t"<<sector<<std::endl;
    }
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
    
