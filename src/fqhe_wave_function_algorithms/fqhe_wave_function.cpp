////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This file defines some common functions called by all classes derived from 
//!		this base wave function class. 
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

#include "fqhe_wave_function.hpp"

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Default constructor for the wave function data struct
//!
//! Sets default parameter values
////////////////////////////////////////////////////////////////////////////////
FQHE::WaveFunctionData::WaveFunctionData()
    :   
    nbr(4),
    nbrUp(4),
    nbrDown(0),
    nbrQh(0),
    jastrowExponent(1),
    nbrClusters(3),
    fillNumerator(1),
    fillDenominator(1),
    fillingFactor(1),
    flux(4),
    wfName("Laughlin"),
    wfFileName("out"),
    path("./"),
    fusionChannels(1),
    geometry(FQHE::_SPHERE_),
    type(FQHE::_LAUGHLIN_),
    statistics(FQHE::_BOSONS_),
    monopoleStrength(1),
    shift(0)
{
    #if _ENABLE_TORUS_GEOMETRY_
    this->torusState = 0;
    this->torusDegeneracy = 1;
    #endif
    #if _ENABLE_DISC_GEOMETRY_
    this->discRadius = 1;	
    #endif
    #if _ENABLE_LATTICE_FQHE_
    this->latticeSpacing = 1;
    this->Lx = 1;
    this->Ly = 1;
    this->maxOccupation = 1;
    #endif
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Sets parameter values using command line arguments
//!
//! \return true if there was an error , false otherwise
////////////////////////////////////////////////////////////////////////////////
void FQHE::WaveFunctionData::InitFromCommandLine(
#if _ENABLE_MPI_
    boost::program_options::variables_map* options, //!<    Command line argument list
    const utilities::MpiWrapper& mpi)               //!<    instance of the mpi wrapper class
#else
    boost::program_options::variables_map* options) //!<    Command line argument list
#endif
{
    std::string tempType;
    std::string tempGeometry;
    bool exitFlag = false;
    #if _ENABLE_MPI_
    if(0 == mpi.m_id)    //  For the master node
    #endif
    {
        this->nbr = (*options)["nbr"].as<int>();
        //this->nbrQh = (*options)["nbr-qh"].as<int>();
        this->jastrowExponent = (*options)["jastrow-exponent"].as<int>();
        this->nbrClusters = (*options)["cluster-size"].as<int>();
        if((*options)["type"].as<std::string>()=="l")
        {	
	        this->type = FQHE::_LAUGHLIN_;
        }
        else if((*options)["type"].as<std::string>()=="cf")
        {
	        this->type = FQHE::_COMPOSITE_FERMION_;
        }
        else if((*options)["type"].as<std::string>()=="bs")
        {	
	        this->type = FQHE::_BONDERSON_SLINGERLAND_;
        }
        else if((*options)["type"].as<std::string>()=="mr")
        {
            this->type = FQHE::_MOORE_READ_;
        }
        else if((*options)["type"].as<std::string>()=="nass")
        {
            this->type = FQHE::_NASS_;
        }
        else if((*options)["type"].as<std::string>()=="para")
        {
            this->type = FQHE::_PARAFERMION_;
        }
        else
        {
	        std::cerr << "\tERROR! unknown type " << tempType << std::endl;
	        exitFlag = true;
        }
        if((*options)["bosonic"].as<bool>())
        {
	        this->statistics = FQHE::_BOSONS_;
        }
        else
        {
	        this->statistics = FQHE::_FERMIONS_;
        }
        if((*options)["geometry"].as<std::string>()=="s")
        {
            this->geometry = FQHE::_SPHERE_;
        }
        else if((*options)["geometry"].as<std::string>()=="d")
        {
            this->geometry = FQHE::_DISC_;
        }
        else if((*options)["geometry"].as<std::string>()=="t")
        {
            this->geometry = FQHE::_TORUS_;
        }
        else
        {
	        std::cerr << "\tERROR! unknown geometry " << tempGeometry << std::endl;
	        exitFlag = true;
        }
        #if _ENABLE_LATTICE_FQHE_
        this->Lx = (*options)["torus-x"].as<int>();
        this->Ly = (*options)["torus-y"].as<int>();
        this->torusState = (*options)["torus-state"].as<int>();
        #endif
    }
    #if _ENABLE_MPI_
    mpi.Sync(&exitFlag, 1, 0);
    #endif
    if(exitFlag)
    {
        exit(EXIT_FAILURE);
    }
    #if _ENABLE_MPI_
    this->MpiSync(0, mpi);
    #endif
    return;
}

//!
//! Set the Sphere and/or disc radius based on the current value of the flux
//!
void FQHE::WaveFunctionData::InitRadius()
{   
    this->sphereRadius=sqrt((double)(this->flux/2.0));
    #if _ENABLE_DISC_GEOMETRY_
        this->discRadius=sqrt((double)(this->flux/2.0));
    #endif
}      

//!
//!	Generates a string containing a file name composed of relevant
//! wave function data e.g of the form bosons_sphere_n_10_flux_10
//!	as appropriate for the geometry
//!
void FQHE::WaveFunctionData::GenerateFileName()
{
    std::stringstream name;
    name.str("");
    if(this->statistics == _BOSONS_)    
    {
        name << "bosons_";
    }
    else				           	    
    {
        name << "fermions_";
    }
    #if _ENABLE_SPHERE_GEOMETRY_
    if(this->geometry == _SPHERE_)	
    {
        name << "sphere";
    }
    #endif
    #if _ENABLE_DISC_GEOMETRY_
    if(this->geometry == _DISC_)		
    {
        name << "disc";
    }
    #endif
    #if _ENABLE_TORUS_GEOMETRY_
    if(this->geometry == _TORUS_) 	
    {
        name << "torus" << "_" << this->torusState << "_" << this->Lx << "_" << this->Ly;
    }
    #endif
    name << "_n_" << this->nbr << "_flux_" << this->flux;
    wfFileName = name.str();
}

//!
//!	Prints out a full summary of current wave function data
//!
void FQHE::WaveFunctionData::CheckAndPrint() const
{
	utilities::cout.DebuggingInfo() << "\tPROGRAM COMPILED WITH THE FOLLOWING FEATURES:\n" << std::endl;
	#if _ENABLE_TORUS_GEOMETRY_
		utilities::cout.DebuggingInfo() << "\t\tTORUS GEOMETRY ENABLED\n" << std::endl;
		#if _WITH_GENERALISED_THETA_FUNCS_
		utilities::cout.DebuggingInfo() << "\t\t(COMPILED GENERALISED THETA FUNCTIONS FORMULATION)\n" << std::endl;
		#endif
		#if _WITH_COM_ZEROES
		utilities::cout.DebuggingInfo() << "\t\t(COMPILED COM ZEROS FORMULATION)\n" << std::endl;
		#endif
	#endif
	#if _ENABLE_DISC_GEOMETRY_
		utilities::cout.DebuggingInfo() << "\t\tDISC GEOMETRY ENABLED\n" << std::endl;
	#endif
	#if _ENABLE_SPHERE_GEOMETRY_
		utilities::cout.DebuggingInfo() << "\t\tSPHERE GEOMETRY ENABLED\n" << std::endl;
	#endif
	utilities::cout.MainOutput() << utilities::cout.HyphenLine() << std::endl 
	                             << utilities::cout.HyphenLine() << std::endl << std::endl;
    utilities::cout.MainOutput() << "\tWorking with the "<<this->wfName << " wave function:\n";
    utilities::cout.MainOutput() << "\n\tNo. electrons\t\t\t" << this->nbr;
    utilities::cout.MainOutput() << "\n\tTotal flux\t\t\t" << this->flux;
    utilities::cout.MainOutput() << "\n\tFilling Factor\t\t\t" << this->fillNumerator 
                                 << "/" << this->fillDenominator;
    #if _ENABLE_SPHERE_GEOMETRY_
    if(this->geometry==_SPHERE_)
    {
	    utilities::cout.MainOutput() << "\n\tUsing sphere geometry" << std::endl;
	    utilities::cout.MainOutput() << "\n\tShift (sphere geometry)\t\t" << this->shift;
	    utilities::cout.MainOutput() << "\n\tMonopole Strength\t\t" << this->monopoleStrength;
    }
    #endif
    #if _ENABLE_TORUS_GEOMETRY_
    if(this->geometry==_TORUS_)
    {
	    utilities::cout.MainOutput() << "\n\tUsing torus geometry" << std::endl;
	    utilities::cout.MainOutput() << "\n\tGround state degeneracy\t\t" << this->torusDegeneracy;
	    utilities::cout.MainOutput() << "\n\tFlux per site \t" << (this->flux)/(this->Lx*this->Ly);
    }
    #endif
    #if _ENABLE_DISC_GEOMETRY_
    if(this->geometry==_DISC_)
    {
	    utilities::cout.MainOutput() << "\n\tUsing sphere geometry" << std::endl;
	    utilities::cout.MainOutput() << "\n\tDisc radius\t\t" << this->discRadius;
    }
    #endif
    utilities::cout.MainOutput() << "\n\tStatistics\t\t\t" 
        << ((this->statistics == _BOSONS_) ? ("Bosons") : ("Fermions")) << std::endl << std::endl;
    return;
}

#if _ENABLE_MPI_                      
//!
//!	MPI - communicate values required on parallel processes to the other
//! nodes
//!                                  
void FQHE::WaveFunctionData::MpiSync(
    const int syncId,                   //!<    ID of node to sync with
    const utilities::MpiWrapper& mpi)   //!<    instance of the mpi wrapper class
{
    mpi.Sync(&this->nbr, 1, syncId);
    mpi.Sync(&this->nbrUp, 1, syncId);
    mpi.Sync(&this->nbrDown, 1, syncId);
    mpi.Sync(&this->nbrQh, 1, syncId);
    mpi.Sync(&this->jastrowExponent, 1, syncId);
    mpi.Sync(&this->nbrClusters, 1, syncId);
    mpi.Sync(&this->fillNumerator, 1, syncId);
    mpi.Sync(&this->fillDenominator, 1, syncId);
    mpi.Sync(&this->fillingFactor, 1, syncId);
    mpi.Sync(&this->flux, 1, syncId);
    mpi.Sync(this->wfName, syncId);
    mpi.Sync(this->wfFileName, syncId);
    mpi.Sync(this->path, syncId);
    mpi.Sync(&this->fusionChannels, 1, syncId);
    mpi.Sync((int*)&this->geometry, 1, syncId);
    mpi.Sync((int*)&this->type, 1, syncId);
    mpi.Sync((int*)&this->statistics, 1, syncId);
    mpi.Sync(&this->monopoleStrength, 1, syncId);
    mpi.Sync(&this->shift, 1, syncId);
    mpi.Sync(&this->sphereRadius, 1, syncId);
    #if _ENABLE_TORUS_GEOMETRY_
    mpi.Sync(&this->torusState, 1, syncId);
    mpi.Sync(&this->torusDegeneracy, 1, syncId);
    #endif
    #if _ENABLE_DISC_GEOMETRY_
    mpi.Sync(&this->discRadius, 1, syncId);
    #endif
    #if _ENABLE_LATTICE_FQHE_
    mpi.Sync(&this->latticeSpacing, 1, syncId);
    mpi.Sync(&this->Lx, 1, syncId);
    mpi.Sync(&this->Ly, 1, syncId);
    mpi.Sync(&this->maxOccupation, 1, syncId);
    #endif
}
#endif

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the WaveFunction base class. This is called 
//!	immediately before a call to the constructor of any derived class.
//!
//!	It simply points the address of the internal wave function data structure
//!	at the externally set data set passes to the class.
////////////////////////////////////////////////////////////////////////////////
FQHE::WaveFunction::WaveFunction(
	FQHE::WaveFunctionData *wavefunctionData)		//!<	Address of a struct of wave function data
	:	m_wfData(wavefunctionData)
{}

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the WaveFunction base class. This is called 
//!	immediately after a call to the destructor of any derived class.
//!
//!	Declared empty.
////////////////////////////////////////////////////////////////////////////////
FQHE::WaveFunction::~WaveFunction()
{}

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Get the particle number
//!
//! \return Number of particles
////////////////////////////////////////////////////////////////////////////////
int FQHE::WaveFunction::GetNbrParticles() const
{
    return m_wfData->nbr;
}

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Get the monopole strength
//!
//! \return Monopole strength
////////////////////////////////////////////////////////////////////////////////
double FQHE::WaveFunction::GetMonopoleStrength() const
{
    return m_wfData->monopoleStrength;
}
