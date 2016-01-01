////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 30/12/2014
//!
//!  \file 
//!		This file contains a declaration of the program options and data
//!     structures to hold option parameters for the entanglement spectrum
//!     code fqhe_entanglement_spectrum
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

#ifndef _PROGRAM_OPTIONS_HPP_INCLUDED_
#define _PROGRAM_OPTIONS_HPP_INCLUDED_

///////		LIBRARY INCLUSIONS		////////////////////////////////////////////

#include <boost/program_options.hpp>        //	Include program options library
#include "../utilities/wrappers/mpi_wrapper.hpp"    
                                            //  Wrapper for MPI functionality

#if _DEBUG_
#include "../utilities/general/debug.hpp"
#endif

//////      DECLARATION OF PROGRAM OPTIONS      ////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a list of entanglement spectrum program options
//!
////////////////////////////////////////////////////////////////////////////////

inline boost::program_options::options_description GetEntanglementSpectrumOptions()
{
    boost::program_options::options_description entanglementSpectrumOpt("Entanglement Spectrum Options");
    entanglementSpectrumOpt.add_options()
    ("path",boost::program_options::value<std::string>()->default_value("./"),
     "Specify the name of the parameters.path where program data are stored\n")
    ("lza",boost::program_options::value<int>()->default_value(3000),
     "An integer label denoting the lowest orbital angular momentum sector used\n")
    ("real-cut",boost::program_options::value<int>()->default_value(22),
     "Set the real space cut: to confine region A to a disk set by l_A<realCut [set to be roughly flux/2]\n")
    ("nbr-fixed",boost::program_options::value<int>()->default_value(0),
    "Set the number of particles in region A.\n")
    ("nbr-sect",boost::program_options::value<int>()->default_value(20),
     "Set the number of angualr momentum sectors to calculate for (starting from the sector set with lzA)\n")
    ("nbr-fourier",boost::program_options::value<int>()->default_value(300),
     "Set the number of Fourier component to use\n")
    ("rows",boost::program_options::value<int>()->default_value(128),
     "Set the number of rows in the reduced density matrix\n")
    ("columns",boost::program_options::value<int>()->default_value(256),
     "Set the number of columns in the reduced density matrix\n")
    ("renorm",boost::program_options::value<int>()->default_value(0),
     "Specfy x, where we renormalise the wave function value by a factor of 10^(x)\n");

    return entanglementSpectrumOpt;
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a list of metropolis sampling program options
//!
////////////////////////////////////////////////////////////////////////////////

inline boost::program_options::options_description GetMetropolisOptions()
{
    boost::program_options::options_description metropOptions("Metropolis Options");
    metropOptions.add_options()
    ("iter,i",boost::program_options::value<int>()->default_value(50),
     "Number of Metropolis steps between each sample\n")
    ("samples",boost::program_options::value<int>()->default_value(5000),
     "Number of Monte Carlo samples\n")
    ("therms,t",boost::program_options::value<int>()->default_value(5000),
     "Number of thermalizing Monte Carlo samples\n")
    ("maxd", boost::program_options::value<double>()->default_value(0.5),
    "Proportion of maximum MC move distance\n");

    return metropOptions;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////      DEFINE DATA STRUCTS TO HOLD PROGRAM OPTION PARAMETERS      /////////

////////////////////////////////////////////////////////////////////////////////
//! \brief This struct contains all parameters that are parsed from the 
//! command line option set. 
//!
////////////////////////////////////////////////////////////////////////////////

struct GeneralOptions
{
    std::string path;		//!<   	Specify path to file directory for input/output
	                        //!     files
	int choice;				//!<	Choose between different operational modes:
	                        //!	 - 1: generate 'indexes' 
	                        //!  - 2: compute density matrix and store in a file
	                        //!	 - 3: calculate SVD from file	of matrix values
	int lzA;				//!<	z component of angular momentum for group A
	                        //!
	int realCut;			//!<	The particle entanglement spectrum cut
	                        //!
	int nbrFixed;			//!<	No. particles in B group (traced over to get ES)
	                        //!
	int nbrSect;			//!<    Number of sectors in the spectrum
	                        //!
	int nbrFourier;			//!<    Number of Fourier components
	                        //!<    
	int rows;		        //!<	No. rows in pseudo-reduced density matrix
	                        //!
	int columns;            //!<	No. columns in pseudo-reduced density matrix
	                        //!
	int renormFactor;		//!<	10^Factor by which to renormalise the wave function
	                        //!     to avoid overflows
	
	////////////////////////////////////////////////////////////////////////////////
    //! \brief Set parameters from command line argument list
	//!
	////////////////////////////////////////////////////////////////////////////////
	                        
	inline void InitFromCommandLine(
	    boost::program_options::variables_map* options, //!<    Command line argument list
	    const utilities::MpiWrapper& mpi)   //!<    MPI wrapper class
	{
	    if(0 == mpi.m_id)    //  For the master node
	    {
	        path         = (*options)["path"].as<std::string>();
	        choice       = (*options)["mode"].as<int>();
	        lzA          = (*options)["lza"].as<int>();
	        realCut      = (*options)["real-cut"].as<int>();
	        nbrFixed     = (*options)["nbr-fixed"].as<int>();
	        nbrSect      = (*options)["nbr-sect"].as<int>();
	        nbrFourier   = (*options)["nbr-fourier"].as<int>();
	        rows         = (*options)["rows"].as<int>();
	        columns      = (*options)["columns"].as<int>();
	        renormFactor = (*options)["renorm"].as<int>();
        }
        
        //  Synchronize all values with those set on the master node
        
        this->MpiSync(0,mpi);
        
        if(0 == choice)
        {
            exit(EXIT_SUCCESS);
        }
        
	    return;
	}
	     
	////////////////////////////////////////////////////////////////////////////////
    //! \brief Function to synchronise all these variables
	//! with the value set on syncNode
	//!  
	////////////////////////////////////////////////////////////////////////////////
                                              
	inline void MpiSync(
	    const int syncNode,                 //!<    Node to sync with
	    const utilities::MpiWrapper& mpi)   //!<    MPI wrapper class
	{   
	    mpi.Sync(path,syncNode);
	    mpi.Sync(&choice,1,syncNode);
	    mpi.Sync(&lzA,1,syncNode);
	    mpi.Sync(&realCut,1,syncNode);
	    mpi.Sync(&nbrFixed,1,syncNode);
	    mpi.Sync(&nbrSect,1,syncNode);
	    mpi.Sync(&nbrFourier,1,syncNode);
	    mpi.Sync(&rows,1,syncNode);
	    mpi.Sync(&columns,1,syncNode);
	    mpi.Sync(&renormFactor,1,syncNode);
	}       
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief This struct contains parameters associated with the metropolis.
//! sampling part of the code (mode 1). These parameters are parsed from the
//! command line and are only required on the master node
//!
////////////////////////////////////////////////////////////////////////////////

struct MetropOptions
{
	double width;			//!<	Max displacement in metropolis. move
	                        //!
	int iter;				//!<	Number of Metropolis moves per sample
	                        //!
	int nbrSamples;			//!<	Number of Monte Carlo samples to take
	                        //!
	int nbrTherms;			//!<	Number of thermalizing samples
	                        //!<

	////////////////////////////////////////////////////////////////////////////////
    //! \brief Set parameters from command line argument list
	//!  
	////////////////////////////////////////////////////////////////////////////////
	                        
	inline void InitFromCommandLine(
	    boost::program_options::variables_map* options, //!<    Command line argument list
	    const utilities::MpiWrapper& mpi)   //!<    MPI wrapper class
	{
	    if(0 == mpi.m_id)    //  For the master node
	    {
	        width        = (*options)["maxd"].as<double>();
	        iter         = (*options)["iter"].as<int>();
	        nbrSamples   = (*options)["samples"].as<int>();
	        nbrTherms    = (*options)["therms"].as<int>();
        }

	    return;
	}                               
};

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#endif
