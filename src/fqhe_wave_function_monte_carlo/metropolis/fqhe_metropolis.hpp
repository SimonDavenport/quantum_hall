////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!		This is the header file for the Metropolis class.
//!
//!	 \todo
//!		Construct derived classes for different type of geometry.
//!		Add pre-compiler commands to enable/disable different parts of 
//!		the compilation based on the geometry required. 
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

#ifndef _FQHE_METROPOLIS_HPP_INCLUDED_
#define _FQHE_METROPOLIS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../utilities/mathematics/mt.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <complex>
#include <vector>
#include "../../fqhe_wave_function_algorithms/fqhe_wave_function.hpp"
#include <boost/program_options.hpp>
#ifdef _ENABLE_MPI_
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#endif
#ifdef _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

namespace FQHE
{
    ///////		STATIC CONSTANT DECLARATIONS		    ////////////////////////////
    static const int SIZE_OF_DOUBLE=sizeof(double);	//!<	Size of double type
    static const int SIZE_OF_INT=sizeof(int);		//!<	Size of int type
    static const int BG_FUNC_DIM=50000;				//!<	Size of disc background function
												    //!		stored in a file

    //!
    //! Generate a list of Metropolis monte carlo program options
    //!
    inline boost::program_options::options_description GetMetropolisOptions()
    {
        boost::program_options::options_description monteCarloOpt("Metropolis Monte Carlo options");
	    monteCarloOpt.add_options()
	    ("nbr-therm,t", boost::program_options::value<unsigned int>()->default_value(100000),
	     "Number of MC steps per particle for thermalisation\n")
	    ("nbr-samples,m", boost::program_options::value<unsigned int>()->default_value(100000),
	     "Number of configuration samples to take\n")
	    ("run", boost::program_options::value<unsigned int>()->default_value(1),
	     "Index of the MC run e.g run 1 produces files indexed with run_1\n")
	    ("maxd", boost::program_options::value<double>()->default_value(1),
	     "Proportion of maximum MC move distance\n")
	    ("sample-freq", boost::program_options::value<unsigned int>()->default_value(2),
	     "Set sampling frequency (in units of nbr particles)\n")
	    ("nbr-move", boost::program_options::value<unsigned int>()->default_value(1),
	     "Number of particles to move at each MC step\n")
	    ("wave",boost::program_options::value<bool>()->default_value(true),
	     "Set to true to store complex values of the wavefunction with each sample\n")
	    ("use-init-file",boost::program_options::value<bool>()->default_value(false),
	     "Set to true to read in initial configurtation from a file.")
	    ("init-file",boost::program_options::value<std::string>()->default_value("init"),
	     "Begin the simulation with a stored particle configuration in this path: \n"
	     "\t  For sphere case require *_phi.dat and *_theta.dat files\n" 
	     "\t  For disc case require *_theta.dat and *_r.dat files\n")
	    ("path",boost::program_options::value<std::string>()->default_value("./"),
	     "Specify the name of the path where data are to be stored\n");
	     return monteCarloOpt;
    }

    //!
    //!	Data structure type to contain Monte Carlo simulation information
    //!
    struct MonteCarloData
    {
        unsigned int nbrTherm;     	//!<    Number of thermalising configurations 
								    //!
        long int nbrThermSteps;	   	//!<	Number of thermalising steps
								    //!
        unsigned int nbrSamples;   	//!< 	Total number of sample configurations to
								    //!	 	take
        long int nbrSteps;		   	//!<	Number of Monte Carlo steps to perform
								    //!
        unsigned int runLabel;     	//!<	A number labelling the parallel MC run
								    //!
        double maxDist;            	//!<	Maximum size of MC moves (in units of sphere
								    //!   	or disc radius or lattice spacing 
        int maxTorusDist;		   	//!<	1+Maximum number of lattice spacings for
								    //!		lattice moves
        unsigned int sampleFreq;   	//!< 	Number of MC moves to perform between sampling
								    //!
        unsigned int nbrToMove;     //!<	Number of moves per Monte Carlo step
								    //!
        bool wave;                 	//!<	Option to store wave function values in a file
								    //!
        std::string	folderName;	   	//!<	Name of Folder where data are to be stored
								    //!
	    bool useInitFile;           //!<    Option to use initial configuration from a file
	                                //! 
	    std::string	initFileName;  	//!	    Name of the path containing initial	
								    //!	    configuration files
	    std::string path;           //!<    Path to output data
	                                //!      							
	    void InitFromCommandLine(
	    #if _ENABLE_MPI_
            boost::program_options::variables_map* options,
            const utilities::MpiWrapper& mpi);
        #else
            boost::program_options::variables_map* options);
        #endif               
	    #if _ENABLE_MPI_                                             
	    void MpiSync(int syncId,const utilities::MpiWrapper& mpi);
        #endif	
    };

    //!
    //!	This class contains functions to perform Metropolis Monte Carlo
    //!	moves for fractional quantum Hall states in various geometries and to 
    //!	calculate and store associated sample data. 
    //!
    class Metropolis
    {
        private:
	    std::ofstream f_wave;	        //!<    To store values of the wavefunction
	    std::ofstream f_theta,f_phi;    //!<    To store spherical polar co-ordinates
	    std::ofstream f_r;		        //!<    To store radial co-ordinate 
	                                    //!     (disc geometry)
	    std::ofstream f_latticeConfig;  //!<    To store torus lattice configuration
	    std::ifstream f_init;	        //!<    Initial configuration file
	    std::ifstream f_elbg;	        //!<	Disc geometry background energy function
	                                    //!<    wave function to a file	
	    std::stringstream fileName;	    //!<    To store name of file to be opened					
	    std::string phiFile;            //!<    Name of spherical polar phi data file
	    std::string thetaFile;          //!<    Name of spherical theta phi data file
	    std::string waveFile;           //!<    Name of wave function amplitude data file
	    std::string rFile;              //!<    Name of circular polar r data file
	    std::string latticeConfigFile;  //!<    Name of lattice configuration data file
	    MersenneTwister mt;		        //!<    Object in Mersenne class
	    double *theta;                  //!<    Address to store spherical polar theta
	    double *phi;		            //!<    Address to store spherical polar phi
	    double *r;				        //!<    Address to store circular polar r
	    double *thetaNew;               //!<    Address to store update for theta values
	    double *phiNew;                 //!<    Address to store update for phi values
	    double *bgFunc;			        //!<    Disc geometry background energy function
	    WaveFunctionData *wfData;       //!<    Object to store wave function metadata
	    MonteCarloData *mcData;         //!<    Object to store composite fermion wave
	                                    //!     function metadata
        public:
	    dcmplx *u;              //!<    Address to store spinor u co-ordinates
	    dcmplx *v;			    //!<    Address to store spinor v co-ordinates
	    dcmplx *uNew;           //!<    Address to store updated spinor u co-ordinates
	    dcmplx *vNew;		    //!<    Address to store updated spinor v co-ordinates
	    dcmplx *uTmp;           //!<    Address to store previous spinor u co-ordinates
	    dcmplx *vTmp;		    //!<    Address to store previous spinor v co-ordinates
	    dcmplx *z;				//!<    Address to store ratio of v to u co-ordiantes
	    dcmplx *zNew;           //!<    Address to store new ratio of v to u co-ordiantes
	    dcmplx *zTmp;		    //!<    Address to store previous ratio of v to u co-ordiantes
	    std::complex<int> *latticeConfig;
	                            //!<    To store torus lattice configuration
	    std::complex<int> *latticeConfigTmp;
	                            //!<    Temporary store of moved lattice co-ordinates
	    std::complex<int> *latticeConfigNew;
	                            //!<    New set of lattice co-ordinates
	    dcmplx wfValOld;		//!<    Value of the wave function from the last accepted move
	    dcmplx wfValNew;		//!<    New trial value of the wave function
	    long int acceptCount;	//!<    Count the number of MC moves which are accepted	
	    int *toMove;			//!<    List of the particles to be moved
	    double meanEnergy;		//!<    Mean coulomb energy
	    double stdDevEnergy;	//!<    Standard deviation of coulomb energy (estimate)	
	    double accumEnergy;		//!<    Keep track of total calculate energy (for averaging)
	    double accumEnergySquared;//!<  Keep track of total calculate energy^2
	    int energyCounter;		//!<    Keep track of how many times the energy is calculated
	    Metropolis(WaveFunctionData*, MonteCarloData*);
	    ~Metropolis();
	    void RandomMoveSphere();
	    void RandomMoveDisc();
	    void RandomMoveTorus();
	    void MetropolisTestSphere();
	    void MetropolisTestDisc();
	    void MetropolisTestTorus();
	    void ConfigurationToFileSphere();
	    void ConfigurationToFileDisc();
	    void ConfigurationToFileTorus();
	    double CoulombEnergy(int, dcmplx*, dcmplx*);
	    double CoulombEnergy(int, dcmplx*);
	    void CalculateRunningMean();
	    void PrintConfigurationSphere();
	    void PrintConfigurationDisc();
	    void PrintConfigurationTorus();
	    bool CheckOccupancy(std::complex<int>*);
    };
}
#endif
