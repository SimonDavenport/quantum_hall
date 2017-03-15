////////////////////////////////////////////////////////////////////////////////
//!
//!                   Header file for fqheMonteCarloAnalyis.cpp
//!
//!						  Author: Simon C. Davenport
//!
//! \file
//!		This is the header file for the analyseMonteCarlo.cpp program
//!		Basic compilation will only include the file i/o functions
//!		to enable all analysis methods, compile with -DENABLE_ANALYSIS_METHODS
//!		This is done to save compilation time on projects that only need to use
//!		the file i.o and re-sampling functions
//!
//!					Copyright (C) Simon C Davenport
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

#ifndef _ANALYSIS_METHODS_HPP_INCLUDED_
#define _ANALYSIS_METHODS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstring>
#if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
#include "../../utilities/wrappers/murmur_hash_wrapper.hpp"
#endif
#if _SPEED_OPTIMIZED_MAP_
#include <sparsehash/dense_hash_map>
#elif _MEMORY_OPTIMIZED_MAP_
#include <sparsehash/sparse_hash_map>
#else
#include <unordered_map>
#endif
#include "../../fqhe_wave_function_algorithms/fqhe_wave_function.hpp"
#include "../../fqhe_wave_function_algorithms/laughlin.hpp"
#include "../../fqhe_wave_function_algorithms/moore_read.hpp"
#include "../../fqhe_wave_function_algorithms/nass.hpp"
#include "../../fqhe_wave_function_algorithms/parafermion.hpp"
#include "../../fqhe_wave_function_algorithms/composite_fermion.hpp"
#include <boost/program_options.hpp>
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/mathematics/mt.hpp"
#include "../../utilities/general/dcmplx_type_def.hpp"
#include "../../utilities/general/pi_const_def.hpp"
#if _ENABLE_MPI_
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#endif
static const int DISC_BG_FUNC_DIM  = 50000;
static const int THICKNESS_POTENTIAL_FILE_SIZE = 88000;
static const int THICKNESS_BG_SIZE = 1000;
static const int INIT_ARRAY_SIZE   = 100;
static const int SIZE_OF_DOUBLE    = sizeof(double);
static const int SIZE_OF_INT       = sizeof(int);

inline boost::program_options::options_description GetAnalysisOptions()
{
    namespace po = boost::program_options;

    po::options_description files_opt("File/Folder options");
	files_opt.add_options()
	("input",po::value<std::string>()->default_value("./"),
	 "Specify the name of the folder where input data are stored\n")
	("output",po::value<std::string>()->default_value("./"),
	 "Specify the name of the folder where output data are to be stored\n")
	 ("get-file-name",po::value<bool>()->default_value(false),
	 "Set this option to true in order to return the data file name, then exit\n")
	("skip",po::value<int>()->default_value(0),
	 "Skip the first arg samples before averaging")
	("check",po::value<bool>()->default_value(false),
	 "Set option to true in order to carefully check data for numerical problems");
	
	po::options_description analysis_opt("Basic Analysis options");	
	analysis_opt.add_options()
	("pair",po::value<bool>()->default_value(true),
	 "Option to calculate pair correlation function")
	("pair-bins",po::value<int>()->default_value(500),
	 "Number of bins to discretize pair correlation data")
	("dens",po::value<bool>()->default_value(true),
	 "Option to calculate density function")
	("dens-bins",po::value<int>()->default_value(500),
	 "Number of bins to discretize density data")
	("coulomb",po::value<bool>()->default_value(true),
	 "Analyse data with Coulomb interaction")
	("second",po::value<bool>()->default_value(false),
	 "Analyse data with 2nd LL type interaction")
	("second-type",po::value<std::string>()->default_value("d"),
	 "Select type of 2nd LL interaction to use \n\nALLOWED TYPES:\n"
	 "  d (Disc geometry derived effective potential)\n"
	 "	  s (Sphere geometry derived effective potential)\n"
	 "	  g (Graphene appropriate effectiuve potential)\n")
	("re-sample",po::value<bool>()->default_value(true),
		 "Set to false to skip the re-sampling procedure\n")
	("nbr-re-sample",po::value<int>()->default_value(1000),
		 "Set the number of re-samples to take\n");
	
	po::options_description finite_opt("Finite thickness correction options");	
	finite_opt.add_options()
	("finite",po::value<bool>()->default_value(false),
	 "Analyse data with finite thickness correction")
	("finite-type",po::value<std::string>()->default_value("h"),
	 "Select type of finite thickness potential to use \nALLOWED TYPES:\n"
	 "  h (Heterojunction/Fang-Howard)\n"
	 "	  w (Square well)\n")
	("min-thickness",po::value<double>()->default_value(0.01),
	 "Select minimum value of the thickness parameter")
	("max-thickness",po::value<double>()->default_value(5),
	 "Select maximum value of the thickness parameter")
	("step-thickness",po::value<double>()->default_value(0.01),
	 "Select thickness parameter step size between min and max");
	
	po::options_description bilayer_opt("Bilayer potential options");	
	bilayer_opt.add_options()
	("bilayer",po::value<bool>()->default_value(false),
	 "Analyse data with bilayer potential")
	("min-bilayer",po::value<double>()->default_value(0.01),
	 "Select minimum value of the bilayer parameter")
	("max-bilayer",po::value<double>()->default_value(5),
	 "Select maximum value of the bilayer parameter")
	("step-bilayer",po::value<double>()->default_value(0.01),
	 "Select bilayer parameter step size between min and max");
	
	po::options_description plate_opt("Charge-plate potential options");	
	plate_opt.add_options()
	("charge-plate",po::value<bool>()->default_value(false),
	 "Analyse data with charge plate (capacitor) potential")
	("min-plate",po::value<double>()->default_value(0.01),
	 "Select minimum value of the charge plate parameter")
	("max-plate",po::value<double>()->default_value(5),
	 "Select maximum value of the charge plate parameter")
	("step-plate",po::value<double>()->default_value(0.01),
	 "Select charge plate parameter step size between min and max")
	("plate-bilayer",po::value<bool>()->default_value(false),
	 "Analyse data with combined bilayer and charge plate potential. Note that the parameter range is set to bilayer thickness 0 to 10 in steps of 0.4 and charge plate separation 0.5 to 10.5  in steps of 0.5 ");
	
	po::options_description partial_opt("Partial Coulomb interaction");	
	partial_opt.add_options()
	("partial-coulomb",po::value<bool>()->default_value(false),
	 "Analyse data with only part of the Coulomb interaction")
	("partial-type",po::value<int>()->default_value(1),
	 "An integer code to select the partial type:\n\nALLOWED TYPES:\n"
	 "  1 - LLL Coulomb (up-up pairs only)\n"
	 "  2 - LLL Coulomb (up-down pairs only)\n"
	 "  3 - Disc 2LL Coulomb (up-up pairs only)\n"
	 "  4 - Disc 2LL Coulomb (up-down pairs only)\n" 
	 "  5 - Graphene 2LL Coulomb (up-up pairs only)\n"
	 "  6 - Graphene 2LL Coulomb (up-down pairs only)\n" 
	 "  7--13 - Effective potential (r)^(arg-6)*e^-r (up-up pairs only)\n"
	 "  14--20 - Effective potential (r)^(arg-13)*e^-r (up-down pairs only)\n");

	po::options_description all("");
	all.add(files_opt).add(analysis_opt).add(finite_opt).add(bilayer_opt).add(plate_opt).add(partial_opt);
	return all;
}

struct AnalysisOptions
{
	bool pairCorrel;            //!<	Option to calculate pair correlation function
	bool dens;                  //!<	Option to calculate density function
	bool coulomb;               //!<    Analyse data with Coulomb interaction
	bool secllCoulomb;          //!<	Analyse data with 2nd LL type interaction	    
	bool finiteThickness;       //!<	Analyse data with finite thickness correrction
	bool bilayer;               //!<    Use bilayer type interaction
	bool chargePlate;           //!<    Use charge plate type interaction
	bool biChargePlate;         //!<	Combine bilayer and charge plate potentials
	bool partialCoulomb;        //!<	Analyse data with only part of the Coulomb interaction
	bool resample;              //!<	Option to resample energy data or not
	int nbrResamples;           //!<    Option to set the number of energy resamples
	std::string dataDirIn;      //!<	String containing folder name for input data
	std::string dataDirOut;     //!<	String containing folder name for output data
	int skip;                   //!<    Skip the first 'skip' samples before averaging
	bool check;                 //!<	Set option to carefully check data for numerical errors
	int nPairBins;              //!<    Number of bins to discretize pair correlation data
	int nDensBins;              //!<	Number of bins to discretize density data
	std::string secllType;      //!<	Select type of 2nd LL interaction to use 
	                            //!<    (disc - d, graphene - g, sphere - s)
	std::string finiteType;     //!<	Select finite thickness potential type
	                            //!<    h for heterojunction, w for square well
	double minThickness;        //!<	Minimum value of thickness parameter to use
	double maxThickness;        //!<	Maximum value of thickness parameter to use
	double thicknessStepSize;   //!<	Step size of thickness parameter between 
	                            //!     minThickness and maxThickness
	double minBilayer;          //!<	Maximum value of bilayer parameter to use
	double maxBilayer;          //!<	Maximum value of bilayer parameter to use
	double bilayerStepSize;     //!<	Step size of bilayer parameter between 
	                            //!     minBilayer and maxBilayer
	double minPlate;            //!<	Maximum value of charge plate parameter to use
	double maxPlate;            //!<	Minimum value of charge plate parameter to use
	double plateStepSize;       //!<	Step size of charge plate parameter between 
	                            //!<    minPlate and maxPlate
	int partialCoulombType;     //!<	An integer to select which part of the coulomb 
	                            //!     interaction to work with
	bool getFileName;           //!<	Option to print out only the file names, 
                                //!<    then exit
	AnalysisOptions();
	void InitFromCommandLine(boost::program_options::variables_map* options);
};

class AnalysisMethods
{
	private:
	#if _ENABLE_ANALYSIS_METHODS_
	double *pairFreqUp;     //!<    Array of sample frequencies for
                            //!     discretised pair chord distances - up-up data
	double *pairFreqDown;   //!<    Array of sample frequencies for
                            //!     discretised pair chord distances - down-down data
	double *pairFreqDiff;   //!<    Array of sample frequencies for
                            //!     discretised pair chord distances - up-down data
	double *densityFreq;    //!<    Array of sample frequencies for
                            //!     discretised densoty function chord distances
	std::ofstream f_pair;   //!<    File stream for pair data
	std::ofstream f_dens;   //!<    File stream for density data
	std::ifstream f_background;
	                        //!<    File stream to file containing values of 2nd LL
                            //!<    background potential for different values of flux
	double *discBgFunc;     //!<    Buffer to store disc background energy function for
                            //!<    different values of chord distance
	double *thicknessBg;    //!<    Buffer to store finite thickness potential background
                            //!<    energy function for different values of chord distance
	double secondBg;        //!<    Background value for the second LL potential
	double *biChargePlateBg;//!<    Buffer to store background energy values for the
                            //!<    combined bilayer and charge plate potentials, as a 
                            //!<    function ofchord distance
	double partialBg;       //!<    Background energy value for the partial Coulomb
                            //!     interaction
	double *finiteThicknessPotential;
	                        //!<    Buffer to store discretized values of the finite 
                            //!     thickness potential function
	std::ifstream f_potential;
	                        //!<    In file stream for finite thickness potential function
	std::ifstream f_coeffs; //!<    File stream pointer for coefficients in the partial
                            //!     Coulomb interaction
	double *coeffList;      //!<	Array to store coefficients of 2nd LL 
	                        //!     effective potential
	std::string partialType;//!<	String to store info on the partial Coulomb 
	                        //!     interaction potential type
	std::ofstream f_energy; //!<    Filestream to store calcualted energy value
	std::ofstream f_second; //!<    Filestream to store calculated 2nd LL energy value
	std::ofstream f_thickness;
	                        //!<    Filestream to store calculated finite thickness
	                        //!     interaction energy
	std::ofstream f_bilayer;//!<    Filestream to store calculated bilayer interaction 
	                        //!     energy
	std::ofstream f_plate;  //!<    Filestream to store charge place interaction energy 
	std::ofstream f_biChargePlate;
	                        //!<    Filestream to store combined charge plate and bilayer
	                        //!     energy
	std::ofstream f_partialCoulomb;
	                        //!<    Filestream to store calculated partial Coulomb 
	                        //!     interaction energy
	#endif
	public:
	#if _ENABLE_ANALYSIS_METHODS_
	AnalysisOptions options;
	double *coulombEnergy;  //!<    Array to store all values of the energy
	double *secondEnergy;   //!<	Array to store all values of the second LL energy
	double *partialEnergy;  //!<	Array to store all values of the partial coulomb energy
	double *thicknessEnergy;//!<	Store cumulative values of the energy for the
					        //!	    finite-thickness type interaction
	double *bilayerEnergy;  //!<	Array to store cumulative values of bilayer energy
	double *chargePlateEnergy;
	                        //!<    Store cumulative values of the energy for the
					        //!	    charge plate type interactions
	double *biChargePlateEnergy;
	                        //!<	Array to store cumulative values of the combined 
					        //!	    bilayer-charge plate interaction energy
	double coulombStdev;    //!<    Error estimate for coulomb energy
	double secondStdev;     //!<    Error estimate for 2nd LL coulomb energy
	double partialStdev;    //!<    Error estimate for the partial Coulomb interaction
	double coulombMeanEnergy;
	                        //!<    Mean Coulomb energy
	double secondMeanEnergy;//!<    Mean 2nd LL Coulomb energy
	double partialMeanEnergy;
	                        //!<    Mean partial Coulomb energy
	#endif
	dcmplx *u;              //!<    Array to store spinor coordinates u (sphere geometry)
	dcmplx *v;              //!<    Array to store spinor coordinates v (sphere geometry)
	dcmplx *z;              //!<    Array to store disc geometry coordinate
	double *phi;            //!<    Array to store angle coordinates phi
	double *theta;          //!<    Array to store angle coordinates theta
	std::complex<int> *latticeConfig;
	                        //!<    Array to store lattice configuration
	std::ifstream f_theta;  //!<    In file stream to read in theta configuration data
	std::ifstream f_phi;    //!<    In file stream to read in phi configuration data
	std::ifstream f_r;      //!<    In file stream to read in r configuration data (disc
                            //!     geometry only)
	std::ifstream f_config; //!<    In file stream to read in lattice configuration data
	int maxRuns;            //!<    Value to store the maximum number of configuration
                            //!     files detected
	long int nbrSamples;    //!<    Total number of samples to be analysed
	FQHE::WaveFunction *wf; //!<    Object wrapping wave function routines
	FQHE::WaveFunctionData wfData;
	                        //!<    Object holding wave function metadata
    FQHE::CompositeFermionData cfData;
                            //!<    Object holding composite fermion wave function metadata
	AnalysisMethods();
	AnalysisMethods(AnalysisOptions& options);
	~AnalysisMethods();
	void InitFromCommandLine(	
	#if _ENABLE_MPI_
        boost::program_options::variables_map* options,
        const utilities::MpiWrapper& mpi);
    #else
        boost::program_options::variables_map* options);
    #endif
	int  OpenConfigurationDataFile(std::string,int);
	void GetConfigFromFileSphere();
	void GetConfigFromFileDisc();
	void GetConfigFromFileTorus();
	void PolarToSpinor();
	void CheckConfig(int);
	void PrintName();
	int  LookForConfigurationDataFiles();	
	void Resample(double*,int,double&,double&);
	#if _ENABLE_ANALYSIS_METHODS_
	void InitializeMethods();
	double CoulombEnergy(int, dcmplx*, dcmplx*);
	double CoulombEnergy(int, dcmplx*);
	double BilayerEnergy(int, int,double ,dcmplx*,dcmplx*);
	double ChargePlateEnergy(int, double, dcmplx*,dcmplx*);
	void PairCorrelation(int,int,dcmplx*,dcmplx*);
	void Density(int,dcmplx*);
	double SecondLLEnergy(int,dcmplx*,dcmplx*);
	double FiniteThicknessEnergy(int, double, dcmplx*,dcmplx*);
	double PartialCoulombEnergy(int,int, int, dcmplx*,dcmplx*);
	double BiChargePlateEnergy(int,int, int, int,int, dcmplx*,dcmplx *);
	void MeanEnergy(double*,double&,double&);
	void FinaliseMethods();
	#endif
};
#endif
