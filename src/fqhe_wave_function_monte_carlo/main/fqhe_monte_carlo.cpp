////////////////////////////////////////////////////////////////////////////////
//!
//!					\author Simon C Davenport
//!
//!  \mainpage 
//!		The purpose of this program is to generate statistically random sample
//!		configurations of particle co-ordinates, based on the probability
//!		distribution |psi|^2 associated with a given ground state fractional 
//!		quantum Hall wave function.	
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

///////		LIBRARY INCLUSIONS		////////////////////////////////////////////
#include "../metropolis/fqhe_metropolis.hpp"
#include "../../fqhe_wave_function_algorithms/fqhe_wave_function.hpp"
#include "../../fqhe_wave_function_algorithms/laughlin.hpp"
#include "../../fqhe_wave_function_algorithms/moore_read.hpp"
#include "../../fqhe_wave_function_algorithms/nass.hpp"
#include "../../fqhe_wave_function_algorithms/composite_fermion.hpp"
#include "../../fqhe_wave_function_algorithms/parafermion.hpp"
#include <boost/program_options.hpp>
#include "../../utilities/mathematics/mt.hpp"
#include "../../utilities/general/cout_tools.hpp"
#if _ENABLE_MPI_  
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#endif
utilities::Cout utilities::cout;
///////		FUNCTION FORWARD DELCATATIONS		    ////////////////////////////
void PrintWelcomeMessage();
boost::program_options::variables_map ParseComandLine(int argc, char *argv[]);
#if _ENABLE_MPI_ 
utilities::MpiWrapper mpi(utilities::cout);
#endif

////////////////////////////////////////////////////////////////////////////////
//! \brief The main function performs a number of tasks:
//!     -   It calls a function to parse the command line and set the various 
//!			parameters stored in global data structs.
//!     -   It initializes a Monte Carlo configuration for the chosen
//!		    geometry and wave function type
//!		- 	It performs thermalization steps on the Monte Carlo configuration
//!			using the Metropolis method.
//!		- 	It gathers sample configurations and writes these to files, optionally
//!			along with a file containing the wave function values. 
//!
//!	With the Metropolis method employed here, sample configurations are taken from 
//!	the probability	distribution |psi|^2 associated with a given ground state fractional 
//!	quantum Hall wave function.	 
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
    #if _ENABLE_MPI_ 
	mpi.Init(argc, argv);
    #endif
	utilities::cout.MainOutput().precision(8);
	PrintWelcomeMessage();
    FQHE::WaveFunction *wf = 0;
    FQHE::WaveFunctionData wfData;
    FQHE::CompositeFermionData cfData;
    FQHE::Metropolis *mc;
    FQHE::MonteCarloData mcData;
    boost::program_options::variables_map optionList;
	optionList = ParseComandLine(argc, argv);
	#if _ENABLE_MPI_
	wfData.InitFromCommandLine(&optionList, mpi);
	cfData.InitFromCommandLine(&optionList, mpi);
    mcData.InitFromCommandLine(&optionList, mpi);
	#else
	wfData.InitFromCommandLine(&optionList);
	cfData.InitFromCommandLine(&optionList);
    mcData.InitFromCommandLine(&optionList);
	#endif
	if(wfData.type==FQHE::_LAUGHLIN_)
	{	
		wf = new FQHE::Laughlin(&wfData);
	}
	else if(wfData.type==FQHE::_MOORE_READ_)
	{	
		wf = new FQHE::MooreRead(&wfData);
	}
	else if(wfData.type==FQHE::_NASS_)
	{	
		wf = new  FQHE::NonAbelianSpinSinglet(&wfData);
	}
	else if(wfData.type==FQHE::_COMPOSITE_FERMION_ || wfData.type==FQHE::_BONDERSON_SLINGERLAND_)
	{
		wf = new FQHE::CompositeFermion(&wfData,&cfData);
	}
	else if(wfData.type==FQHE::_PARAFERMION_)
	{
	    wf = new FQHE::Parafermion(&wfData);
	}
	mc = new FQHE::Metropolis(&wfData, &mcData);
	//	Do first evaluation of the wave function
	switch (wfData.geometry)
	{
		case FQHE::_SPHERE_:	
			mc->wfValOld = wf->EvaluateWfSphere(wfData.nbr, mc->u, mc->v);
			break;
		case FQHE::_DISC_:
			mc->wfValOld = wf->EvaluateWfDisc(wfData.nbr, mc->z);
			break;
		case FQHE::_TORUS_:
			mc->wfValOld = wf->EvaluateWfTorus(wfData.nbr, mc->latticeConfig, wfData.torusState);
			break;		
	}			
	mc->accumEnergy = 0;
	mc->accumEnergySquared = 0;
	mc->energyCounter = 0;
	if(!mcData.useInitFile)		//	Run thermalization unless an initial configuration file is found
	{
		double timer;		    //	Keep track of the time taken during the simulation
		double energy=0;		//	Temporary energy value
		utilities::cout.MainOutput() << "\n\tThermalising system...\n\n" << utilities::cout.HyphenLine() 
		                             << std::endl;
		utilities::cout.SecondaryOutput() << "\tProgress\tRunning Average\t\tEstimated Time Remaining (hours)\n" 
		                                  << utilities::cout.HyphenLine() << std::endl;
		fflush(stdout);
		timer = clock();
		utilities::cout.SecondaryOutput().precision(6);
		for(long unsigned int k=1; k<=mcData.nbrTherm; ++k)
		{
			switch (wfData.geometry)
			{
				case FQHE::_SPHERE_:
					mc->RandomMoveSphere();
					mc->wfValNew = wf->EvaluateWfSphere(wfData.nbr, mc->u, mc->v);
					mc->MetropolisTestSphere();
					energy=mc->CoulombEnergy(wfData.nbr, mc->u, mc->v);
					break;
				case FQHE::_DISC_:
					mc->RandomMoveDisc();
					mc->wfValNew = wf->EvaluateWfDisc(wfData.nbr, mc->z);
					mc->MetropolisTestDisc();
					energy = mc->CoulombEnergy(wfData.nbr, mc->z);
					break;
				case FQHE::_TORUS_:
					mc->RandomMoveTorus();
					mc->wfValNew = wf->EvaluateWfTorus(wfData.nbr, mc->latticeConfig, wfData.torusState);
					mc->MetropolisTestTorus();
					energy = real(mc->wfValOld);
					break;	
			}
			mc->accumEnergy += energy;
			mc->accumEnergySquared += energy*energy;
			mc->energyCounter++;
			//	Progress monitor:		
			if((k)%((int)ceil((double)mcData.nbrTherm/10))==0)
			{
				mc->CalculateRunningMean();
				utilities::cout.SecondaryOutput() << "\t" << (double)100*(k)/mcData.nbrTherm
				                                  << "%\t\t" << mc->meanEnergy << " +-";
				utilities::cout.SecondaryOutput().precision(3);
				utilities::cout.SecondaryOutput() << mc->stdDevEnergy;
				utilities::cout.SecondaryOutput().precision(6);
				utilities::cout.SecondaryOutput() << "\t\t" << (double)1/3600*(((double)mcData.nbrTherm/(k))-1)
							 *(double) (clock() - timer) / CLOCKS_PER_SEC <<std::endl;
				mc->accumEnergy = 0;
				mc->energyCounter = 0;
				mc->accumEnergySquared = 0;
			}
		}
		utilities::cout.MainOutput().precision(8);
		utilities::cout.MainOutput() << utilities::cout.HyphenLine() << "\n\tDone!" << std::endl;
	}
	utilities::cout.MainOutput() << "\n\tSampling...\n\n" << utilities::cout.HyphenLine() << std::endl;
	utilities::cout.SecondaryOutput() << "\tProgress\tCombined Average\tEstimated Time Remaining (hours)\n" 
	                                  << utilities::cout.HyphenLine() << std::endl;
	double timer;		//	Keep track of the time taken during the simulation
	timer = clock();
	mc->acceptCount = 0;	//	reset the count of accepted moves
	mc->accumEnergy = 0;
	mc->accumEnergySquared = 0;
	mc->energyCounter = 0;
	for(long int k=1; k<=mcData.nbrSteps; ++k)
	{
		double energy = 0;		//	Temporary energy value
		switch (wfData.geometry)
		{
			case FQHE::_SPHERE_:
				mc->RandomMoveSphere();
				mc->wfValNew = wf->EvaluateWfSphere(wfData.nbr, mc->u, mc->v);
				mc->MetropolisTestSphere();
				break;
			case FQHE::_DISC_:
				mc->RandomMoveDisc();
				mc->wfValNew = wf->EvaluateWfDisc(wfData.nbr, mc->z);
				mc->MetropolisTestDisc(); 
				break;
			case FQHE::_TORUS_:
				mc->RandomMoveTorus();
				mc->wfValNew = wf->EvaluateWfTorus(wfData.nbr, mc->latticeConfig, wfData.torusState);
				mc->MetropolisTestTorus();
				break;		
		}
		if(k%(mcData.sampleFreq*wfData.nbr)==0)
		{
			//	store sample configuration
			switch (wfData.geometry)
			{
				case FQHE::_SPHERE_:
					mc->ConfigurationToFileSphere();
					energy = mc->CoulombEnergy(wfData.nbr, mc->u, mc->v);
					break;
				case FQHE::_DISC_:
					mc->ConfigurationToFileDisc();
					energy = mc->CoulombEnergy(wfData.nbr, mc->z);
					break;
				case FQHE::_TORUS_:
					mc->ConfigurationToFileTorus();
					energy = real(mc->wfValOld);
				break;	
			}
			mc->accumEnergy += energy;
			mc->accumEnergySquared += energy*energy;
			mc->energyCounter++;
		}
		//	Progress monitor:		
		if((k)%((int)ceil((double)mcData.nbrSteps/100))==0)
		{
			mc->CalculateRunningMean();
			utilities::cout.SecondaryOutput() << "\t" << (double)100*(k)/mcData.nbrSteps << "%\t\t" 
			                                  << mc->meanEnergy;
			utilities::cout.SecondaryOutput() << "\t\t\t" << (double)1/3600*(((double)mcData.nbrSteps/(k))-1)
			    *(double) (clock() - timer) / CLOCKS_PER_SEC << std::endl;
		}
	}
	utilities::cout.MainOutput() << "\n" << utilities::cout.HyphenLine() << std::endl;
	utilities::cout.MainOutput() << "\tSimulation Completed. \n\n\tTotal time taken for the sampling steps was: ";
	utilities::cout.MainOutput() << (double)1/3600*(double)(clock() - timer) / CLOCKS_PER_SEC
	                             << " Hours.\n" << std::endl;
	utilities::cout.MainOutput().precision(4);
	utilities::cout.MainOutput() << "\tProportion of moves accepted:\t" 
	                             << 100*(double)mc->acceptCount/mcData.nbrSteps << "%" << std::endl << std::endl;
	utilities::cout.MainOutput().precision(8);
	mc->CalculateRunningMean();
	if(wfData.geometry==FQHE::_TORUS_)
	{
		utilities::cout.MainOutput() << "\tMean (real part of) wave function:\t" << mc->meanEnergy;
	}
	else
	{
		utilities::cout.MainOutput() << "\tMean Coulomb energy:\t" << mc->meanEnergy;
	}
	utilities::cout.MainOutput() << "\n\tError estimate (+/-):\t ";
	utilities::cout.MainOutput().precision(4);
	utilities::cout.MainOutput() << mc->stdDevEnergy << std::endl;
	utilities::cout.MainOutput().precision(8);
	if(100*(double)mc->acceptCount/mcData.nbrSteps<49)
	{
		utilities::cout.MainOutput() << "\n\tNOTE: low acceptance ratio (sampling accuracy might be poor)" << std::endl;
		utilities::cout.MainOutput() << "\tConsider decreasing the number of moves per MC step,\n\t";
		utilities::cout.MainOutput() << "or decreasing the maximum move distance.\n" << std::endl;
	}
	if(100*(double)mc->acceptCount/mcData.nbrSteps>85)
	{
		utilities::cout.MainOutput() << "\n\tNOTE: high acceptance ratio (check for sufficient thermalisation)" << std::endl;
		utilities::cout.MainOutput() << "\tConsider increasing the number of moves per MC step,\n\t";
		utilities::cout.MainOutput() << "or increasing the maximum move distance.\n" << std::endl;
	}
	utilities::cout.MainOutput() << utilities::cout.HyphenLine() << std::endl 
	                             << utilities::cout.HyphenLine() << std::endl;
	delete wf;
	delete mc;
	return 0;
}

//!
//! Print a message to display at start of program
//!
void PrintWelcomeMessage()
{
	std::string message=
	"\n\n---------------------------------------------"
	"---------------------------------------------\n"
	"---------------------------------------------"
	"---------------------------------------------\n\n"
	"\t\t   Fractional Quantum Hall Effect Monte Carlo\n\n"
	"\t\t\t   Author: Simon C Davenport	\n\n"
	"\t\t\t\t Version 3.1 \n\n"
	"\t\t\t   Last Modified: 30/11/2012\n\n"
	"\tThe purpose of this program is to generate statistically random\n"
	"\tsample configurations of particle co-ordinates, based on the  \n"
	"\tprobability distribution |psi|^2 associated with a given \n"
	"\tquantum Hall wavefunction.\n\n"
	"\t\t\tCopyright (C) 2012 SIMON C DAVENPORT.\n\n" 
	"\tThis program comes with ABSOLUTELY NO WARRANTY (see readme.txt).  \n"
	"\tThis is free software, and you are welcome to redistribute it \n"
	"\tunder certain conditions (see readme.txt)\n\n"
	"---------------------------------------------"
	"---------------------------------------------\n"
	"---------------------------------------------"
	"---------------------------------------------\n";
	utilities::cout.MainOutput() << message << std::endl;
	return;
}

//!
//! Converts command line arguments into parameters values in global
//!	data structures. 
//!
boost::program_options::variables_map ParseComandLine(
	int argc,	    //!<	Number of characters to parse
	char *argv[])   //!<	Character array to parse
{
	namespace po = boost::program_options;
	po::options_description general("General options");
	general.add_options()
	("help,h", 
	 "Display this message\n")
	("verbose,v",po::value<int>()->default_value(1),
	 "Set a value for the verbosity level:\n\t 0 output off (after command line parsed) \n\t 1 print brief information \n\t 2 print more detialed information \n\t 4 print debugging messages");
	po::options_description all("\tThe program input options are as follows");
	all.add(general).add(FQHE::GetWaveFunctionOptions()).add(FQHE::GetCompositeFermionOptions()).add(FQHE::GetMetropolisOptions());
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv,all), vm);
	po::notify(vm);
	if(vm.count("help")) 
	{
		utilities::cout.MainOutput() << all << "\n";
		exit(EXIT_SUCCESS);
	}
	utilities::cout.SetVerbosity(vm["verbose"].as<int>());
	return vm;
}
