//////////////////////////////////////////////////////////////////////////////////
//!
//!							Author: Simon Davenport
//!
//!		The purpose of this program is to analyse data produced by the Monte
//!		Carlo program fqheMonteCarlo.cpp
//!
//! 				Copyright (C) Simon C Davenport
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
#include "../analysis/fqhe_monte_carlo_analysis.hpp"
#if _ENABLE_MPI_
#include "../../utilities/wrappers/mpi_wrapper.hpp"
#endif
utilities::Cout utilities::cout;
///////		FUNCTION FORWARD DELCATATIONS		    ////////////////////////////
void PrintWelcomeMessage();
boost::program_options::variables_map ParseComandLine(int argc,char *argv[]);
#if _ENABLE_MPI_
utilities::MpiWrapper mpi(utilities::cout);
#endif

int main(int argc, char *argv[])
{
    #if _ENABLE_MPI_ 
	mpi.Init(argc,argv);
    #endif
	utilities::cout.MainOutput().precision(8);
	PrintWelcomeMessage();
	boost::program_options::variables_map optionList;
	optionList = ParseComandLine(argc, argv);
    AnalysisMethods analyse;
    #if _ENABLE_MPI_
    analyse.InitFromCommandLine(&optionList, mpi);
    #else
    analyse.InitFromCommandLine(&optionList);
    #endif
	if(analyse.options.getFileName)
	{
		analyse.PrintName();
		return 0;
	}
	if(analyse.LookForConfigurationDataFiles()==EXIT_FAILURE) 
	{
	    return EXIT_FAILURE;
	}
	std::cout << std::endl << utilities::cout.HyphenLine() << std::endl 
	          << utilities::cout.HyphenLine() << std::endl;
	std::cout << "\n\tAnalyzing data..." << std::endl;
	analyse.InitializeMethods();
	std::cout << std::endl << utilities::cout.HyphenLine() << std::endl#
	          << utilities::cout.HyphenLine() << std::endl;
	int maxk;
	size_t size=0;
	double timer=clock();	
	int sample = 0;
	int runNo = 1;
	std::cout.precision(15);
	while(runNo<=analyse.maxRuns)
	{
		std::cout << "\n\tImporting data set " << runNo << " of " << analyse.maxRuns << std::endl;
		analyse.OpenConfigurationDataFile(analyse.options.dataDirIn, runNo);
		std::cout << std::endl;
		if(analyse.wfData.geometry == FQHE::_TORUS_)
		{
			analyse.f_config.seekg(0, std::ios::end);
			size = analyse.f_config.tellg();
			maxk = size/(analyse.wfData.nbr*sizeof(std::complex<int>))-analyse.options.skip;
			analyse.f_config.seekg (0, std::ios::beg);
		}
		else
		{
			analyse.f_phi.seekg (0, std::ios::end);
			size = analyse.f_phi.tellg();
			maxk = size/(analyse.wfData.nbr*sizeof(double))-analyse.options.skip;
			analyse.f_phi.seekg (0, std::ios::beg);
		}
		//	TODO set the file stream pointer rather than needlessly reading values
		for(int k=0; k<analyse.options.skip; ++k)
		{
			analyse.GetConfigFromFileSphere();
		}
		//	if skip>0, skip samples
		std::cout << "\tContains " << maxk+analyse.options.skip 
		          << " sample configurations of which " << analyse.options.skip
		          << " are skipped" << std::endl << std::endl;
		if(analyse.wfData.geometry==FQHE::_TORUS_ || analyse.wfData.geometry==FQHE::_DISC_)
		{
			std::cerr << "\tFATAL ERROR: torus and disc goemetry not completed yet!" << std::endl;
			return EXIT_FAILURE;
			//	TODO fix this
		}
		for(int k=0; k<maxk; ++k)
		{
			analyse.GetConfigFromFileSphere();
			analyse.PolarToSpinor();
			if(analyse.options.check)	
			{
			    analyse.CheckConfig(k);
			}
			if(analyse.options.pairCorrel)
			{
			    analyse.PairCorrelation(analyse.wfData.nbrUp, analyse.wfData.nbrDown, 
			                            analyse.u, analyse.v);
			}
			if(analyse.options.dens)			
			{
			    analyse.Density(analyse.wfData.nbr, analyse.u);
		    }
			if(analyse.options.coulomb)			
			{
			    analyse.coulombEnergy[sample] = analyse.CoulombEnergy(analyse.wfData.nbr, 
			                                                          analyse.u, analyse.v);
			if(analyse.options.secllCoulomb)	
			{
				analyse.secondEnergy[sample] = analyse.SecondLLEnergy(analyse.wfData.nbr, 
				                                                      analyse.u, analyse.v);
			}
			if(analyse.options.finiteThickness) 
			{
				double thicknessValue = analyse.options.minThickness;
				int counter = 0;
				while(thicknessValue <= analyse.options.maxThickness)
				{
					analyse.thicknessEnergy[counter] += analyse.FiniteThicknessEnergy(
					    analyse.wfData.nbr, thicknessValue, analyse.u, analyse.v);
					thicknessValue += analyse.options.thicknessStepSize;
					++counter;
				}
			}
			if(analyse.options.bilayer)
			{
				double bilayerValue = analyse.options.minBilayer;
				int counter = 0;
				while(bilayerValue <= analyse.options.maxBilayer)
				{
					analyse.bilayerEnergy[counter] += analyse.BilayerEnergy(
					    analyse.wfData.nbrUp, analyse.wfData.nbrDown, bilayerValue, analyse.u, analyse.v);
					bilayerValue += analyse.options.bilayerStepSize;
					++counter;
				}
			}
			if(analyse.options.chargePlate)
			{
				double chargePlateValue = analyse.options.minPlate;
				int counter = 0;
				while(chargePlateValue <= analyse.options.maxPlate)
				{
					analyse.chargePlateEnergy[counter] += analyse.ChargePlateEnergy(
					    analyse.wfData.nbr, chargePlateValue, analyse.u, analyse.v);
					chargePlateValue += analyse.options.plateStepSize;
					++counter;
				}
			}
			if(analyse.options.biChargePlate)
			{
				int counter = 0;
				for(int d=0; d<25; ++d)
				{
					for(int D=0; D<20; ++D)
					{
						analyse.biChargePlateEnergy[counter] += analyse.BiChargePlateEnergy(
						    analyse.wfData.nbrUp, analyse.wfData.nbrDown, D, d, counter, analyse.u, analyse.v);
						++counter;	
					}
				}
			}
			if(analyse.options.partialCoulomb)
			{
				analyse.partialEnergy[sample] = analyse.PartialCoulombEnergy(analyse.wfData.nbrUp, 
				    analyse.wfData.nbrDown, analyse.options.partialCoulombType, analyse.u, analyse.v);
			}
			if((k+1)%(maxk/10)==0)
			{
				std::cout.precision(3);
				std::cout << "\t" << (double)(100*(k+1))/maxk << "% completed\t" << "Time remaining (hours):\t";
				std::cout.precision(4);
				std::cout << (((double)maxk/(k+1))-1)/3600*(double)(clock() - timer) / CLOCKS_PER_SEC << std::endl;
				fflush(stdout);
				std::cout.precision(12);
			}
			++sample;
		}
		++runNo;
	}
//////      Resample energy data in order to determine proper standard      //////
//////      deviation for correlated MC sampling.	                        /////
	if(analyse.options.resample)
	{
		std::cout << "\n\tDone!" << std::endl;
		std::cout << std::endl << utilities::cout.HyphenLine() << std::endl 
		          << utilities::cout.HyphenLine() << std::endl;	
		if(analyse.options.coulomb)			
		{
			std::cout << std::endl << "\tRe-sampling coulomb energy...\n" << std::endl;
			analyse.Resample(analyse.coulombEnergy,analyse.options.nbrResamples,analyse.coulombStdev,analyse.coulombMeanEnergy);
		}
		if(analyse.options.secllCoulomb)
		{
			std::cout<<std::endl<<"\tRe-sampling 2nd LL coulomb energy...\n"<<std::endl;
			analyse.Resample(analyse.secondEnergy, analyse.options.nbrResamples, 
			                 analyse.secondStdev, analyse.secondMeanEnergy);
		}
		if(analyse.options.partialCoulomb)
		{
			std::cout << std::endl << "\tRe-sampling partial coulomb potential energy...\n" << std::endl;
			analyse.Resample(analyse.partialEnergy, analyse.options.nbrResamples, 
			                 analyse.partialStdev, analyse.partialMeanEnergy);
		}
	}
	else
	{
		std::cout << "\n\tDone!\n\n\tSkipping resampling.\n" << std::endl;
		if(analyse.options.coulomb)			
		{
			analyse.MeanEnergy(analyse.coulombEnergy, analyse.coulombStdev, analyse.coulombMeanEnergy);
		}
		if(analyse.options.secllCoulomb)
		{
			analyse.MeanEnergy(analyse.secondEnergy, analyse.secondStdev, analyse.secondMeanEnergy);
		}
		if(analyse.options.partialCoulomb)
		{
			analyse.MeanEnergy(analyse.partialEnergy,analyse.partialStdev,analyse.partialMeanEnergy);
		}
	}
	std::cout << std::endl << utilities::cout.HyphenLine() << std::endl << utilities::cout.HyphenLine() 
	          << std::endl << std::endl;
	analyse.FinaliseMethods();
	std::cout << utilities::cout.HyphenLine() << std::endl << utilities::cout.HyphenLine() << std::endl;	
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
    "\t\t\t     analyse_monte_carlo.cpp\n\n"
    "\t\t\t   Author: Simon C Davenport	\n\n"
    "\t\t\t\t Version 3.0 \n\n"
    "\t\t\t   Last Modified: 10/10/2012\n\n"
    "\tThe purpose of this program is to analyse data produced by the Monte\n"
    "\tCarlo program fqheMonteCarlo.cpp.\n\n"
    "\t\t\tCopyright (C) 2012 SIMON C DAVENPORT.\n\n" 
    "\tThis program comes with ABSOLUTELY NO WARRANTY (see ReadMe.txt).  \n"
    "\tThis is free software, and you are welcome to redistribute it \n"
    "\tunder certain conditions (see ReadMe.txt)\n\n"
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
	all.add(general).add(FQHE::GetWaveFunctionOptions()).add(FQHE::GetCompositeFermionOptions()).add(GetAnalysisOptions());
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, all), vm);
	po::notify(vm);
	if(vm.count("help"))
	{
		utilities::cout.MainOutput() << all << "\n";
		exit(EXIT_SUCCESS);
	}
	utilities::cout.SetVerbosity(vm["verbose"].as<int>());
	return vm;
}
