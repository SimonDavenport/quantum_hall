////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 31/01/2015
//!
//!  \file 
//!		This program is designed to calculate an entanglement energy spectrum
//!     due to a single particle entanglement energy model
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

#include "../entanglement_energy_model/entanglement_energy_model.hpp"
#include "../../program_options/general_options.hpp"
#include "../../program_options/rses_options.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include "../../utilities/wrappers/mpi_wrapper.hpp"

///////     GLOBAL DATA STRUCTURES      ////////////////////////////////////////

//  Declare an instance of the global utilities::cout struct for cout_tools
utilities::Cout utilities::cout;

////////////////////////////////////////////////////////////////////////////////

///////      FUNCTION PRE_DECLARATIONS      ////////////////////////////////////
//	See below MAIN for function declarations and description

boost::program_options::variables_map ParseComandLine(int argc,char *argv[]);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	The main function performs three tasks:
//!	 - It parses the command line specified options and puts that information into the
//!	   optionSet variables
//!	 - It generates a basis Hilbert space of  single particle orbital occupation patterns
//!	 - It generates the entanglement energies associated with that basis
//!
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	///////		PARSE THE COMMAND LINE OPTIONS		////////////////////////////

    double timer = clock();

    //  Declare a command line variables map
    boost::program_options::variables_map optionList;

	optionList = ParseComandLine(argc,argv);
	
    //  Diagonalize the specified number of sectors

    int sector       = optionList["sector"].as<int>(); 
	int minSector    = optionList["min-sector"].as<int>();
	int maxSector    = optionList["max-sector"].as<int>();
    bool multiSector = optionList["multi-sector"].as<bool>();

    if(maxSector!=minSector)
	{
	    multiSector = true;
	}

    if(multiSector)
    {
        if(maxSector<minSector)
        {
            std::cerr<<"WARNING: max-sector < min-sector. Automatically swapped!"<<std::endl;
            
            maxSector = minSector;
        }
    }
    else
    {
        maxSector = sector;
        minSector = sector;
    }
    
    //	Print summary of program options

	utilities::cout.MainOutput()<<utilities::cout.EqualsLine()<<std::endl<<utilities::cout.EqualsLine()<<std::endl;
	utilities::cout.MainOutput()<<"\n\tPROGRAM OPTIONS"<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tNUMBER OF PARTICLES\t"<<optionList["nbr-n"].as<int>()<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tN_A PARTICLE IN THE CUT\t"<<optionList["nbr-a"].as<int>()<<std::endl;
	if(!multiSector)    utilities::cout.MainOutput()<<"\n\t\tSECTOR\t\t\t"<<optionList["sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMAX SECTOR\t\t"<<optionList["max-sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMIN SECTOR\t\t"<<optionList["min-sector"].as<int>()<<std::endl;

    //////      Construct basic entanglement energy model
    
    EntanglementEnergyModel model;
    
    model.BuildEntanglementEnergyModel(&optionList,minSector,maxSector);
    
    utilities::cout.MainOutput()<<"\n\tPROGRAM TERMINATED ";
	utilities::cout.MainOutput()<<"\t\tTIME ELAPSED "<<( clock() - timer ) / CLOCKS_PER_SEC<<" SECONDS.\n"<<std::endl;
	utilities::cout.MainOutput()<<utilities::cout.EqualsLine()<<std::endl<<utilities::cout.EqualsLine()<<std::endl;

    return 0;
}

///////      FUNCTION DECLARATIONS      ////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief A function to parse the command line arguments
//!
//! \return An instance of the boost program options variables map
//! containing the parsed command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

boost::program_options::variables_map ParseComandLine(
	int argc,							//!<	Number of characters to parse
	char *argv[])						//!<	Character array to parse
{
	namespace po = boost::program_options;

    //	Combine all the option groups into one
	po::options_description allOpt("\n\n\tThis program is designed to calculate an entanglement energy spectrum due to a single particle entanglement energy model. \n\n\tThe program input options are as follows");
    //allOpt.add_options()
    //("use-perturbation",po::value<bool>()->default_value(false),
	//"Set to 1 in order to include 1st order perturbation correction\n");
    
	allOpt.add(myOptions::GetGeneralOptions()).add(myOptions::GetEntanglementSpectrumOptions());

	po::variables_map vm;
    
    try
    {
        //	Map the command line argument onto the option set
        po::store(po::parse_command_line(argc, argv,allOpt), vm);
        
        //	Respond to help option declaration
        if (vm.count("help"))
        {
            utilities::cout.MainOutput() << allOpt << "\n";
            exit(EXIT_SUCCESS);
        }
        
        po::notify(vm);
    }
    catch(po::error& e)
    {
        utilities::cout.MainOutput()<<allOpt<<std::endl;
        
        std::cerr<<utilities::cout.HyphenLine()<<std::endl;
        
        std::cerr<<std::endl<<"\tERROR:\t"<<e.what()<<std::endl;
        
        std::cerr<<std::endl<<utilities::cout.HyphenLine()<<std::endl;
        
        exit(EXIT_FAILURE);
    }
			
	//  Set global verbosity level
	
	utilities::cout.SetVerbosity(vm["verbose"].as<int>());

	return vm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
