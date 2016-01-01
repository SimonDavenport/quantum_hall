////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 31/01/2015
//!
//!  \file
//!		This program calculates the exact RSES for integer quantum Hall 
//!     states using a numerically exact entanglement energy method.
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

//  Functions to calculate the RSES
#include "../iqhe_rses/iqhe_rses.hpp"

//  Functions to import program options
#include "../../program_options/general_options.hpp"
#include "../../program_options/rses_options.hpp"

//  Functions to manipulate std::cout output
#include "../../utilities/general/cout_tools.hpp"

//  For clock()
#include <time.h>                       

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

utilities::Cout utilities::cout;

///////      FUNCTION PRE_DECLARATIONS      ////////////////////////////////////
//	See below MAIN for function declarations and description

boost::program_options::variables_map ParseComandLine(int argc,char *argv[]);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	The main function imports command line arguments, then calls functions
//! to calculate the RSES spectrum for the selected IQHE state
//!
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	double timer = clock();
	
	///////		PARSE THE COMMAND LINE OPTIONS		////////////////////////////

    //  Declare a command line variables map
    boost::program_options::variables_map optionList;

	optionList = ParseComandLine(argc,argv);

    //  Diagonalize the specified number of sectors

    int sector       = optionList["sector"].as<int>(); 
	int minSector    = optionList["min-sector"].as<int>();
	int maxSector    = optionList["max-sector"].as<int>();
    bool multiSector = false;

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
	utilities::cout.MainOutput()<<"\n\t\tNBR LLs\t\t"<<optionList["nbr-ll"].as<int>()<<std::endl;
	if(!multiSector)    utilities::cout.MainOutput()<<"\n\t\tSECTOR\t\t\t"<<optionList["sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMAX SECTOR\t\t"<<optionList["max-sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMIN SECTOR\t\t"<<optionList["min-sector"].as<int>()<<std::endl;

    ///////     GENERATE RSES       ////////////////////////////////////////////

    IqheRses iqhe;
    
    iqhe.GenerateIqheSpectrum(&optionList,minSector,maxSector);
    
    //////      END PROGRAM     ////////////////////////////////////////////////
    
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

	po::options_description allOpt("\n\tThis program calculates the exact RSES for integer quantum Hall states using a numerically exact entanglement energy method.\n\n\tThe program input options are as follows");
	
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
