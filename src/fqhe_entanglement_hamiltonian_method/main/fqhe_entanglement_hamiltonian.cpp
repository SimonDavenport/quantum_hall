////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 31/01/2015
//!
//!  \file
//!		This program is designed to calculate and diagonalize a generic
//!		'entanglement Hamiltonian' written in terms of Virosoro algebra operators.
//!		See Phys. Rev. B 86, 245310 (2012).
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

//	Colour labels _U1_LABEL_1_, _U1_LABEL_2_, _U1_LABEL_3_ and _MAJORANA_CURRENT_
//  are defined in term.h

#include "../entanglement_hamiltonian/term.hpp"
#include "../entanglement_hamiltonian/list_of_terms.hpp"
#include "../entanglement_hamiltonian/hilbert.hpp"
#include "../../program_options/general_options.hpp"
#include "../../program_options/rses_options.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include <time.h>   //  For clock()

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

utilities::Cout utilities::cout;

////////////////////////////////////////////////////////////////////////////////
#if _BENCHMARK_MODE_==1	

int g_callsCtor=0;
int g_callsCommuteRight=0;
int g_callsSameTest=0;
int g_callsZeroTest=0;
int g_callsCommutator=0;

void PrintCalls()
{
	utilities::cout.DebuggingInfo()<<"\n\tCALLS TO TERM CTOR: "<<g_callsCtor<<std::endl;
	utilities::cout.DebuggingInfo()<<"\n\tCALLS TO TERM COMMUTE RIGHT: "<<g_callsCommuteRight<<std::endl;
	utilities::cout.DebuggingInfo()<<"\n\tCALLS TO TERM SAME TEST: "<<g_callsSameTest<<std::endl;
	utilities::cout.DebuggingInfo()<<"\n\tCALLS TO TERM ZERO TEST: "<<g_callsZeroTest<<std::endl;
	utilities::cout.DebuggingInfo()<<"\n\tCALLS TO COMMUTATOR: "<<g_callsCommutator<<std::endl;
}

#endif
////////////////////////////////////////////////////////////////////////////////

///////      FUNCTION PRE_DECLARATIONS      ////////////////////////////////////
//	See below MAIN for function declarations and description

boost::program_options::variables_map ParseComandLine(int argc,char *argv[]);

ListOfTerms TestTerm();

ListOfTerms DubailReadRezai431(boost::program_options::variables_map* optionList,
                               const int lzA);

ListOfTerms CompositeFermionModel(boost::program_options::variables_map* optionList,
                                  const int lzA);

ListOfTerms DubailReadRezai431(const double alpha,const double beta,
				  const double gamma,const double fillingFactor,
				  const int lzA,const char colour,const int nbr);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief	The main function performs three tasks:
//!	 - It parses the command line specified options and puts that information into the
//!	   optionSet variables
//!	 - It generates a basis Hilbert space in which the entanglement Hamiltonian 
//!	   can be expressed as a matrix
//!	 - It generates the entanglement Hamiltonian selected from the program options
//!	   and calls functions to diagonalize it. 
//!
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	double timer = clock();
	
	///////		PARSE THE COMMAND LINE OPTIONS		////////////////////////////

    //  Declare a command line variables map
    boost::program_options::variables_map optionList;

	optionList = ParseComandLine(argc,argv);

    if(abs(optionList["statistics"].as<int>()) != 1)
    {
        std::cerr<<"\n\tERROR: Invalid statistics option: "<<optionList["statistics"].as<int>()<<". Choose +1 for fermions or -1 for bosons only."<<std::endl;
        exit(EXIT_FAILURE);    
    }

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
	utilities::cout.MainOutput()<<"\n\t\tOPTIMIZING FOR BLOCKS?\t"<<(optionList["blocks"].as<bool>() ? "true" : "false")<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tNUMBER OF PARTICLES\t"<<optionList["nbr-n"].as<int>()<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tN_A PARTICLE IN THE CUT\t"<<optionList["nbr-a"].as<int>()<<std::endl;
	if(!multiSector)    utilities::cout.MainOutput()<<"\n\t\tSECTOR\t\t\t"<<optionList["sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMAX SECTOR\t\t"<<optionList["max-sector"].as<int>()<<std::endl;
	if(multiSector)     utilities::cout.MainOutput()<<"\n\t\tMIN SECTOR\t\t"<<optionList["min-sector"].as<int>()<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tNBR U(1) CURRENTS\t"<<optionList["nbr-ll"].as<int>()<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tINCLUDE MAJORANA MODES\t"<<(optionList["majorana-modes"].as<bool>() ? "YES" : "NO")<<std::endl;
	utilities::cout.MainOutput()<<"\n\t\tSTATISTICS\t\t"<<(optionList["statistics"].as<int>() == 1 ? "bosons" :"fermions"  )<<std::endl;
    utilities::cout.MainOutput()<<"\n\t\tPARAMETERS FILE\t\t"<<(optionList["parameters-file"].as<std::string>())<<std::endl;

    const std::string hamiltonian = optionList["hamiltonian"].as<std::string>();

    ///////     GENERATE ENTANGLEMENT HAMILTONIAN     //////////////////////////////

    for(sector = minSector;sector<=maxSector;sector++)
    {
        utilities::cout.MainOutput()<<"\n\t============================    SECTOR "<<sector<<"    ============================\n"<<std::endl;
         
        //  Build Hilbert space

        HilbertSpace basis;
    
        basis.BuildHilbertSpace(&optionList,sector);

        //  Generate Hamiltonian

        if("DRR431" == hamiltonian)
        {
            basis.GenerateMatrixElements(DubailReadRezai431(&optionList,sector));
        }
        else if("extended_DRR" == hamiltonian)
        {
            basis.GenerateMatrixElements(CompositeFermionModel(&optionList,sector));
        }
        else
        {
            std::cerr<<"\n\tERROR: unknown type specified: "<<hamiltonian<<std::endl;
        }
 
        //  Diagonalize the Hamiltonian
        
        basis.Diagonalize();
        
        basis.EigenvaluesToFile(&optionList,sector);
    }

	utilities::cout.MainOutput()<<"\n\tPROGRAM TERMINATED ";
	utilities::cout.MainOutput()<<"\t\tTIME ELAPSED "<<( clock() - timer ) / CLOCKS_PER_SEC<<" SECONDS.\n"<<std::endl;
	utilities::cout.MainOutput()<<utilities::cout.EqualsLine()<<std::endl<<utilities::cout.EqualsLine()<<std::endl;
	
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1	
	
	PrintCalls();
	
	#endif
	////////////////////////////////////////////////////////////////////////////////
	
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

	//	Entanglement Hamiltonian options
	po::options_description hamiltonianOpt("Hamiltonian Options");
	hamiltonianOpt.add_options()
	("hamiltonian",po::value<std::string>()->default_value("DRR431"),
	"Specify the entanglement Hamiltonian to diagonalize: \n\t DRR431 - Dubail-Read-Rezai Eq. 4.31\n\t extended_DRR - Composite Fermion model.\n")
    ("parameters-file",po::value<std::string>()->default_value("parameters.dat"),
	 "Specify the name of the file containing the model parameters\n")
	("power,p",po::value<int>()->default_value(3),
	"Exponent of Jastrow factor i.e. (z_i - z_j)^p (used only in calculation of filling factor).\n")
	("blocks,b", po::value<bool>()->default_value(0),
	"Set to 1 to optimize for the block structure of the Hamiltonian, if available.\n")
	("statistics", po::value<int>()->default_value(-1),
	"Set to -1 for fermion statistics or +1 for boson statistics (used only in calculation of filling factor).\n")
	("majorana-modes", po::value<bool>()->default_value(false),
	"Set to 1 in order to include the majorana modes in the entanglement Hilbert space, in addition to U(1) modes.\n");
	
	//	Combine all the option groups into one
	po::options_description allOpt("\n\tThis program is designed to calculate and diagonalize a generic 'entanglement Hamiltonian' written in terms of Virosoro algebra operators. See Phys. Rev. B 86, 245310 (2012). Note that parameters in the Hamiltonian are set from the parameters.dat file. \n\n\tThe program input options are as follows");
    
    allOpt.add(myOptions::GetGeneralOptions()).add(myOptions::GetEntanglementSpectrumOptions()).add(hamiltonianOpt);

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

	return vm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Allows us to define and calculate the matrix element of a single term
//! 
//!	\return  A vector array containing a representation of the Hamiltonian
//! in terms of Term objects
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms TestTerm()
{
	ListOfTerms op;
	std::vector<short int> indices;
	std::vector<char> labels;
	double coefficient;

	//	currently set for J_2 J_2 J_-2 J_-2
	
	indices.push_back(2);
	indices.push_back(3);
	indices.push_back(-2);
	indices.push_back(-2);

	for(unsigned int i=0;i<indices.size();i++)
	{
		labels.push_back(_U1_LABEL_1_);
	}

	coefficient=1.0;

	op += Term(coefficient,indices,labels);

    op.PrintOperator();

	return op;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Returns a representation of the entanglement Hamiltonian given in
//!	Eq 4.31 of PRB 86,245310. (Reads parameters from a file).
//!
//!	Number of terms to include in Fourier expansion of terms is the same as
//!	the sector since terms containing J_n>m will commute with bras/kets
//! containing only up to J_m.
//!
//!	\return     A vector array containing a representation of the Hamiltonian
//! in terms of Term objects
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms DubailReadRezai431(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzA)                                      //!<    Angular momentum
{
	double timer = clock();
	
	//	Parameters in DRR 4.31 Hamiltonian

	double alpha;
	double beta; 
	double gamma;
	double fillingFactor;

	fillingFactor = (double)1.0/((*optionList)["power"].as<int>());
	
	std::stringstream fileName;		//	string stream to build file name
	
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==0	
	
	//	read values from a file

	std::ifstream f_param;			//	file stream to read in parameter file

	fileName.str("");
	
	//fileName<<(*optionList)["in-path"].as<std::string>()<<"parameters_DRR431_Model_n_"
	//        <<(*optionList)["nbr-n"].as<int>()<<"_nA_"<<(*optionList)["nbr-a"].as<int>()<<".dat";
	
    fileName<<(*optionList)["in-path"].as<std::string>()<<(*optionList)["parameters-file"].as<std::string>();
    
	f_param.open(fileName.str().c_str());
	
	if(f_param==NULL)
	{
		std::cerr<<"\tFile: "<<fileName.str()<<" not found !!! "<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}

	//	skip the comment lines in the file

	std::string line;

	while(getline(f_param,line))
	{
		if(line[0]!='#')	break;
	}

	f_param>>alpha;
	f_param>>beta;
	f_param>>gamma;

	f_param.close();
	
	#endif
    ////////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1

			alpha  = 1.20935651933;
			beta   = -0.240908896371;
			gamma  = 1.11000998789;

	#endif
	////////////////////////////////////////////////////////////////////////////////
		
	//  Call function to generate the Hamiltonian
	
	ListOfTerms hamiltonian(DubailReadRezai431(alpha,beta,gamma,fillingFactor,lzA,_U1_LABEL_1_,(*optionList)["nbr-n"].as<int>()));

	//	Print out the final result:

	utilities::cout.MainOutput()<<"\n\t- USING READ-DUBAIL-REZAYI EQ 4.31 with:\n"<<std::endl;
	utilities::cout.MainOutput()<<"\t\talpha = "<<alpha<<std::endl;
	utilities::cout.MainOutput()<<"\t\tbeta = "<<beta<<std::endl;
	utilities::cout.MainOutput()<<"\t\tgamma = "<<gamma<<std::endl;
	
	utilities::cout.MainOutput()<<"\n\t\t(READING PARAMETER DATA FROM FILE: "<<fileName.str()<<" )."<<std::endl;

    utilities::cout.SecondaryOutput()<<"\n\t- U(1) OPERATOR REPRESENTATION OF THE HAMILTONIAN:\n\n\t";
    
    hamiltonian.PrintOperator();

	utilities::cout.MainOutput()<<"\n\t\tTIME TAKEN TO BUILD HAMILTONIAN "<<( clock() - timer ) / CLOCKS_PER_SEC<<" SECONDS.\n"<<std::endl;
	
	return hamiltonian;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief	Returns a test representation of the entanglement Hamiltonian used
//!	to model the Composite Fermion entanglement spectrum. 
//!
//! In this model the entanglement Hamiltonian includes both the set of two/three
//! identical copies of the DubailReadRezai431() entanglement Hamiltonian,
//! one copy for each effective CF LL in the model, and an additional
//! composite fermion "kinetic energy" term.
//!
//!	Number of terms to include in Fourier expansion of terms is the same as
//!	the sector since terms containing J_n>m will commute with bras/kets
//! containing only up to J_m.
//!
//!	\return     A vector array containing a representation of the Hamiltonian
//! in terms of Term objects
//! 
////////////////////////////////////////////////////////////////////////////////

ListOfTerms CompositeFermionModel(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzA)                                      //!<    Angular momentum
{
	double timer = clock();

	//	In general we can specify one set of DRR 4.31 Hamiltonian parameters
	//  for each copy of the DRR 4.31 wave function in this construction

    const int nbrCurrents = (*optionList)["nbr-ll"].as<int>();
	double alpha[nbrCurrents];
	double beta[nbrCurrents]; 
	double gamma[nbrCurrents];
	double cfCyclotronEnergy;           //  An effective scale for hbar omega_c
	double cfCyclotronEnergy2;          //  A quadratic correction to the cyclotron energy
	double fillingFactor;

	fillingFactor = (double)nbrCurrents/(2.0*nbrCurrents-(*optionList)["statistics"].as<int>());
	
	std::stringstream fileName;		//	string stream to build file name
	
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==0	
	
	//	read values from a file

	std::ifstream f_param;			//	file stream to read in parameter file
	
	fileName.str("");
	
	//fileName<<(*optionList)["in-path"].as<std::string>()<<"parameters_CF_Model_n_"<<(*optionList)["nbr-n"].as<int>()
	//        <<"_nA_"<<(*optionList)["nbr-a"].as<int>()<<"_LL_"<<nbrCurrents<<".dat";
	
    fileName<<(*optionList)["in-path"].as<std::string>()<<(*optionList)["parameters-file"].as<std::string>();
    
	f_param.open(fileName.str().c_str());
	
	if(f_param==NULL)
	{
		std::cerr<<"\tFile: "<<fileName.str()<<" not found !!! "<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}

	//	skip the comment lines in the file

	std::string line;

	while(getline(f_param,line))
	{
		if(line[0]!='#')	break;
	}

    for(int i=0;i<nbrCurrents;i++)
    {
	    f_param>>alpha[i];
	    f_param>>beta[i];
	    f_param>>gamma[i];
	}

    f_param>>cfCyclotronEnergy;
    f_param>>cfCyclotronEnergy2;

	f_param.close();

	#endif
    ////////////////////////////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////////////////////////////
	#if _BENCHMARK_MODE_==1

		alpha[0] = 1;
		beta[0]  = 1;
		gamma[0] = 1;
		
		alpha[1] = 1;
		beta[1]  = 1;
		gamma[1] = 1;
		
		if(nbrCurrents>2)
		{
			alpha[2] = 1;
			beta[2]  = 1;
			gamma[2] = 1;
		}
		
		cfCyclotronEnergy=1;
		
	#endif
	////////////////////////////////////////////////////////////////////////////////
	
	//  Call function(s) to generate the Hamiltonian
	
	//  In this case we want to construct an independent copy of the DRR_431 Hamiltonian 
	//  with parameters corresponding to the occupations of each CF LL
	//  Mostly we just need to get the correct cut-off sector and put the right particle
	//  cuts in the CF LLs. So start with, these are set to be identical to each other. 
	//	The Dubail_Read_Rezayi function is called as:
	//	std::vector<Term> DubailReadRezai431(const double alpha,
	//	const double beta,const double gamma,const int nbrA,const int lzA,const char colour)

	ListOfTerms hamiltonian(DubailReadRezai431(alpha[0],beta[0],gamma[0],fillingFactor,lzA,_U1_LABEL_1_,(*optionList)["nbr-n"].as<int>()));
	
	if(nbrCurrents>1)
	{
	    hamiltonian += DubailReadRezai431(alpha[1],beta[1],gamma[1],fillingFactor,lzA,_U1_LABEL_2_,(*optionList)["nbr-n"].as<int>());
    }
    
    if(nbrCurrents>2)
	{
	    hamiltonian += DubailReadRezai431(alpha[2],beta[2],gamma[2],fillingFactor,lzA,_U1_LABEL_3_,(*optionList)["nbr-n"].as<int>());
    }
    
    if(nbrCurrents>3)
    {
        std::cerr<<"\n\tERROR in construction of extended_DRR entanglement Hamiltonian: You asked it to build a Hamiltonian for >3 CF LLS. "<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    //  Explicitly add in a CF "kinetic energy" term i.e. Energy*J_0 + (Energy+1)*K_0
    //  in units of hbar omega_c
    
    if(nbrCurrents>1)
	{
        std::vector<short int> indices(1);
	    std::vector<char> labels(1);
	    double coefficient;

	    indices[0]=0;

	    labels[0]=_U1_LABEL_1_;
	
	    coefficient = (0.5)*cfCyclotronEnergy;///sqrt(fillingFactor);

        hamiltonian += Term(coefficient,indices,labels);
        
        indices[0]=0;

	    labels[0]=_U1_LABEL_2_;
	
	    coefficient = (1.0+0.5)*cfCyclotronEnergy;///sqrt(fillingFactor);

        hamiltonian += Term(coefficient,indices,labels);
    
        if(nbrCurrents>2)
	    {
	        indices[0]=0;

	        labels[0]=_U1_LABEL_3_;
	
	        coefficient = (2.0+0.5)*cfCyclotronEnergy;///sqrt(fillingFactor);

            hamiltonian += Term(coefficient,indices,labels);
	    } 
    }
    
    //  Include a term that is quadratic in J_0
    
    if(nbrCurrents>1)
	{
        std::vector<short int> indices(2);
	    std::vector<char> labels(2);
	    double coefficient;

	    indices[0]=0;
	    indices[1]=0;

	    labels[0]=_U1_LABEL_1_;
	    labels[1]=_U1_LABEL_1_;
	
	    coefficient = (0.5)*cfCyclotronEnergy2;///sqrt(fillingFactor);

        hamiltonian += Term(coefficient,indices,labels);
        
        indices[0]=0;
        indices[1]=0;

	    labels[0]=_U1_LABEL_2_;
	    labels[1]=_U1_LABEL_2_;
	
	    coefficient = (1.0+0.5)*cfCyclotronEnergy2;///sqrt(fillingFactor);

        hamiltonian += Term(coefficient,indices,labels);
    
        if(nbrCurrents>2)
	    {
	        indices[0]=0;
	        indices[1]=0;

	        labels[0]=_U1_LABEL_3_;
	        labels[1]=_U1_LABEL_3_;
	
	        coefficient = (2.0+0.5)*cfCyclotronEnergy2;///sqrt(fillingFactor);

            hamiltonian += Term(coefficient,indices,labels);
	    }
    
    }
    
    #if 0
    //	Add in the possible terms that mix number operators
    
    {
    
    std::vector<short int> indices;
	std::vector<char> labels;
	double coefficient;

	indices.push_back(0);
	indices.push_back(0);

	labels.push_back(_U1_LABEL_1_);
	labels.push_back(_U1_LABEL_2_);
	
	coefficient = 1.0;

    hamiltonian.push_back(Term(coefficient,indices,labels));
    
    }
    #endif

	//	Print out the final result:

	utilities::cout.MainOutput()<<"\n\t- USING COMPOSITE FERMION MODEL WITH:\n"<<std::endl;

	for(int i=0;i<nbrCurrents;i++)
    {
	    utilities::cout.MainOutput()<<"\t\talpha["<<i<<"] = "<<alpha[i]<<std::endl;
	    utilities::cout.MainOutput()<<"\t\tbeta["<<i<<"] = "<<beta[i]<<std::endl;
	    utilities::cout.MainOutput()<<"\t\tgamma["<<i<<"] = "<<gamma[i]<<std::endl;
	}
	
	utilities::cout.MainOutput()<<"\t\tCyclotron energy = "<<cfCyclotronEnergy<<std::endl;
	
	//utilities::cout.MainOutput()<<"\n\t\t(READING PARAMETER DATA FROM FILE: "<<fileName.str()<<" )."<<std::endl;

    utilities::cout.SecondaryOutput()<<"\n\t- U(1) OPERATOR REPRESENTATION OF THE HAMILTONIAN:\n\n\t";

    hamiltonian.PrintOperator();

	utilities::cout.MainOutput()<<"\n\t\tTIME TAKEN TO BUILD HAMILTONIAN "<<( clock() - timer ) / CLOCKS_PER_SEC<<" SECONDS.\n"<<std::endl;
	
	return hamiltonian;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief  This is function overload for the DubailReadRezai431 function which 
//! Can be used to produce a Hamiltonian of the DRR 4.31 form for a specified
//! parameter set and a given colour.
//!
//!	Returns a representation of the entanglement Hamiltonian given in Eq 4.31
//!	of PRB 86,245310.
//!
//!	Number of terms to include in Fourier expansion of terms is the same as
//!	the sector since terms containing J_n>m will commute with bras/kets 
//! containing only up to J_m.
//!
//!	\return     A vector array containing a representation of the Hamiltonian
//! in terms of Term objects
//!
////////////////////////////////////////////////////////////////////////////////

ListOfTerms DubailReadRezai431(
    const double alpha,  //!<    Coefficient of the energy-tensor term
    const double beta,   //!<    Coefficient of first order dispersion term
    const double gamma,  //!<    Coefficient of second order dispersion
    const double fillingFactor,   //!<    Filling factor to normalize the number operators
    const int lzA,       //!<    The angular momentum sector (used for selecting an appropriate cut-off)
    const char colour,   //!<    The colour label of the current operator
    const int nbr)       //!<    Number of particles in the state
{	
    //  Set the truncation index in Fourier summations

	const int fourierTruncation = lzA;
	
	//  Set the N_A value of the lowest pseudo-energy cut
	
	int nbrA0 = (int)floor(nbr/2);

	std::vector<Term> hamiltonian;

	//	Generate the L_0 operator:
	//
	//	alpha/sqrt(N_A) * [ 0.5 J_0^2  + sum_k>=1  J_-k J_k ]

	for(int k=0;k<=fourierTruncation;k++)
	{
		std::vector<short int> indices;
		std::vector<char> labels;
		double coefficient;

		indices.push_back(-k);
		indices.push_back(k);

		labels.push_back(colour);
		labels.push_back(colour);

		if(k==0)
		{		
			coefficient = 0.5/fillingFactor;
		}
		else
		{	
			coefficient = 1.0;
		}

		coefficient *= alpha/sqrt(nbrA0);
		
		hamiltonian.push_back(Term(coefficient,indices,labels));
	}

	//	Generate O_0^1 operator:
	//
	//	beta/N_A^(3/2) * [ 0.5 J_0^2 + sum_k>=2 (1 - k^2) J_-k J_k ]

	for(int k=0;k<=fourierTruncation;k++)
	{
		if(k!=1)
		{
			std::vector<short int> indices;
			std::vector<char> labels;
			double coefficient;
	
			indices.push_back(-k);
			indices.push_back(k);
	
			labels.push_back(colour);
			labels.push_back(colour);
	
			if(k==0)
			{
				coefficient = 0.5/fillingFactor;
			}
			else
			{
				coefficient = (1.0-k*k);
			}
			
			coefficient *= beta/pow(nbrA0,1.5);
	
			hamiltonian.push_back(Term(coefficient,indices,labels));
		}
	}

	//	Generate O_0^2 operator:
	//
	//	gamma/N_A^(3/2) * [ 1/4! sum_{k1,k2,k3} : J_k1 J_k2 J_k3 J_{-k1-k2-k3} : ]

	//	k1+k2+k3<0 terms

	for(int k1=-fourierTruncation;k1<=fourierTruncation;k1++)
	{
		for(int k2=-fourierTruncation;k2<=fourierTruncation;k2++)
		{
			for(int k3=-fourierTruncation;k3<=fourierTruncation;k3++)
			{
				int sum=k1+k2+k3;

				if(abs(sum)<=fourierTruncation && sum!=0 && k1!=0 && k2!=0 && k3!=0)
				{
					std::vector<short int> indices;
					std::vector<char> labels;
					double coefficient;

					//	Ensure here that the term is normal ordered
					//	Meaning put in descending order

					int list[4];

					list[0]=k1;
					list[1]=k2;
					list[2]=k3;
					list[3]=-sum;

					//utilities::cout.DebuggingInfo()<<"UNSORTED LIST"<<std::endl;

					//utilities::cout.DebuggingInfo()<<"0: "<<list[0]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"1: "<<list[1]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"2: "<<list[2]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"3: "<<list[3]<<std::endl;

					for(int i=0;i<4;i++)
					{
						for(int j=i+1;j<4;j++)
						{
							if(list[i]>list[j])
							{

								//	swap elements i and j
								int temp;

								temp=list[i];

								list[i]=list[j];

								list[j]=temp;
							}
						}
					}

					//utilities::cout.DebuggingInfo()<<"SORTED LIST"<<std::endl;

					//utilities::cout.DebuggingInfo()<<"0: "<<list[0]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"1: "<<list[1]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"2: "<<list[2]<<std::endl;
					//utilities::cout.DebuggingInfo()<<"3: "<<list[3]<<std::endl;

					//getchar();

					indices.push_back(list[0]);
					indices.push_back(list[1]);
					indices.push_back(list[2]);
					indices.push_back(list[3]);

					labels.push_back(colour);
					labels.push_back(colour);
					labels.push_back(colour);
					labels.push_back(colour);

					coefficient=(gamma/(24.0))/pow(nbrA0,1.5);

					hamiltonian.push_back(Term(coefficient,indices,labels));
				}
			}
		}
	}

	//	Call routine to combine terms together in the Hamiltonian
	//	if they have the same operator content
	
	ListOfTerms temp(hamiltonian);

	temp.CombineSameTerms();

	return temp;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
