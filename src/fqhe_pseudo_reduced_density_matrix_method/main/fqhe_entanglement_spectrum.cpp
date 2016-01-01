////////////////////////////////////////////////////////////////////////////////
//!
//!					\authors Ivan D Rodriguez, Simon C Davenport
//!
//!                         \date Last Modified: 25/01/2015
//!
//!  \file 
//!		This program is designed to generate (approximately) the real space 
//!		entanglement spectrum of a given fractional quantum Hall wave function.
//!
//!		The calculation performed by this code is roughly as follows.
//!
//!		We want to determine an approximation to the reduced density
//!		matrix for a quantum Hall system subdivided into parts A and B
//!		in real space or in particle space. The entanglement spectrum of
//!		the associated reduced matrix is a plot of \f$ - \f$ log(eigenvalues) vs.
//!		the angular momentum sector.
//!	
//!		To do that, we can obtain an approximation to the spectrum from a 
//!		certain "pseudo" reduced density matrix. First define:
//!
//!			\f[ \Phi(L_z^A,Z_A,Z_B) = \int d \phi exp(-2I\pi L_z^A \phi) \psi(z_1 exp(2I\pi \phi),\ldots,z_k exp(2I\pi \phi), z_k+1,\ldots z_N) \f]
//!	
//!		The Fourier transform in \f$ \Phi \f$ selects a specific angular momentum 
//!		sector in the wave function. The reduced density matrix is given by:
//!																				
//!		\f[ \rho_A(Z_A;Z_A')=\int dZ_B \Phi(L_z^A,Z_A,Z_B) \Phi^*(L_Z^A,Z'_A,Z_B) \f]
//!	
//!		where \f$ Z_A \f$ and \f$ Z_B \f$ are the sets of coordinates in the A and B subsystems
//!
//!		For the particle entanglement spectrum \f$ Z_A \f$ and \f$ Z_B \f$ can take any value
//!		but for a real space entanglement spectrum (RSES) the \f$ Z_A \f$ and \f$ Z_B \f$ 
//!		co-ordinates must lie in separate regions.
//!
//!		An approximation to rho_A is given by generating a set of rho_A 
//!		using a selected list of sample coordinate values. We shall refer to 
//!		this list as the "indexes". The same indexes are used for \f$ Z_A \f$ and \f$ Z_A' \f$,
//!		but there is no restriction for the number of rows and columns in
//!		the \f$ \rho_A \f$ to be equal. The number of indexes must be greater than both
//!		the rows and columns of the reduced density matrix. As the size of
//!		the effective reduced density matrix is increased, so its eigenvalues
//!		tend to converge to those of the true reduced density matrix. The 
//!		convergences expected to be faster for a rectangular effective \f$ \rho_A \f$
//!		(for which we then need to calculate the SVD).
//!
//!		When using this method, it is very important to select a good set of
//!		indexes. The co-ordinate samples used must have a good overlap with
//!		the angular momentum sectors for which we want to sample.
//!		The indexes are selected by taking \f$ l^i_A \f$ to be an integer closest to
//!		\f$ |z_i|^2 \f$ and then only selecting samples that satisfy:	
//!	
//!			\f$ \sum_i l^i_A = L_z^A \f$
//!
//!		Further, when two such samples share the same set of \f$ l^i_A \f$ we can
//!		discard one of them without loss of rank. Also, we don't need to
//!		include all possible \f$ l_z^A \f$ sets to get good convergence.
//!
//!		Of particular interest is the low lying part of the RSES (the largest)
//!		eigenvalues of \f$ \rho_A \f$. In that case, we can apply a further constraint
//!		that the particles in subsystem A are located in a disk of discRadius \f$ r_A \f$
//!		(since for the RSES calculation, the particles are confined to some
//!		region of space).
//!
//!			Copyright (C) 2012 - 2013 Ivan D Rodriguez, Simon C Davenport
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

#include "../program_options.hpp"
//  Declaration of program options
#include "../file_name_generator.hpp"
//  Declaration of file name generation functions
#include <boost/program_options.hpp>
//	Include program options library
#include <fftw3.h>
//	"Fastest Fourier Transform in the West" library functions
#include "../../utilities/mathematics/mt.hpp"
//	Mersenne Twistor random number generator
#include "../../utilities/wrappers/mpi_wrapper.hpp"
//  Wrapper for MPI functionality
#include "../../utilities/general/cout_tools.hpp"
//  Functions to manipulate std::cout output
//	MPI functions and variables for parallelization
#include "../../fqhe_wave_function_algorithms/fqhe_wave_function.hpp"
//	Base class for the wave function algorithms
#include "../../fqhe_wave_function_algorithms/laughlin.hpp"
//	Laughlin wave function algorithms
#include "../../fqhe_wave_function_algorithms/composite_fermion.hpp"
//	Jain and BS wave function algorithms
#include "../../fqhe_wave_function_algorithms/moore_read.hpp"
//	Moore-Read wave function algorithms
#include "../../utilities/wrappers/lapack_wrapper.hpp"
//  Wrapper for LAPACK's SVD subroutine
#include "../../utilities/general/load_bar.hpp"
//  For load bar
#include "../../utilities/general/run_script.hpp"
//  For running python scripts

#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif

///////		FUNCTION FORWARD DELCATATIONS		    ////////////////////////////

void PrintWelcomeMessage(const utilities::MpiWrapper& mpi);

boost::program_options::variables_map ParseComandLine(int argc,char *argv[],utilities::MpiWrapper& mpi);

void SeedRandomNumberGenerator(MersenneTwister* mt);

double Modulus(double a,double b);

double FindStartSector(const GeneralOptions* parameters,const FQHE::WaveFunctionData* wfData,
const FQHE::CompositeFermionData* cfData);	 

void InitializeIndexes(MersenneTwister* mt,const FQHE::WaveFunctionData* wfData,
const FQHE::CompositeFermionData* cfData,dcmplx* sphereCoords,const int nbrOrbitals);
 
int DoMetropolisSampling(const GeneralOptions* parameters,MersenneTwister* mt,
const MetropOptions* metropolis,FQHE::WaveFunction* wf,dcmplx *sphereCoords,const int nbrOrbitals); 
	   
void OrderIndexes(const int nbrSamples,const int nbr,const int nbrFixed,
int* listLzA,dcmplx* listCoord);

int ApplyLowLyingConstraint(const GeneralOptions* parameters,const int nbrSamples,
const int nbr,const int nbrFixed,int* listLzA,dcmplx* listCoord,const int nbrOrbitals);

int FilterEmptyRegionA(const GeneralOptions* parameters,const int nbrSamples,
const int nbr,const int nbrFixed,int* listLzA,dcmplx* listCoord);

void EvaluateFourierTransform(const GeneralOptions* parameters,FQHE::WaveFunction* wf,
const dcmplx* sphereCoords,dcmplx* waveFunctionValue);

void MakeRoughPlot(const GeneralOptions* parameters,const FQHE::WaveFunctionData* wfData,
const FQHE::CompositeFermionData* cfData,
const utilities::MpiWrapper& mpi);

//  Declare an instance of the global utilities::cout class
utilities::Cout utilities::cout;

//  Declare an instance of the global mpi wrapper class
utilities::MpiWrapper mpi(utilities::cout);

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief The main function performs a number of tasks:
//!     -   It parses the command line and sets all internal parameters
//!         including synchronizing between different nodes when using multiple
//!         processes.
//!     -   It calls a function to determine the lowest lzA configuration associated
//!         with a certain FQHE state with a certain particle cut. 
//!     -   It optionally calculates a sample set of co-ordinates with which
//!         to approximately evaluate the reduced density matrix of a FQHE
//!         state with a specified real space cut
//!     -   It optionally evaluates a pseudo reduced density matrix for a specified
//!         FHQE state (one matrix for each sector).
//!     -   It optionally determines the singular values of the pseudo-reduced
//!         density matrix, which are approximations to the true spectrum
//!         of the reduced density matrix of that system.
//!
//!	SUB SECTIONS:
//!
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	//////  	START PARALLEL PROCESSES      //////////////////////////////////////

	mpi.Init(argc,argv);
	
	//////  	PRINT WELCOME MESSAGAE      ////////////////////////////////////////

    PrintWelcomeMessage(mpi);
	
	//////		DECLARE DATA SCTURCUTRES/CLASSES      //////////////////////////
    //  (these declarations are replicated on all nodes) 

    MersenneTwister mt;				//	Random number generator object called mt
    FQHE::WaveFunction *wf=0;	    //  Object to hold wave function data/routines
    FQHE::WaveFunctionData wfData;	//	Data struct to store wave function data	
    FQHE::CompositeFermionData cfData;
                                    //	Data struct to store composite fermion data		
    GeneralOptions parameters;      //  General command line options list
    MetropOptions metropolis;       //  Metropolis parameters                                                             
    boost::program_options::variables_map optionList;   
                                    // Command line argument list

    //////      PARSE COMMAND LINE OPTIONS      ////////////////////////////////////

    optionList = ParseComandLine(argc,argv,mpi);

	//////  	Initialize and synchronize variables      //////////////////////

    wfData.InitFromCommandLine(&optionList,mpi);
	cfData.InitFromCommandLine(&optionList,mpi);
	parameters.InitFromCommandLine(&optionList,mpi);
	metropolis.InitFromCommandLine(&optionList,mpi);

    //	Default nbrFixed to half the number of particles
    if(parameters.nbrFixed==0)	parameters.nbrFixed=(int)(0.5*wfData.nbr);

	//////  	SELECT WAVE FUNCTION TYPE (BASED ON PROGRAM OPTIONS)      //////////
	
	if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
		if(FQHE::_COMPOSITE_FERMION_ == wfData.type|| FQHE::_BONDERSON_SLINGERLAND_ == wfData.type)
		{
			wf = new FQHE::CompositeFermion(&wfData,&cfData);
		}
		else if(FQHE::_LAUGHLIN_ == wfData.type)
		{
			wf = new FQHE::Laughlin(&wfData);	
		}
		else if(FQHE::_MOORE_READ_ == wfData.type)
		{
			wf = new FQHE::MooreRead(&wfData);	
		}
    
        //  Print out program parameters information

		utilities::cout.MainOutput()<<"\tFile in/out path is \t\t "<<parameters.path<<std::endl<<std::endl;

		utilities::cout.MainOutput()<<"\tES calculation data"<<std::endl<<std::endl;

		if(parameters.choice==1)
		{
			utilities::cout.MainOutput()<<"\tNo. thermalizing Metropolis function calls\t"<<metropolis.nbrTherms<<std::endl;
			utilities::cout.MainOutput()<<"\tNo. sampling Metropolis calls\t\t\t"<<metropolis.nbrSamples<<std::endl;
			utilities::cout.MainOutput()<<"\tMetropolis steps between each sample\t\t"<<metropolis.iter<<std::endl;
			utilities::cout.MainOutput()<<"\tNumber fixed\t\t\t\t\t"<<parameters.nbrFixed<<std::endl;
			utilities::cout.MainOutput()<<"\tReal cut\t\t\t\t\t"<<parameters.realCut<<std::endl;
			utilities::cout.MainOutput()<<"\tLzA\t\t\t\t\t\t"<<parameters.lzA<<std::endl;
			utilities::cout.MainOutput()<<"\tRenorm factor\t\t\t\t\t"<<parameters.renormFactor<<std::endl;
		}
		else
		{
			utilities::cout.MainOutput()<<"\tNumber fixed\t\t\t"<<parameters.nbrFixed<<std::endl;
			utilities::cout.MainOutput()<<"\tReal cut\t\t\t"<<parameters.realCut<<std::endl;
			utilities::cout.MainOutput()<<"\tLzA\t\t\t\t"<<parameters.lzA<<std::endl;
			utilities::cout.MainOutput()<<"\tNo. sectors\t\t\t"<<parameters.nbrSect<<std::endl;
			utilities::cout.MainOutput()<<"\tNo. Fourier components\t\t"<<parameters.nbrFourier<<std::endl;
			utilities::cout.MainOutput()<<"\tRho_A dimensions\t\t"<<parameters.rows<<"x"<<parameters.columns<<std::endl;
			utilities::cout.MainOutput()<<"\tRenorm factor\t\t\t"<<parameters.renormFactor<<std::endl<<std::endl;;
		}
	}	
	else	// FOR ALL OTHER NODES
	{
	    //  Turn off command line output
		utilities::cout.SetVerbosity(0);

		if(FQHE::_COMPOSITE_FERMION_ == wfData.type|| FQHE::_BONDERSON_SLINGERLAND_ == wfData.type)
		{
			wf = new FQHE::CompositeFermion(&wfData,&cfData);
		}
		else if(FQHE::_LAUGHLIN_ == wfData.type)
		{
			wf = new FQHE::Laughlin(&wfData);	
		}
		else if(FQHE::_MOORE_READ_ == wfData.type)
		{
			wf = new FQHE::MooreRead(&wfData);	
		}
	}

	//  MPI sync verbosity levelz
    
    utilities::cout.MpiSync(0,mpi.m_comm);

	//	Barrier waits for all nodes to catch up to this point
	MPI_Barrier(mpi.m_comm);
	
	//////      EXECUTE DIFFERENT PROGRAM MODES     ////////////////////////////////

switch(parameters.choice)
{
	case 1: 
	{	
        ////////////////////////////////////////////////////////////////////////////////
        //!	- GENERATE INDEXES OF RHO_A.
		//!		
        //!		In this section of the code we generate a file called "indexes" that
        //!		contains a list of coordinate samples. These samples are chosen to have
        //!		a high overlap with a particular angular momentum sector set by lzA.
        //!		This is done by metropolis. sampling and selecting only those 
        //!		samples where \f$ l^i_A = |z|^2 \f$ (rounded to the nearest integer) such 
        //!		that \f$ \sum_i l^i_A=l_z^S \f$.
        //!	
        //!		Duplicate sets of \f$ l^i_A \f$ can be discarded.
        //!
        ////////////////////////////////////////////////////////////////////////////////

        //////      OPEN INDEXES FILE FOR WRITING      /////////////////////

        std::ofstream f_coord;

        if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
			f_coord.open(GenerateIndexesFileName(&parameters,&cfData,&wfData).c_str(), std::ios::out | std::ios::binary);

			if(!f_coord.is_open())
			{
				std::cerr<<"\tFile: "<<GenerateIndexesFileName(&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
				mpi.m_exitFlag = true;
			}
			else
			{
			    utilities::cout.MainOutput()<<"\n\tFind indexes file at "<<GenerateIndexesFileName(&parameters,&cfData,&wfData)<<std::endl<<std::endl;
			}
	    }
	    
	    mpi.ExitFlagTest();

		if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
            //////      CALCULATE LOWEST LZA STATE      ////////////////////////////////////

			utilities::cout.SecondaryOutput()<<"\tLowest lzA EWF state\t\t\t\t"<<FindStartSector(&parameters,&wfData,&cfData)<<std::endl;

            //  Number of orbitals
            const int nbrOrbitals = (int)(2.0*wfData.monopoleStrength)+1;

            utilities::cout.SecondaryOutput()<<"\tNo. orbitals\t\t\t\t\t"<<nbrOrbitals<<std::endl;
    
			//////       SEED THE RANDOM NUMBER GENERATOR      /////////////////

            SeedRandomNumberGenerator(&mt);

			//////      ALLOCATE MEMORY TO STORE INDEXES      //////////////////

			dcmplx* sphereCoords = new dcmplx[wfData.nbr];
			dcmplx* listCoord    = new dcmplx[wfData.nbr*metropolis.nbrSamples];
			int* listLzA         = new int[wfData.nbr*metropolis.nbrSamples];

			//////      GENERATE INDEXES      //////////////////////////////////
			
			generate_indexes:   //  GoTo marker

            utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;

			utilities::cout.MainOutput()<<"\n\tGenerating indexes for calculation of rho_A\n\n\tPress any key to continue"<<std::endl;
			getchar();

            //////      SPECIFY STARTING CONFIGURATION      ////////////////////

			InitializeIndexes(&mt,&wfData,&cfData,sphereCoords,nbrOrbitals);

			utilities::cout.DebuggingInfo()<<"print init co-ords"<<std::endl;
			
			for(int i=0;i<wfData.nbr;i++)
			{
		        utilities::cout.DebuggingInfo()<<sphereCoords[i];
			}
			
			utilities::cout.DebuggingInfo()<<std::endl;

			//////      THERMALIZE MONTE CARLO      ////////////////////////////
			
			int restartFlag=0;	//	use to restart the iterations if many steps are rejected

            utilities::cout.MainOutput()<<"\tThermalizing ... "<<std::endl;

			for (int i=0; i<metropolis.nbrTherms; i++)
			{
				utilities::cout.SecondaryOutput()<<"\tThermalization "<<i;
				DoMetropolisSampling(&parameters,&mt,&metropolis,wf,sphereCoords,nbrOrbitals);
			}

			//////      MAIN MONTE CARLO SAMPLING      /////////////////////////
			
			utilities::cout.MainOutput()<<"\tPerforming sampling... "<<std::endl;
			
			for (int i=0; i<metropolis.nbrSamples; i++)
			{
				utilities::cout.SecondaryOutput()<<"\tIteration "<<i;
				
				restartFlag += DoMetropolisSampling(&parameters,&mt,&metropolis,wf,sphereCoords,nbrOrbitals);

				if(i%20==0)
				{
					//	check for all rejected steps every 20 iterations

					if(restartFlag==0)
					{
						utilities::cout.MainOutput()<<"\tRESTART METROPOLIS SAMPLING"<<std::endl;
						goto generate_indexes;  
					}

					restartFlag=0;
				}

                //  Store sample in a list

				for (int j=0;j<wfData.nbr;j++)
				{  
					listCoord[i*wfData.nbr+j]=sphereCoords[j];

					//	Extract the integer closest to |z_i|^2
					listLzA[i*wfData.nbr+j]=(int)round(wfData.monopoleStrength*(1.0 + cos(real(sphereCoords[j]))) );
				}
			}        

			//////      PUT INDEXES IN ORDER OF THEIR LZA      /////////////////
		
		    OrderIndexes(metropolis.nbrSamples,wfData.nbr,parameters.nbrFixed,listLzA,listCoord);

			//////      APPLY EXTRA CONSTRAINT FOR LOW-LYING RSES CASE      ////
            
            int remainingNbrSamples = ApplyLowLyingConstraint(&parameters,metropolis.nbrSamples,
                               wfData.nbr,parameters.nbrFixed,listLzA,listCoord,nbrOrbitals);
			
			//////      LOOK AT SAMPLES WHERE REGION A IS EMPTY      ///////////

            remainingNbrSamples = FilterEmptyRegionA(&parameters,remainingNbrSamples,wfData.nbr,
                                    parameters.nbrFixed,listLzA,listCoord);

			//////      WRITE INDEXES TO A FILE      ///////////////////////////

			//	The first entry in the file is the number of sets of samples

			f_coord.write(reinterpret_cast<char*>(&remainingNbrSamples),sizeof(int));

			//	Write the co-ordinate samples to the file

            f_coord.write(reinterpret_cast<char*>(listCoord),remainingNbrSamples*wfData.nbr*sizeof(dcmplx));
            f_coord.close();

			for (int i=0;i<remainingNbrSamples;i++)
			{
				for (int j=0;j<wfData.nbr;j++)
				{
					utilities::cout.SecondaryOutput()<<"\t"<<listLzA[i*wfData.nbr+j]<<" ";
				}

				utilities::cout.SecondaryOutput()<<std::endl;
			}

			utilities::cout.MainOutput()<<"\n\tNumber of usable samples: "<<remainingNbrSamples<<std::endl;
			
			//	deallocate memory
            delete[] listLzA;
            delete[] listCoord;
            delete[] sphereCoords;
            
            utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;
		}
		
	break;
	}
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//	
	case 2:
	{ 
	    //////////////////////////////////////////////////////////////////////////////////
	    //!  - COMPUTE APPROXIMATION TO RHO_A.
		//!
	    //!  	Evaluate the pseudo-reduced density matrix for each sector.
	    //////////////////////////////////////////////////////////////////////////////////

        std::ifstream f_coord;

		if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
		    utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;
		
			utilities::cout.MainOutput()<<"\n\tComputing rho_A\n\n\tPress any key to continue"<<std::endl;
			getchar();

			//	Read in indexes

			f_coord.open(GenerateIndexesFileName(&parameters,&cfData,&wfData).c_str(),std::ios::in | std::ios::binary);

			if(!f_coord.is_open())
			{
				std::cerr<<"\tFile: "<<GenerateIndexesFileName(&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
				mpi.m_exitFlag = true;
			}
        }
        
        mpi.ExitFlagTest();
        	
	    int dim1;
        
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
		    utilities::cout.MainOutput()<<"\tReading in index file..."<<std::endl;
		
			f_coord.read(reinterpret_cast<char*>(&dim1),sizeof(int));

			if(parameters.rows>dim1 || parameters.columns>dim1)
			{
				std::cerr<<"\tThe dimension of rho_A ("<<parameters.rows<<","<<parameters.columns<<") exceeds the number of indexes ("<<dim1<<") in "<<GenerateIndexesFileName(&parameters,&cfData,&wfData)<<std::endl;
				mpi.m_exitFlag=true;
			}
        }
        
        mpi.ExitFlagTest();
		
	    //	ALL NODES - Allocate memory to store pseudo reduced density matrix

		dcmplx* rowsRho     = new dcmplx[(parameters.rows+2)*parameters.nbrFixed];
		dcmplx* columnsRho  = new dcmplx[(parameters.columns+2)*(wfData.nbr-parameters.nbrFixed)];
		
		if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{	
			//	Read in column and row samples
			
			int counter1=0,counter2=0;
		   
			for(int i=0; i<dim1; i++)
			{
				for(int j=0; j<wfData.nbr; j++)
				{
                    dcmplx temp;

					f_coord.read(reinterpret_cast<char*>(&temp),sizeof(dcmplx));

					if(Modulus(i,(int) (dim1/(parameters.rows-1)))==0 && j<parameters.nbrFixed && counter1 != parameters.rows)
					{           
						rowsRho[counter1*parameters.nbrFixed+j]=temp;
					} 

					if(Modulus(i,(int) (dim1/(parameters.columns-1)))==0 && j>=parameters.nbrFixed && counter2 != parameters.columns)
					{                             
						columnsRho[counter2*(wfData.nbr-parameters.nbrFixed)+(j-parameters.nbrFixed)]=temp;
					}
				}
				
				if(Modulus(i,(int) (dim1/(parameters.rows-1)))==0 && counter1 != parameters.rows - 1)		counter1++;
				if(Modulus(i,(int) (dim1/(parameters.columns-1)))==0 && counter2 != parameters.columns-1)	counter2++;
			}
			
			f_coord.close();
		}

		//	Barrier waits for all nodes to catch up to this point
		MPI_Barrier(mpi.m_comm);
		
		//	MPI - SYNCHRONISE VARIABLE VALUES WITH MASTER NODE VALUES
		MPI_Bcast(rowsRho,(parameters.rows+2)*parameters.nbrFixed,MPI_DOUBLE_COMPLEX,0,mpi.m_comm);
		MPI_Bcast(columnsRho,(parameters.columns+2)*(wfData.nbr-parameters.nbrFixed),MPI_DOUBLE_COMPLEX,0,mpi.m_comm);
		
		//	Barrier waits for all nodes to catch up to this point
		MPI_Barrier(mpi.m_comm);

        //////	    CALCULATE FOURIER TRANSFORMS      //////////////////////////

		//  Determine the distribution of tasks amongst the available processors

		mpi.DivideTasks(mpi.m_id,parameters.rows,mpi.m_nbrProcs,
		                       &mpi.m_firstTask,&mpi.m_lastTask,true);

		MPI_Barrier(mpi.m_comm);

		//	Allocate memory to store co-ordinates and wave function values

		dcmplx* waveFunctionValue = new dcmplx[parameters.nbrSect];
		dcmplx* sphereCoords      = new dcmplx[wfData.nbr];

		//	Output file for sector data

        std::ofstream f_sector;

		f_sector.open(GenerateNodeFileName(mpi.m_id,&parameters,&cfData,&wfData).c_str(), std::ios::out | std::ios::binary);
        
        utilities::cout.MainOutput()<<"\tCreating output file at "<<GenerateNodeFileName(mpi.m_id,&parameters,&cfData,&wfData)<<std::endl;
        
        MPI_Barrier(mpi.m_comm);
        
        utilities::LoadBar progress;
        
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
		    utilities::cout.AdditionalInfo()<<std::endl;
		
            progress.Initialize(mpi.m_lastTask - mpi.m_firstTask+1);
        }
        
		for(int i=mpi.m_firstTask; i<=mpi.m_lastTask; ++i)
		{		
			for(int k=0; k<parameters.nbrFixed; ++k)
			{
				sphereCoords[k]=rowsRho[i*parameters.nbrFixed+k];

				//std::cout<<"NODE "<<mpi.m_id<<" "<<sphereCoords[k]<<std::endl;
			}
			
			for(int j=0; j<parameters.columns; ++j)
			{
				for (int k=parameters.nbrFixed;k<wfData.nbr;++k)
				{
					sphereCoords[k]=columnsRho[j*(wfData.nbr-parameters.nbrFixed)+k-parameters.nbrFixed];

					//std::cout<<"NODE "<<mpi.m_id<<" "<<sphereCoords[k]<<std::endl;
				}
				
				//std::cout<<"----------------------"<<std::endl;
				
				EvaluateFourierTransform(&parameters,wf,sphereCoords,waveFunctionValue);
				
				for(int k=0;k<parameters.nbrSect;++k)
				{
					//std::cout<<"i "<<i<<"j "<<j<<"k "<<k<<" "<<waveFunctionValue[k]<<std::endl;getchar();
					
					//	write reduced density matrix data to a file

					if(std::isnan(real(waveFunctionValue[k])) || std::isnan(imag(waveFunctionValue[k])))
					{
						std::cerr<<"\tWARNING: nan waveFunctionValue value detected on node "<<mpi.m_id<<std::endl;
						getchar();
					}
				}
				
				f_sector.write(reinterpret_cast<char*>(waveFunctionValue),parameters.nbrSect*sizeof(dcmplx));
			}
			
			//utilities::cout.SecondaryOutput()<<"\tEvaluated row: "<<i<<" of "<<mpi.m_firstTask<<" to "<<mpi.m_lastTask<<std::endl;
			
            //	Barrier waits for all nodes to catch up to this point

            MPI_Barrier(mpi.m_comm);
            
			if(0 == mpi.m_id)	// FOR THE MASTER NODE
            {
                progress.Display(i+1);
            }
		}
        
        if(0 == mpi.m_id)	// FOR THE MASTER NODE
        {
            utilities::cout.SecondaryOutput()<<"\n\tWAITING FOR ALL NODES TO CATCH UP TO THIS POINT..."<<std::endl;
        }

		//	Barrier waits for all nodes to catch up to this point

		MPI_Barrier(mpi.m_comm);

		//	close file
		f_sector.close();
		
		if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
            utilities::cout.AdditionalInfo()<<std::endl;
        }
        
        MPI_Barrier(mpi.m_comm);

		utilities::cout.MainOutput()<<"\tCompleted output file: "<<GenerateNodeFileName(mpi.m_id,&parameters,&cfData,&wfData)<<std::endl;

        MPI_Barrier(mpi.m_comm);

		//	Deallocate memory
		delete[] rowsRho;	
		delete[] columnsRho;
		delete[] waveFunctionValue;
		delete[] sphereCoords;
		
		if(0 == mpi.m_id)	// FOR THE MASTER NODE
	    {
		    utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;
		}
		
		break;
	}
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
	case 3:
	{
	    //////////////////////////////////////////////////////////////////////////////////
	    //! - COMPUTE RHO_A SPECTRUM.
	    //////////////////////////////////////////////////////////////////////////////////

        //////      REORGANIZE .tmp FILE DATA      /////////////////////////////
        //  Generate new files containing the pseudo-reduced density matrix
        //  for each sector
		
		if(0 == mpi.m_id)	// FOR THE MASTER NODE	//
		{
		    utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;
		
			utilities::cout.MainOutput()<<"\n\tChecking for process.tmp files..."<<std::endl<<std::endl;
			
			//	determine number of available files to reconstruct
			
			int nbrProcs=0;
			std::fstream f_sector;

			while(true)
			{
				f_sector.open(GenerateNodeFileName(nbrProcs,&parameters,&cfData,&wfData).c_str(), std::ios::in | std::ios::binary);
				
				if(!f_sector.is_open())
				{	 
						utilities::cout.SecondaryOutput()<<"\tFile: "<<GenerateNodeFileName(nbrProcs,&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
						break;
				}
				else
				{
					utilities::cout.SecondaryOutput()<<"\tDetected file: "<<GenerateNodeFileName(nbrProcs,&parameters,&cfData,&wfData)<<std::endl;
					nbrProcs++;
					//	close file
					f_sector.close();
				}
			}
			
			if(nbrProcs==0)
			{
				utilities::cout.MainOutput()<<"\tNo valid .tmp files found. Attempting to use sector.dat files instead."<<std::endl;
				goto calculateSVD;
			}
			
			utilities::cout.MainOutput()<<"\tReading in and re-organising sector data files (takes a few minutes!)."<<std::endl;

            utilities::cout.SecondaryOutput()<<std::endl;

			//	Allocate memory to store data from all sectors

			utilities::cout.SecondaryOutput()<<"\tAllocating "<<(long int)parameters.rows*parameters.columns*parameters.nbrSect*sizeof(dcmplx)/(1024*1024)<<" mb of memory to store data\n"<<std::endl;

		    dcmplx* reducedDensityMatrixList = new dcmplx[(long int)parameters.rows*parameters.columns*parameters.nbrSect];

			//	Read in the full set of data

			for(int proc=0;proc<nbrProcs;proc++)
			{
				//	read in file generated by each process to reconstruct one sector at a time

				f_sector.open(GenerateNodeFileName(proc,&parameters,&cfData,&wfData).c_str(), std::ios::in | std::ios::binary);

				if(!f_sector.is_open())
				{
						std::cerr<<"\tFile: "<<GenerateNodeFileName(proc,&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
						exit(EXIT_FAILURE);
				}

				//	determine mpi.m_firstTask and mpi.m_lastTask for a given proc

				mpi.DivideTasks(proc,parameters.rows,nbrProcs,&mpi.m_firstTask,&mpi.m_lastTask,false);

                for(int i=mpi.m_firstTask;i<=mpi.m_lastTask;++i)
				{
					for (int j=0;j<parameters.columns;++j)
					{
					    f_sector.read(reinterpret_cast<char*>(reducedDensityMatrixList+(long int)(j*parameters.rows+i)*parameters.nbrSect),parameters.nbrSect*sizeof(dcmplx));
					}
				}

				f_sector.close();

				utilities::cout.SecondaryOutput()<<"\tCompleted reading in file: "<<GenerateNodeFileName(proc,&parameters,&cfData,&wfData)<<std::endl;
			}

			utilities::cout.SecondaryOutput()<<std::endl;

			//	generate a new file to contain the reduced density matrix for each sector

			for(int k=0;k<parameters.nbrSect;++k)
			{
				f_sector.open(GenerateSectorFileName(k,&parameters,&cfData,&wfData).c_str(), std::ios::out | std::ios::binary);

                dcmplx* p_matrix = reducedDensityMatrixList + k;

				for(int i=0; i<parameters.rows;++i)
				{
					for(int j=0; j<parameters.columns;++j)
					{
						f_sector.write(reinterpret_cast<char*>(p_matrix+(long int)(i*parameters.columns+j)*parameters.nbrSect),sizeof(dcmplx));
					}
				}

				f_sector.close();

				utilities::cout.SecondaryOutput()<<"\tCompleted writing file: "<<GenerateSectorFileName(k,&parameters,&cfData,&wfData)<<std::endl;
			}

			//	deallocate memory

			delete[] reducedDensityMatrixList;
		}

        //////  	CALCULATE SVD     //////////////////////////////////////////

		calculateSVD:   //  GoTo marker

		if(0 == mpi.m_id)	// FOR THE MASTER NODE
		{
		    utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;
		
			utilities::cout.MainOutput()<<"\n\tCalculating SVD"<<std::endl;
			
			utilities::cout.SecondaryOutput()<<"\n\tRho_A size: \t"<<parameters.rows<<"x"<<parameters.columns<<std::endl<<std::endl;

			dcmplx* reducedDensityMatrix = new dcmplx[parameters.rows*parameters.columns];
			double* singularValues       = new double[parameters.rows];

			for(int k=0;k<parameters.nbrSect;k++)
			{
				utilities::cout.SecondaryOutput()<<"\t  Sector "<<k+parameters.lzA<<std::endl<<std::endl;
				
				//	read in reduced density matrix for each sector

                std::ifstream f_sector;

				f_sector.open(GenerateSectorFileName(k,&parameters,&cfData,&wfData).c_str(), std::ios::in | std::ios::binary);
				
				if(!f_sector.is_open())
				{	 
					std::cerr<<"\tFile: "<<GenerateSectorFileName(k,&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
					exit(EXIT_FAILURE);
				}
				
				//	read in reduced density matrix of a given sector
				
				f_sector.read(reinterpret_cast<char*>(reducedDensityMatrix),parameters.rows*parameters.columns*sizeof(dcmplx));

                std::ofstream f_spectrum;

				f_spectrum.open(GenerateEigenvaluesFileName(k,&parameters,&cfData,&wfData).c_str(),std::ios::out);

				if(!f_spectrum.is_open())
				{
					std::cerr<<"\tFile: "<<GenerateEigenvaluesFileName(k,&parameters,&cfData,&wfData)<<" not found !!! "<<std::endl<<std::endl;
					exit(EXIT_FAILURE);
				}

				utilities::linearAlgebra::SingularValueDecomposition<dcmplx>(
				    reducedDensityMatrix,singularValues,parameters.rows,parameters.columns);

				f_spectrum.precision(15);

                //  Output spectrum data to a file

				for (int i=0;i<parameters.rows;i++)
				{
					utilities::cout.SecondaryOutput()<<"\t"<<2.0*log(singularValues[i])<<std::endl;
					f_spectrum<<2.0*log(singularValues[i])<<"\n";
				}
				
				utilities::cout.SecondaryOutput()<<std::endl;
				
				f_spectrum.close();
				f_sector.close();

			}

			//	Deallocate memory
			delete[] reducedDensityMatrix;
			delete[] singularValues;
			
			//////  	PLOT RSES     //////////////////////////////////////////
			
			MakeRoughPlot(&parameters,&wfData,&cfData,mpi);
		}

		break;
        
    }
	default:
			printf("\tInvalid Option !!! \n"); 
			break;
	} 

    //////  	TERMINATE PROGRAM      /////////////////////////////////////////////

    delete wf;

    return 0;
}

//////	    FUNCTION DECLARATIONS      /////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Print a message to display at start of program
//!
////////////////////////////////////////////////////////////////////////////////

void PrintWelcomeMessage(const utilities::MpiWrapper& mpi)
{
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    std::string message=
	    "\n\n---------------------------------------------"
	    "---------------------------------------------\n"
	    "---------------------------------------------"
	    "---------------------------------------------\n\n"
	    "\t\t\tEntanglement Spectrum Calculator\n\n"
	    "\t\t Authors: Ivan Rodriguez & Simon C Davenport	\n\n"
	    "\t\t\t   Last Modified: 25/01/2015\n\n"
	    "\tThis program is designed to generate the real space entanglement\n"
	    "\tspectrum of a given fractional quantum Hall trial wave function.  \n"
	    "\tVersion 1.1 onwards uses MPI parallelization \n\n"
	    "\t\tCopyright (C) 2013--2015 IVAN D. RODRIGUEZ & SIMON C DAVENPORT.\n\n"
	    "\tThis program comes with ABSOLUTELY NO WARRANTY.  \n"
	    "\tThis is free software, and you are welcome to redistribute it \n"
	    "\tunder certain conditions. See http://www.gnu.org/licenses/.\n\n"
	    "---------------------------------------------"
	    "---------------------------------------------\n"
	    "---------------------------------------------"
	    "---------------------------------------------\n";

	    utilities::cout.MainOutput()<<message<<std::endl;
    }
    
	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Converts command line arguments into parameters values in global
//!	data structures. 
//!
//!	Note: this function sets values for the global parameters in the optionSet
//!	and metropOptions namespaces. 
//!
////////////////////////////////////////////////////////////////////////////////

boost::program_options::variables_map ParseComandLine(
	int argc,		//!<	Number of characters to parse
	char *argv[],	//!<	Character array to parse
	utilities::MpiWrapper& mpi)   //!<    Instance of the mpi wrapper class
{
    namespace po = boost::program_options;

    po::variables_map vm;

    if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    //  Declare general program options
	
	    po::options_description generalOpt("General Options");
        generalOpt.add_options()
        ("help,h", 
         "Display this message\n")
        ("mode,m", po::value<int>()->default_value(0),
         "Set the program mode: \n\t  1: GENERATE INDEXES FOR CALCULATION OF RHO_A \n\t  2: COMPUTE RHO_A FOR SELECTED SECTORS \n\t  3: ARRANGE OUTPUT FILES AND COMPUTE SPECTRUM\n")
        ("verbose,v",po::value<int>()->default_value(1),
         "Set a value for the verbosity level:\n\t 0 output off (after command line parsed) \n\t 1 print brief information \n\t 2 print more detailed information \n\t 3 also print loading bars\n\t 4 also print some debugging messages.");
	
	    //	Combine all the option groups into one
	    
	    po::options_description allOpt("\tThe program input options are as follows");
	    allOpt.add(generalOpt).add(GetEntanglementSpectrumOptions()).add(GetMetropolisOptions()).add(FQHE::GetWaveFunctionOptions()).add(FQHE::GetCompositeFermionOptions());
	
        try
        {
    
            //	Map the command line argument onto the option set
            po::store(po::parse_command_line(argc, argv,allOpt), vm);
	
            //	Respond to help option declaration
            if (vm.count("help")) 
            {
                utilities::cout.MainOutput() << allOpt << "\n";
                mpi.m_exitFlag=true;
            }
            
            po::notify(vm);
	    
        }
        catch(po::error& e)
        {
            utilities::cout.MainOutput()<<allOpt<<std::endl;
            
            std::cerr<<utilities::cout.HyphenLine()<<std::endl;
            
            std::cerr<<std::endl<<"\tERROR:\t"<<e.what()<<std::endl;
            
            std::cerr<<std::endl<<utilities::cout.HyphenLine()<<std::endl;
            
            mpi.m_exitFlag = true;
        }
    }
    
    mpi.ExitFlagTest();
    
    //  Set global verbosity level
	
	if(0 == mpi.m_id)	// FOR THE MASTER NODE
	{
	    utilities::cout.SetVerbosity(vm["verbose"].as<int>());
    }
    
    //  MPI sync verbosity level
    
    utilities::cout.MpiSync(0,mpi.m_comm);

    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n\tRun with -h option to see program options"<<std::endl;
        utilities::cout.MainOutput()<<"\n\tSet verbosity level with -v option\n"<<std::endl;
    }
    
	return vm;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Seed the random number generator
//!
////////////////////////////////////////////////////////////////////////////////

void SeedRandomNumberGenerator(
    MersenneTwister* mt)    //!<    Address of Mersenne Twister random number 
                            //!     generator object
{
    //	First seed the standard random number generator with PBS job id or
	//  system time

    int seed;

	if(getenv("PBS_JOBID")==NULL)	seed=time(0);
	else							seed=atoi(getenv("PBS_JOBID"));

	srand(seed);

	utilities::cout.MainOutput()<<"\tRandom seed\t\t\t\t\t"<<seed<<std::endl;

	//	Now seed the Mersenne Twister random number generator using
	//	an array of these standard random numbers
	{
		const int INIT_ARRAY_SIZE=200;

		long unsigned int initArray[INIT_ARRAY_SIZE];

		for(int i=0;i<INIT_ARRAY_SIZE;i++)
		{
			initArray[i]=rand();
		}

		mt->init_by_array(initArray,INIT_ARRAY_SIZE);
		//	Now the mt random number generator is initialized
		//	generate random doubles in [0,1) with mt.random()
		//	generate random integers with mt.genrand_int31();
	}
	
	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief This function returns a modulo b (the same as the % operator)
//!
//! \return a modulo b
////////////////////////////////////////////////////////////////////////////////

double Modulus(
    double a,   //!<  Maximum value that can be returned
    double b)   //!<  Value that you want to take modulo of  
{
	return a - b *floor(a/b);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief This function determines the value of lzA associated with the
//! lowest lzA sector that appears in the entanglement spectrum (RSES).
//!
//! This is deduced by considering the orbital angular momentum
//!	of the CF-EWF quasihole wave function with the lowest orbitals occupied
//!
//! - Laughlin case
//! 
//!     Example quasihole diagram of lowest lzA state
//!	\verbatim
//!     X X X X X _ _ _ _ _ * X X X X X _ _ _ _ _ * ...
//!	\endverbatim
//!
//! - Composite fermion case
//!
//!     Example quasihole diagram of lowest lzA state
//!	\verbatim
//!     X X X _ _ _  * X X X X X _ _ _ _ _ * X X X X X _ _ _ _ _
//!       X X _ _
//!	\endverbatim
//!
//! \return The lzA value of the lowest sector in the RSES
//!
//! \todo Update to include Moore-Read and BS cases
//!
////////////////////////////////////////////////////////////////////////////////

double FindStartSector(
    const GeneralOptions* parameters,           //!<    Pointer to parameter set
    const FQHE::WaveFunctionData* wfData,       //!<    A given wave function data set
    const FQHE::CompositeFermionData* cfData)   //!<    A given composite fermion data set
{
    int branchLza=0;       //   This variable contains an integer which is twice
                           //   the orbital angular momentum

    if(wfData->type==FQHE::_LAUGHLIN_)		//	Laughlin case
    {
	    for(int i=0;i<parameters->nbrFixed;i++)
	    {
		    branchLza+=-(wfData->nbr-1)+2*i;
	    }

	    branchLza*=wfData->jastrowExponent;
    }

    if(wfData->type==FQHE::_COMPOSITE_FERMION_)		//	Composite fermion case
    {
	    //	Flux attached part

	    for(int i=0;i<parameters->nbrFixed;i++)
	    {
		    branchLza+=-(wfData->nbr-1)+2*i;
	    }

	    branchLza*=2*wfData->jastrowExponent;

	    //	Composite fermion LLs part

	    int lzList[wfData->nbr];    //	list of lzA values of each CF nbrLLs orbital
	    int *p_list=lzList;         //  pointer to a value in the list

	    double Q=(double)((int)wfData->nbr-cfData->LLup*cfData->LLup-cfData->LLdown*cfData->LLdown)/(2*cfData->LLup+2*cfData->LLdown);

	    for(int l=0;l<cfData->LLup;l++)
	    {
		    int degen=(int)(2.0*(Q+l))+1;

		    for(int i=0;i<degen;i++)
		    {
			    *(p_list)=-(degen-1)+2*i;
			    p_list++;
		    }
	    }

	    //	Sort the list in terms of angular momentum values

	    for(int i=0;i<wfData->nbr;i++)
	    {
		    for(int j=i;j<wfData->nbr;j++)
		    {
			    if(lzList[i]>lzList[j])
			    {
				    int swap=lzList[i];
				    lzList[i]=lzList[j];
				    lzList[j]=swap;
			    }
		    }
	    }

	    //	Occupy the lowest lzA set of orbitals

	    if(cfData->NEF)
	    {
		    for(int i=0;i<parameters->nbrFixed;i++)
		    {
			    branchLza-=lzList[i];
		    }
	    }
	    else
	    {
		    for(int i=0;i<parameters->nbrFixed;i++)
		    {
			    branchLza+=lzList[i];
		    }
	    }
    }

    return (double)branchLza/2.0;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

////////////////////////////////////////////////////////////////////////////////
//! \brief This function contains an initial set of co-ordinates associated
//! with each type of wave function up to a given number of particles.
//!
//!	To obtain a set of samples with a good overlap with a particular 
//!	angular momentum sector, it is necessary to carefully select an 
//!	initial set of co-ordinates. These choices depend on the wave 
//! function type (initial co-ordinate arrays are explicitly written
//! into the code).
//!
//! If a particular case is not found here then this function instead
//! checks for a set of initial co-ordinates contained in the file
//! called "initCoeffs.dat". If that file does not exist then an error 
//! message is produced.
//!
//! \todo Update to include starting co-ordinates for Moore-Read and 
//! all positive effective field BS wave functions
//!
////////////////////////////////////////////////////////////////////////////////

void InitializeIndexes(
    MersenneTwister* mt,                     //!<    Pointer to random number generator object
    const FQHE::WaveFunctionData* wfData,    //!<    A given wave function data set
    const FQHE::CompositeFermionData* cfData,//!<    A given composite fermion data set      
    dcmplx* sphereCoords,        	   		 //!<    Memory address of a list of theta,phi co-ordinates
    const int nbrOrbitals)                   //!<    Number of orbitals (defined as 2*monopoleStrength+1)
{
    if(wfData->type==FQHE::_LAUGHLIN_)
    {
	    //  Initialize electron configuration - Laughlin

	    sphereCoords[0]=dcmplx(2.0*acos(pow((1.0/(2.0*wfData->monopoleStrength)),0.5)),Modulus(mt->random(), 2.0*PI)) ;
	
	    for (int i=1;i<wfData->nbr-1;i++)
	    {
		    sphereCoords[i]=dcmplx(2.0*acos(pow(((i*wfData->jastrowExponent)/(2.0*wfData->monopoleStrength)),0.5)),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    sphereCoords[wfData->nbr-1]=dcmplx(2.0*acos(pow(((nbrOrbitals-2.0)/(2.0*wfData->monopoleStrength)),0.5)),Modulus(mt->random(), 2.0*PI)) ;

	    utilities::cout.DebuggingInfo()<<"_LAUGHLIN_"<<std::endl;

    }
    else if(wfData->type==FQHE::_COMPOSITE_FERMION_ && abs(cfData->LLup)==3 && cfData->NEF==false && wfData->statistics==FQHE::_BOSONS_)
    {
	    ///// 3/4 Jain state /////

	    //	The empirical coefficients to give the correct co-ordinate initialization
	    double initCoeffs[36]={0.1,1.0,2.0,3.0,4.0,5.0,7.0,8.0,8.8,10.0,11.0,11.8,14.0,15.0,15.8,17.0,18.0,18.8,
			    21.0,22.0,22.8,24.0,25.0,25.8,28.0,29.0,29.8,31.0,32.0,32.8,35.0,36.0,36.8,38.0,39.0,39.8};

	    for (int i=0;i<wfData->nbr;i++)
	    {
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    utilities::cout.DebuggingInfo()<<"JAIN 3/4"<<std::endl;

    }
    else if(wfData->type==FQHE::_BONDERSON_SLINGERLAND_ && abs(cfData->LLup)==2  && cfData->NEF==true && wfData->nbr<=36)
    {
	     ///// 2/5 BS state /////

        //  FIXME

	    //	The empirical coefficients to give the correct co-ordinate initialization 
	    double initCoeffs[36]={0.1, 1.0, 4.0, 5.0, 8.0, 9.0, 12.0, 13.0, 13.8, 16.0, 16.8, 17.0,20.0,20.8,22.8,
							   23.0, 28.0, 29.0, 29.8,32.0,32.8,33.0,55.8,60.0,60.8,65.0,65.8,70.0,70.8,75.0,75.8,80.0,80.8,85.0,85.8};
	
	    //double initCoeffs[36]={0.1, 1.0, 5.0, 6.0, 10.0, 11.0, 15.0, 16.0, 20.0, 21.0, 25.0, 25.8,30.0, 31.0, 35.0, 35.8,
		//					   40.0, 40.8, 45.0, 45.8,50.0,50.8,55.0,55.8,60.0,60.8,65.0,65.8,70.0,70.8,75.0,75.8,80.0,80.8,85.0,85.8};
	
	    for (int i=0;i<wfData->nbr;i++)
	    {
	        std::cout<<"initCoeffs[i] = "<<initCoeffs[i]<<" arg = "<<((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength)<<std::endl;
	    
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    utilities::cout.DebuggingInfo()<<"BS 2/5"<<std::endl;
    }
    else if(wfData->type==FQHE::_COMPOSITE_FERMION_ && abs(cfData->LLup)==2 && cfData->NEF==false && wfData->nbr<=50)
    {
	    //	The empirical coefficients to give the correct co-ordinate initialization 
	
        ///// 2/5 Jain state /////

	    double initCoeffs[50]={0.1, 1.0, 5.0, 5.8, 10.0, 10.8, 15.0, 15.8, 20.0,                  
        20.8,25.0,25.8,30.0,30.8,35.0,35.8,40.0,40.8,45.0,45.8,50.0,50.8,55.0,55.8,60.0,60.8,65.0,65.8,70.0,70.8,75.0,
        75.8,80.0,80.8,85.0,85.8,90.0,90.8,95.0,95.8,100.0,100.8,105.0,105.8,110.0,110.8,115.0,115.8,120.0,120.8};
	
	    for (int i=0;i<wfData->nbr;i++)
	    {
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    utilities::cout.DebuggingInfo()<<"JAIN 2/5"<<std::endl;
    }
    else if(wfData->type==FQHE::_COMPOSITE_FERMION_ && abs(cfData->LLup)==3 && cfData->NEF==false && wfData->nbr<=51)
    {
	    //	The empirical coefficients to give the correct co-ordinate initialization 
	
         ///// 3/7 Jain state /////
	
	    double initCoeffs[51]={0.1, 1.0, 2.0, 7.0, 8.0, 8.8, 14.0, 15.0, 15.8,                  
        21.0,22.0,22.8,28.0,29.0,29.8,35.0,36.0,36.8,42.0,43.0,43.8,49.0,50.0,50.8,56.0,57.0,57.8,63.0,64.0,64.8,70.0,
        71.0,71.8,77.0,78.0,78.8,84.0,85.0,85.8,91.0,92.0,92.8,98.0,99.0,99.8,105.0,106.0,106.8,112.0,113.0,113.8};
	
	    for (int i=0;i<wfData->nbr;i++)
	    {
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    utilities::cout.DebuggingInfo()<<"JAIN 3/7"<<std::endl;
    }
    else if(wfData->type==FQHE::_COMPOSITE_FERMION_ && abs(cfData->LLup)==2 && cfData->NEF==true && wfData->nbr<=40)
    {
	    //	The empirical coefficients to give the correct co-ordinate initialization 
	
	    //////  Jain 2/3  //////
	
	    double initCoeffs[40]={0.1, 1.0, 3.0, 3.8, 6.0, 6.8, 9.0, 9.8, 12.0, 
	    12.8,15.0,15.8,18.0,18.8,21.0,21.8,24.0,24.8,27.0,27.8,30.0,30.8,33.0,33.8,36.0,36.8,39.0,39.8,42.0,42.8,45.0,45.8,48.0,48.8,51.0,51.8,53.0,53.8,56.0,56.8};
	
	    for (int i=0;i<wfData->nbr;i++)
	    {
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    utilities::cout.DebuggingInfo()<<"JAIN 2/3"<<std::endl;
    }
    else if(wfData->type==FQHE::_COMPOSITE_FERMION_ && abs(cfData->LLup)==3 && cfData->NEF==true && wfData->nbr<=30)
    {
	    //	The empirical coefficients to give the correct co-ordinate initialization 

	    //////  Jain 3/5  //////

	    double initCoeffs[30]={0.1, 1.0, 2.0, 5.0, 6.0, 6.8, 10.0, 11.0, 11.8, 15.0, 16.0, 16.8, 20.0, 21.0, 21.7,
	    25.0,26.0,26.8,30.0,31.0,31.8,35.0,36.0,36.8,40.0,41.0,41.8,45.0,46.0,46.8};

	    for (int i=0;i<wfData->nbr;i++)
	    {
		    sphereCoords[i]=dcmplx(acos((initCoeffs[i]-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }
	
	    utilities::cout.DebuggingInfo()<<"JAIN 3/5"<<std::endl;
    }
    else
    {
	    //	read values from a file called initCoeffs.dat
	
	    std::ifstream f_init;
	
	    f_init.open("initCoeffs.dat");
	
	    if(f_init.is_open()==0)
	    {
		    std::cerr<<"\tERROR: attempted to read file initCoeffs.dat but file was not found"<<std::endl<<std::endl;
		    exit(EXIT_FAILURE);
	    }
	
	    double initCoeff;

	    for (int i=0;i<wfData->nbr;i++)
	    {
		    f_init>>initCoeff;
		    sphereCoords[i]=dcmplx(acos((initCoeff-wfData->monopoleStrength)/wfData->monopoleStrength),Modulus(mt->random(), 2.0*PI)) ;
	    }

	    f_init.close();
    }
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/


////////////////////////////////////////////////////////////////////////////////
//! \brief This function generates a set of Metropolis moves
//!
//! \return Number of accepted Metropolis moves
//!
////////////////////////////////////////////////////////////////////////////////

int DoMetropolisSampling(
    const GeneralOptions* parameters, //!<    Pointer to parameter set
    MersenneTwister* mt,              //!<    Pointer to random number generator object
    const MetropOptions* metropolis,  //!<    Pointer to metropolis parameters
    FQHE::WaveFunction* wf,           //!<    Pointer to wave function object
    dcmplx *sphereCoords,   	      //!<    An array of spherical coordinates stored as theta + I phi
    const int nbrOrbitals)            //!<    Number of orbitals (defined as 2*monopoleStrength+1)
{
    const int nbr = wf->GetNbrParticles();
    const double monopoleStrength = wf->GetMonopoleStrength();

	//	Allocate memory

    //	new sample set of coordinates (theta and phi stored)
	dcmplx* new_sphereCoords = new dcmplx[nbr];
	
	//	spinor coordinates
	dcmplx* u = new dcmplx[nbr];
	dcmplx* v = new dcmplx[nbr];
   
	//	Convert "sphereCoords" (theta+I*phi) to spinor co-ordinates
	
	for(int i=0;i<nbr;i++)
	{
		u[i]= cos(real(sphereCoords[i]/2.0))*exp(dcmplx(0,imag(sphereCoords[i])/2.0));
		v[i]= sin(real(sphereCoords[i]/2.0))*exp(dcmplx(0,-imag(sphereCoords[i])/2.0));
	}

	//	calculate wave function value

	dcmplx oldwfValue = wf->EvaluateWfSphere(nbr,u,v);

	//	determine the disk area covered by the wave function

	double oldDiskArea = 1.0;

	for(int i=0; i<nbr; i++)
	{
		oldDiskArea *= sin(real(sphereCoords[i]));
	}

    int accepted=0;
    int rejected=0;

	//	do some number metropolis steps
	for(int k=0;k<metropolis->iter;k++)
	{
	    //  Local variables
	
	    int t1,t2,t3,t4;			        //	labels of particles to be moved
	    int	lA,lB;					        //	angular momentum of a given A or B particle
	    int	nbrA,nbrB;				        //	number of A/B subsystem particles in the sample
	    int lzA_sample;				        //	total angular momentum of the sample

		for(int i=0;i<nbr;i++)
		{
			new_sphereCoords[i]=sphereCoords[i];
		}

		//	generate four non-equal integers t1,t2,t3,t4

		t1=mt->genrand_int31()%nbr;

		do
		{
			t2=mt->genrand_int31()%nbr;
			t3=mt->genrand_int31()%nbr;
			t4=mt->genrand_int31()%nbr;
		}
		while(t2==t1 || t2==t3 || t2==t4|| t3==t1|| t3==t4|| t4==t1);

		new_sphereCoords[t1] = dcmplx(Modulus(real(new_sphereCoords[t1]) + (metropolis->width)*(mt->random()-0.5),PI),
		Modulus(imag(new_sphereCoords[t1]) + (metropolis->width)*(mt->random()-0.5), 2.0*PI));
		
		new_sphereCoords[t2] = dcmplx(Modulus(real(new_sphereCoords[t2]) + (metropolis->width)*(mt->random()-0.5),PI),
		Modulus(imag(new_sphereCoords[t2]) + (metropolis->width)*(mt->random()-0.5), 2.0*PI));
		
		new_sphereCoords[t3] = dcmplx(Modulus(real(new_sphereCoords[t3]) + (metropolis->width)*(mt->random()-0.5),PI),
		Modulus(imag(new_sphereCoords[t3]) + (metropolis->width)*(mt->random()-0.5), 2.0*PI));
		
		new_sphereCoords[t4] = dcmplx(Modulus(real(new_sphereCoords[t4]) + (metropolis->width)*(mt->random()-0.5),PI),
		Modulus(imag(new_sphereCoords[t4]) + (metropolis->width)*(mt->random()-0.5), 2.0*PI));

		nbrA=0;
		nbrB=0;
		lzA_sample=0;

		for(int i=0;i<parameters->nbrFixed;i++)
		{
			//	determine l_i_A as an integer nearest to |z_i|^2
			lA = (int) round(monopoleStrength*(1.0 + cos(real(new_sphereCoords[i]))) );

			lzA_sample += lA;

			//	count the number of particles in the A subsystem
			if (lA>=0 && lA <= parameters->realCut)	nbrA++;
		}

		for(int i=parameters->nbrFixed;i<nbr;i++)
		{
			//	determine l_i_B as an integer nearest to |z_i|^2
			lB = (int) round(monopoleStrength*(1.0 + cos(real(new_sphereCoords[i]))) );

			//	count the number of particles in the B subsystem
			if (lB>=parameters->realCut && lB <=nbrOrbitals)	nbrB++;
		}

		//	Only accept the move if the angular momentum is less than the selected lzA
		//	and the A/B partition is equal to that selected by nbrFixed

		if(lzA_sample<=parameters->lzA && nbrA==parameters->nbrFixed && nbrB==(int) (nbr-parameters->nbrFixed))
		{

			//	Convert new_sphereCoords to spinor co-ordinates
			for(int i=0;i<nbr;i++)
			{
				u[i]= cos(real(new_sphereCoords[i]/2.0))*exp(dcmplx(0,imag(new_sphereCoords[i])/2.0));
				v[i]= sin(real(new_sphereCoords[i]/2.0))*exp(dcmplx(0,-imag(new_sphereCoords[i])/2.0));
				//std::cout<<u[i]<<std::endl;
			}//getchar();

			//	calculate wave function value

			dcmplx newwfValue = wf->EvaluateWfSphere(nbr,u,v);

			//	determine the disk area covered by the wave function

			double newDiskArea = 1.0;
			for(int i=0; i<nbr; i++)
			{
				newDiskArea *= sin(real(new_sphereCoords[i]));
			}

			//	weight the metropolis. sampling test by the abs of wave function value
			//	and by the disk area (N.B NOT the abs squared of wave function)
			
			double ratio=abs(exp(newwfValue-oldwfValue))*(newDiskArea/oldDiskArea);

			if((ratio>=1.0 ) || ratio>mt->random())
			{
				accepted++;
				//	update old values if the move is accepted
				sphereCoords[t1] = new_sphereCoords[t1];
				sphereCoords[t2] = new_sphereCoords[t2];
				sphereCoords[t3] = new_sphereCoords[t3];
				sphereCoords[t4] = new_sphereCoords[t4];
				oldwfValue = newwfValue;
				oldDiskArea = newDiskArea;
			}
			else
			{
				rejected++;
			}
		}
		else
		{
			rejected++;
		}
	}

	utilities::cout.SecondaryOutput()<<"\tAccepted "<<accepted<<"\tRejected "<<rejected<<"\tTypical log wf value "<<parameters->renormFactor*log(10)+oldwfValue<<std::endl;

	//	deallocate memory

	delete[]	new_sphereCoords;
	delete[]	u;
	delete[]	v;

	//	return number of accepted moves
	
	return accepted;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief This function sorts each separate sample such that the co-ordinates
//! are in order of their lzA value
//!
////////////////////////////////////////////////////////////////////////////////

void OrderIndexes(
    const int nbrSamples,  //!<	Number of coordinate samples  
    const int nbr,         //!<	Number of coordinate in each sample
    const int nbrFixed,    //!<	Number of coordinates no in region A
    int* listLzA,          //!<	List of the lzA values of each co-ordinate in the 
						   //!	sample (of dimension samples*nbr)
    dcmplx* listCoord)     //!<	List of co-ordinate samples (of dimension samples*nbr)
{
    for (int i=0;i<nbrSamples;i++)
    {
	    for(int j=0;j<nbrFixed-1;j++)
	    {
		    for(int k=j+1;k<nbrFixed;k++)

		    {
			    if (listLzA[i*nbr+j]<=listLzA[i*nbr+k])
			    {
				    int swap_lzA=listLzA[i*nbr+j];

				    listLzA[i*nbr+j]=listLzA[i*nbr+k];
				    listLzA[i*nbr+k]=swap_lzA;

				    dcmplx swap_coord=listCoord[i*nbr+j];

				    listCoord[i*nbr+j]=listCoord[i*nbr+k];
				    listCoord[i*nbr+k]=swap_coord;
			    }
		    }
	    }
    }
	
	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

////////////////////////////////////////////////////////////////////////////////
//! \brief This function filters a set of coordinate samples by applying an
//! additional constraint to select samples with high overlap with the low
//! lying states in the real space entanglement spectrum
//!
//! \return The dimension of the filtered list
//!
////////////////////////////////////////////////////////////////////////////////

int ApplyLowLyingConstraint(
    const GeneralOptions* parameters, //!<    Pointer to parameter set
    const int nbrSamples,       //!<	Number of coordinate samples  
    const int nbr,              //!<	Number of coordinate in each sample
    const int nbrFixed,         //!<	Number of coordinates no in region A
    int* listLzA,         		//!<	List of the lzA values of each co-ordinate 
								//!	    in the sample (of dimension samples*nbr)
    dcmplx* listCoord,    		//!<	List of co-ordinate samples (of dimension samples*nbr)
	const int nbrOrbitals)      //!<	Number of orbitals (2*monopoleStrength+1)
{
    utilities::cout.MainOutput()<<"\tApplying constraint..."<<std::endl;

    int newDim=0;
    int counter;
			
	for(int i=0;i<nbrSamples;i++)
	{
		counter=0;
		
		for (int j=nbrFixed;j<nbr;j++)
		{ 
				//if (listLzA[i*nbr+j] >=  parameters->realCut && listLzA[i*nbr+j] <=  nbrOrbitals)  ;
				counter++;
		}
		
		if(counter==nbr-nbrFixed)
		{
			for(int k=0;k<nbr;k++)
			{ 
				listLzA[newDim*nbr+k]=listLzA[i*nbr+k];
				listCoord[newDim*nbr+k]=listCoord[i*nbr+k];
			}
			newDim++;
		}
	}  

    return newDim;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

////////////////////////////////////////////////////////////////////////////////
//! \brief This function filters a set of coordinate samples by applying an
//! additional constraint to select samples with no particles in the A region.
//!
//! \return The dimension of the filtered list
//!
////////////////////////////////////////////////////////////////////////////////

int FilterEmptyRegionA(
    const GeneralOptions* parameters, //!<    Pointer to parameter set
    const int nbrSamples, //!<	Number of coordinate samples  
    const int nbr,        //!<	Number of coordinate in each sample
    const int nbrFixed,   //!<	Number of coordinates no in region A
    int* listLzA,         //!<	List of the lzA values of each co-ordinate
						  //!	in the sample (of dimension samples*nbr)
    dcmplx* listCoord)    //!<	List of co-ordinate samples (of dimension samples*nbr)
{
    int newDim=0;
    int counter1,counter2;	
			
	for(int i=0;i<nbrSamples-1;i++)
	{
		counter1=0;
		for(int k=0;k<nbrFixed;k++)
		{
			if (listLzA[i*nbr+k]==0)	counter1++;
		}
		
		for(int j=i+1;j<nbrSamples;j++)
		{
			counter2=0; 
			for (int k=0;k<parameters->nbrFixed;k++)
			{
				if(listLzA[i*nbr+k]==listLzA[j*nbr+k] && counter1!=nbrFixed)	counter2++;
			}
			
			if(counter2==nbrFixed)
			{
				for (int k=0;k<nbrFixed;k++)
				{
					listLzA[j*nbr+k]=0;
					listCoord[j*nbr+k]=0;
				}
			}
		}
		
		if(counter1 != nbrFixed)
		{
			for (int k=0;k<nbr;k++)
			{ 
				listLzA[newDim*nbr+k]=listLzA[i*nbr+k];
				listCoord[newDim*nbr+k]=listCoord[i*nbr+k];
			}
			newDim++;
		}
	} 
	
	return newDim;	
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

////////////////////////////////////////////////////////////////////////////////
//! \brief This function evaluates a fast Fourier transform
//!
////////////////////////////////////////////////////////////////////////////////

void EvaluateFourierTransform(
    const GeneralOptions* parameters, //!<    Pointer to parameter set
    FQHE::WaveFunction* wf,     //!<    Pointer to wave function object
    const dcmplx *sphereCoords, //!<   An array of spherical coordinates stored as theta + I phi
    dcmplx *waveFunctionValue)  //!<   An array to store a set of wave function values for each sector
{
	//////////////////////////////////////////////
	//	Allocate memory to store spinor co-ordinate arrays
	//  Fourier component list and FFTW input/output arrays

    const int nbr = wf->GetNbrParticles(); 

	dcmplx* u = new dcmplx[nbr];
	dcmplx* v = new dcmplx[nbr];
	dcmplx* fourierComp = new dcmplx[parameters->nbrFourier];

	fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * parameters->nbrFourier);
	fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * parameters->nbrFourier);

	//////////////////////////////////////////////
	//	Declare FFTW "plans"
	
	//   fftw_plan forward=fftw_plan_dft_1d(parameters->nbrFourier,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_plan  backward=(fftw_plan)fftw_plan_dft_1d(parameters->nbrFourier,in,out,FFTW_BACKWARD,FFTW_PATIENT);
	
	//	convert sphereCoords to spinor v and u values
	for(int i=parameters->nbrFixed;i<nbr;i++)
	{
		u[i]= cos(real(sphereCoords[i]/2.0))*exp(dcmplx(0,imag(sphereCoords[i])/2.0));
		v[i]= sin(real(sphereCoords[i]/2.0))*exp(dcmplx(0,-imag(sphereCoords[i])/2.0));
	}
	
	//	run for each Fourier component
	for(int l=0;l<parameters->nbrFourier;l++)
	{ 
		
		for(int i=0;i<parameters->nbrFixed;i++)
		{
			//	implement angular momentum constraint in the Fourier transform
			u[i]= cos(real(sphereCoords[i]/2.0))*exp(dcmplx(0,imag(sphereCoords[i])/2.0))*exp(dcmplx(0,-2.0*l*PI/parameters->nbrFourier)) ;
			v[i]= sin(real(sphereCoords[i]/2.0))*exp(dcmplx(0,-imag(sphereCoords[i])/2.0));//*exp(dcmplx(0,l*PI/parameters->nbrFourier)) ;
			
			//std::cout<<sphereCoords[i]<<"\t"<<u[i]<<" "<<v[i]<<std::endl;
			
			if(std::isnan(real(sphereCoords[i])) || std::isnan(imag(sphereCoords[i])))
			{
				std::cerr<<"\tWARNING: nan wave z[i]"<<std::endl;
				getchar();
			}
		}

		//	Evaluate wave function
		
		dcmplx wfValue;  //	wave function value
		
		wfValue=wf->EvaluateWfSphere(nbr, u, v);	//	calculate the log of the wave function

		wfValue+=parameters->renormFactor*log(10);		//	renormalize the value by 10^(renormFactor*nbr) to avoid overflows

		//std::cout<<"NODE "<<mpi.m_id<<" renorm factor: "<<parameters->renormFactor<<std::endl;getchar();
		//std::cout<<u[0]<<" "<<v[0]<<" "<<wfValue<<std::endl;getchar();
		
		wfValue=exp(wfValue);			//	exponentiate
		
		if(std::isnan(real(wfValue)) || std::isinf(real(wfValue)) || real(wfValue)==0)
		{
			std::cerr<<"\tWARNING: 0 nan or inf wave function value detected, renorm is "<<parameters->renormFactor<<std::endl;	
			getchar();
		}

	    //std::cout<<u[0]<<" "<<v[0]<<" "<<wfValue<<std::endl;
	
		fourierComp[l] = wfValue*exp(dcmplx(0,(l*2.0*PI/parameters->nbrFourier)*parameters->lzA)) ;

	}

	for(int i=0;i<parameters->nbrFourier;i++)
	{
		in[i][0] = real(fourierComp[i]);
		in[i][1] = imag(fourierComp[i]);
	}

	//  fftw_execute(forward);
	fftw_execute(backward);

	for (int i=0;i<parameters->nbrSect;i++)
	{        
		waveFunctionValue[i]=(1.0/parameters->nbrFourier)*dcmplx(out[i][0],out[i][1]);
	}

	//	deallocate memory

	delete[]	u;
	delete[]	v;
	delete[]	fourierComp;
	
	fftw_destroy_plan(backward);
	fftw_free(in); 
	fftw_free(out);

	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/

////////////////////////////////////////////////////////////////////////////////
//! \brief This function generates and executes a python script to make
//! a rough plot of the output spectrum
//!
////////////////////////////////////////////////////////////////////////////////

void MakeRoughPlot(
    const GeneralOptions* parameters,       //!<    Pointer to parameter set
    const FQHE::WaveFunctionData* wfData,   //!<    Pointer to wave function data
    const FQHE::CompositeFermionData* cfData,
                                        //!<    Pointer to additional CF wave function data
    const utilities::MpiWrapper& mpi)
{
    if(0 == mpi.m_id)	// FOR THE MASTER NODE
    {
        utilities::cout.MainOutput()<<"\n"<<utilities::cout.HyphenLine()<<"\n"<<utilities::cout.HyphenLine()<<std::endl;

        utilities::cout.MainOutput()<<"\n\tGenerating rough plot"<<std::endl;
    
        //  Generate a python script to read in and plot the entanglement spectrum
    
        std::stringstream pythonScript;
        
        pythonScript.str();
        
        pythonScript<<"#! //usr//bin//env python"<<_PYTHON_VERSION_<<"\n"
        "import matplotlib                              \n"
        "import numpy as np                             \n"
        "import math                                    \n"
        "matplotlib.use('Agg')                        \n\n"
        "import matplotlib.pyplot as plt              \n\n"
        "eigenvalues = []                               \n"
        "sectors = []                                   \n"
        "sector = "<<parameters->lzA<<"                \n\n";
        
        for(int k=0;k<parameters->nbrSect;k++)
        {
            pythonScript<<"fin=open('"<<GenerateEigenvaluesFileName(k,parameters,cfData,wfData)<<"')\n\n"
            "for line in fin:                           \n"
            "   eigenvalue = -float(line)               \n"
            "   eigenvalues.append(eigenvalue)          \n"
            "   sectors.append(sector)                \n\n"
            "fin.close()                                \n"
            "sector += 1                              \n\n";
        }
        
        pythonScript<<"main_axes = plt.axes([0.1,0.1,0.75,0.75])\n"
        "main_axes.get_xaxis().tick_bottom()            \n"
        "main_axes.get_yaxis().tick_left()              \n"
        "main_axes.spines['right'].set_visible(False)   \n"
        "main_axes.spines['top'].set_visible(False)     \n"
        "main_axes.set_ylabel(r'$\\xi$',rotation='horizontal')\n"
        "main_axes.yaxis.set_label_coords(0.02,1.01)    \n"
        "plt.xticks(range(0,sector+1),rotation='vertical',fontsize=4)\n"
        "main_axes.scatter(sectors,eigenvalues,s=40,c='black',marker='_',lw=0.4)\n"
        "main_axes.get_figure().savefig('"<<GeneratePlotFileName(parameters,cfData,wfData)<<"',bbox_inches='tight')\n";
        
        //  Execute the script
    
        utilities::Script myScript;

        myScript.SetScript(pythonScript.str());
        
        myScript.Execute();
        
        utilities::cout.MainOutput()<<"\n\tDONE! See it at "<<GeneratePlotFileName(parameters,cfData,wfData)<<std::endl;
    }
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
