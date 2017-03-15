////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!		This is the header file for the wave function base class
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

#ifndef _WAVEFUNCTION_HPP_INCLUDED_
#define _WAVEFUNCTION_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <utilities/general/dcmplx_type_def.hpp>
#include <utilities/general/pi_const_def.hpp>
#include <utilities/general/i_const_def.hpp>
#include <utilities/general/cout_tools.hpp>
#if _ENABLE_MPI_
#include <utilities/wrappers/mpi_wrapper.hpp>
#endif
#if _DEBUG_
#include <utilities/general/debug.hpp>
#endif

//!
//!	A namespace containing all functions and classes used in 
//!	fractional quantum Hall effect calculations. 
//!
namespace FQHE
{
    #ifndef _GEOMETRY_TYPE_DEFINED_
    #define _GEOMETRY_TYPE_DEFINED_
    //!
    //!	Define a variable specifying the possible geometries
    //!
    enum geometry_t {_SPHERE_, _DISC_, _TORUS_};
    #endif
    #ifndef _WAVE_FUNCTION_TYPE_DEFINED_
    #define _WAVE_FUNCTION_TYPE_DEFINED_
    //!
    //!	Define a variable specifying the possible wave function types
    //!
    enum waveFunction_t {_LAUGHLIN_, _NASS_, _COMPOSITE_FERMION_, _MOORE_READ_, 
                         _BONDERSON_SLINGERLAND_, _PARAFERMION_};	
    #endif				
    #ifndef _STATISTICS_TYPE_DEFINED_
    #define _STATISTICS_TYPE_DEFINED_	
    //!
    //!	Define a variable specifying the possible statistics		
    //!			
    enum statistics_t {_BOSONS_, _FERMIONS_};	
    #endif

    //!
    //! Generate a list of composite fermion wave function program options
    //!
    inline boost::program_options::options_description GetWaveFunctionOptions()
    {
        boost::program_options::options_description waveFunctionOptions("Wave Function Options");
        waveFunctionOptions.add_options()
        ("type", boost::program_options::value<std::string>()->default_value("l"),
         "Type of wave function to use \n\nALLOWED TYPES:\n"
         "  l (Laughlin)\n  cf (composite fermion)\n"
         "  bs (Bondersan-Slingerland)\n"
         "  nass (Non-abelian spin-singlet)\n"
         "  para (Z-k parafermions)\n")
        ("geometry,g", boost::program_options::value<std::string>()->default_value("s"),
	     "Select geometry: d: Disc, s: Sphere or t: Torus\n")
        ("nbr,n", boost::program_options::value<int>()->default_value(16),
         "Number of particles\n")
        ("bosonic",boost::program_options::value<bool>()->default_value(false),
         "Set to divide the CF wavefunction by a Jastrow factor (z_i-z_j)\n")
        ("jastrow-exponent,p",boost::program_options::value<int>()->default_value(3),
	     "Function depends on wf type:\n  l type: \tSet Jastrow^p factor in the wavefunction\n"
	     "  NASS type: \tSet Jastrow^p factor in the wavefunction\n"
	     "  mr type: \tSet Jastrow^p factor in the wavefunction\n"
	     "  cf/bs type: \tSet Jastrow^2p factor in the wavefunction (number of flux attatched)\n"
	     "  para type:  \tSet Jastrow^M factor in the cluster group\n")
	     ("cluster-size,k",boost::program_options::value<int>()->default_value(3),
	     "Cluster size (k) in parafermion wave function.\n");
	    #if _ENABLE_LATTICE_FQHE_
	    waveFunctionOptions.add_options()
	    ("torus-x,x", boost::program_options::value<int>()->default_value(10),
	     "Select x dimension of rectangular torus\n")
	    ("torus-y,y", boost::program_options::value<int>()->default_value(10),
	     "Select y dimension of rectangular torus\n")
        ("torus-state",boost::program_options::value<int>()->default_value(0),
         "pick the index of the torus state (a value from 0 to torus degeneracy-1)\n");
        #endif
        return waveFunctionOptions;
    };

    //!
    //!	Define a data structure to contain a very general set of variables 
    //!	associated with a FQHE wave function.
    //!
    struct WaveFunctionData
    {
        int nbr;		        //!<	Total number of electrons
							    //!
	    int nbrUp;		        //!<	Number of spin up species
							    //!
	    int nbrDown;	        //!<	Number of spin down species
							    //!
	    int nbrQh;				//!<	Number of quasi-holes
							    //!
	    int jastrowExponent;	//!<	Function depends on wf type:
							    //!	    - l type: Set Jastrow^p factor
							    //!     - NASS type: Set Jastrow^p factor
							    //!     - mr type: Set Jastrow^p factor 
							    //!     - cf/bs type: Set Jastrow^2p factor
							    //!     - para type: Set Jastrow^M factor 
							    //!     in cluster group
	    int nbrClusters;        //!     k parameter in parafermion wave functions
        int fillNumerator;		//!<	Filling fraction numerator	
							    //!
        int fillDenominator;	//!<	Filling fraction denominator
							    //!
        double fillingFactor;	//!		Value of the filling fraction
							    //!
        double flux;			//!<	Value of the flux 
							    //!
        std::string wfName;     //!<  	String to store wave function name 
							    //!		(for printing data to screen)
        std::string wfFileName;	//!<	String to store file name identifier 
							    //!		That will be associated with any output data file
	    std::string	path;		//!<	Name of file path where data are stored
							    //!
	    int fusionChannels;     //!<  	Number of fusion channels for quasi-hole state
							    //!
	    geometry_t geometry;	//!<	Specify a geometry
							    //!
	    waveFunction_t type;	//!<	Specify a type of wave function
							    //!
	    statistics_t statistics;//<	Specify the statistics (bosons or fermions)
							    //!
	    double monopoleStrength;//!<	The magnitude of the monopole strength
							    //!		At the centre of the sphere
        double shift;			//!<	Flux shift in the sphere geometry
                                //!
        double sphereRadius;    //!<    Sphere radius in units of l_0
                                //!
	    #if _ENABLE_TORUS_GEOMETRY_
	    int torusState;			//!<	An index to specify which torus state we evaluate
							    //!
        int torusDegeneracy;	//!<	Ground state degeneracy on a torus
							    //!
	    #endif
	    #if _ENABLE_DISC_GEOMETRY_
	    double discRadius;		//!<	Disc radius in units of l_0
							    //!
	    #endif
	    #if _ENABLE_LATTICE_FQHE_
	    double latticeSpacing;	//!<	Lattice spacing for lattice type wave functions
							    //!
        int Lx;		            //!<	x dimension of the lattice
							    //!
	    int Ly;		            //!<	y dimension of the lattice
							    //!
        int maxOccupation;		//!<	Maximum allowed occupation of lattice sites
        #endif
        WaveFunctionData();
        void InitFromCommandLine(
        #if _ENABLE_MPI_
            boost::program_options::variables_map* options,
            const utilities::MpiWrapper& mpi);
        #else
            boost::program_options::variables_map* options);
        #endif
        void InitRadius();
        void GenerateFileName(); 
	    void CheckAndPrint() const;                       
        #if _ENABLE_MPI_                                
        void MpiSync(const int syncId, const utilities::MpiWrapper& mpi);
        #endif
    };

    //!
    //!	An abstract base class to store an instance of the basic wave 
    //!	data structure (WaveFunctionData) and virtual function definitions
    //!	used for the implementation of wave function algorithms in any derived
    //!	classes.
    //!
    class WaveFunction
    {
        protected:
        WaveFunctionData *m_wfData;

        public:
        WaveFunction(WaveFunctionData*);
	    virtual ~WaveFunction()=0;
	    //////		Variables used only in the sphere geometry		////////////////
	    #if _ENABLE_SPHERE_GEOMETRY_
	    virtual dcmplx EvaluateWfSphere(const int, dcmplx*,dcmplx*) const {return 0.0;};
	    virtual dcmplx EvaluateQuasiholeWfSphere(const int, dcmplx*, dcmplx*, const int, 
	                                             dcmplx*, dcmplx*) const {return 0.0;};
	    #endif
	    #if _ENABLE_DISC_GEOMETRY_
	    virtual dcmplx EvaluateWfDisc(const int, dcmplx*) const {return 0.0;};
	    #endif
	    #if _ENABLE_TORUS_GEOMETRY_
        virtual dcmplx EvaluateWfTorus(const int, std::complex<int>*, const int) const {return 0.0;};
	    virtual dcmplx EvaluateQuasiholeWfTorus(const int, std::complex<int>*, const int, 
	                                            const int, dcmplx*) const {return 0.0;};
	    #endif
	    #if _ENABLE_LATTICE_FQHE_
	    virtual void LatticeDataPrecalculation(const int, const int, dcmplx*){};
	    #endif
        int GetNbrParticles() const;
        double GetMonopoleStrength() const;
	    #if _DEBUG_
	    dcmplx m_prevWfValue;
	    dcmplx m_testValA;
	    dcmplx m_testValB;
	    #endif
    };
}   //  End FQHE namespace
#endif
