////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                           
//!                                                                             
//!                      \date Last Modified: 10/12/2014                        
//!                                                                             
//!	 \file
//!     This file implements an entanglement energy type fitting function
//!     as an alternative to the CFT form of the entanglement Hamiltonian           
//!                                                        
//!                    Copyright (C) 2014 Simon C Davenport
//!                                                                             
//!     This program is free software: you can redistribute it and/or modify
//!     it under the terms of the GNU General Public License as published by
//!     the Free Software Foundation, either version 3 of the License,
//!     or (at your option) any later version.
//!                                                                             
//!     This program is distributed in the hope that it will be useful, but
//!     WITHOUT ANY WARRANTY; without even the implied warranty of
//!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//!     General Public License for more details.
//!                                                                             
//!     You should have received a copy of the GNU General Public License
//!     along with this program. If not, see <http://www.gnu.org/licenses/>.
//!                                                                             
////////////////////////////////////////////////////////////////////////////////

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////

#include "entanglement_energy_model.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate a list of entanglement energies of the form
//!
//! H = sum_{m,sigma}  epsilon(m,sigma)  n_{m,sigma} + normalization.
//!
//! Where n_{m,sigma} is the occupation number of the m th angular momentum 
//! orbital and sigma is the Landau level index. 
//!
//! With epsilon(m,sigma) = Log[(1-alpha)/alpha] we propose:
//!
//! epsilon_{m,sigma} = a_0 + a_1 m + a_2 m^2 + a_3 m^3  + a_5 m^5
//! + b_0 sigma + b_1 sigma m + b_2 sigma m^2 + b_3 sigma m^3
//!
//! Where a_0,b_0 etc. are possible fitting parameters.
//!
//! The corresponding normalization is given by the sum of Log[1-alpha].
//!
//! This function calculates the values of epsilon(m,sigma) and the value
//! of the normalization
//!
////////////////////////////////////////////////////////////////////////////////

void EntanglementEnergyModel::GenerateEntanglementEnergies(
    double* energyLevels,           //!<    Set of diagonalized energy levels generated here
    const double q,                 //!<    Monopole strength
    double& normalization,          //!<    Value for the normalization to be set
    const int nbrLevels,            //!<    Number of Landau levels
    const int nbrStates,            //!<    Number of single particle orbitals)
    std::vector<double>& mp)   
                                    //!<    List of model parameters
    const
{
    //  Levels are labelled from most negative to most positive m value
    
    utilities::cout.SecondaryOutput()<<"\n\tGenerating energy levels:\n"<<std::endl;
    
    std::cout.precision(15);
    
    utilities::cout.SecondaryOutput()<<"\n\ta_0="<<mp[0]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\ta_1="<<mp[1]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\ta_2="<<mp[2]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\ta_3="<<mp[3]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\ta_5="<<mp[4]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\tb_0="<<mp[5]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\tb_1="<<mp[6]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\tb_2="<<mp[7]<<std::endl;
    utilities::cout.SecondaryOutput()<<"\n\tb_3="<<mp[8]<<std::endl<<std::endl;
    
    double z = -q-(nbrLevels-1);         //  z-component of angular momentum
    double* p_energy = energyLevels;
    
    normalization = 0.0;

    utilities::cout.DebuggingInfo()<<"Lz labels"<<std::endl;

    for(int m=0;m<nbrStates;m++,z+=1.0)
    {
        for(int l=0;l<nbrLevels;l++,p_energy++)
        {
            double sigma = 2.0*l-1;     //  MODIFICATION FOR THE 2 LEVEL CASE ONLY

            *(p_energy) = 
            mp[0] +
            (mp[1] + (mp[2] + (mp[3] + mp[4]*z*z)*z)*z)*z +
            sigma*(mp[5] + (mp[6] + (mp[7] + mp[8]*z)*z)*z);

            #if _ENABLE_HIGH_PRECISION_
            
            //  Calculate the normalization using high precision arithmetic
            
            const int P = 512;    //  Set precision
            
            mpfr_t temp;
            
            mpfr_init2(temp,P);
            
            mpfr_set_d(temp,-*(p_energy),MPFR_RNDN);
            
            mpfr_exp(temp,temp,MPFR_RNDN);
            
            mpfr_add_ui(temp,temp,1,MPFR_RNDN);
            
            mpfr_log(temp,temp,MPFR_RNDN);
            
            normalization -= mpfr_get_d(temp,MPFR_RNDN);
            
            mpfr_clear(temp);
            
            #else
            
            normalization -= log(exp(-*(p_energy))+1);
            
            #endif

            if(std::isinf(normalization) || std::isnan(normalization))
            {
                std::cerr<<"\n\tERROR: normalization of entanglement energy is nan or inf"<<std::endl;
                
                std::cerr<<"\t(Argument of log was "<<1.0-*(p_energy)<<")"<<std::endl;
                
                getchar();
                exit(EXIT_FAILURE);
            }
            
            if(std::isinf(*(p_energy)) || std::isnan(*(p_energy)))
            {
                std::cerr<<"\n\tERROR: entanglement energy is nan or inf"<<std::endl;

                getchar();
                exit(EXIT_FAILURE);
            }

        }
        
        utilities::cout.DebuggingInfo()<<"\t"<<z;
    }
    
    utilities::cout.DebuggingInfo()<<std::endl;
    
    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Function to set model parameters based on a set of command line arguments
//! describing a file name
//!
////////////////////////////////////////////////////////////////////////////////

void EntanglementEnergyModel::ParametersFromFile(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    std::vector<double>& mp,                            //!<    List of single particle energy 
                                                        //!     parameters (to be populated)
    std::vector<double>& qModelParameters)              //!<    List of q function parameters
                                                        //!     (to be populated)
    const
{
    std::stringstream fileName;
    
    //  Generate file name
    
    fileName.str("");

	fileName<<(*optionList)["in-path"].as<std::string>()<<"parameters_entanglement_energy_model_n_"
	        <<(*optionList)["nbr-n"].as<int>()<<"_na_"<<(*optionList)["nbr-a"].as<int>()<<"_ll_"<<(*optionList)["nbr-ll"].as<int>()<<".dat";
	        
    std::ifstream f_param;
    
    f_param.open(fileName.str().c_str(),std::ios::in);
	
	if(!f_param.is_open())
	{
		std::cerr<<"\n\tERROR File: "<<fileName.str()<<" not found !!! "<<std::endl<<std::endl;
		exit(EXIT_FAILURE);
	}

	//	skip the comment lines in the file

	std::string line;

	while(getline(f_param,line))
	{
		if(line[0]!='#')	break;
	}
	
	//  read in all remaining lines as separate parameters
	
	utilities::cout.DebuggingInfo()<<"\n\tReading in parameters from "<<fileName.str()<<std::endl;
	
	while(!f_param.eof())
    {
        double temp;
        
        f_param>>temp;
        
        utilities::cout.DebuggingInfo()<<"\t\t"<<temp<<std::endl;
        
        mp.push_back(temp);
    }
    
    f_param.close();
            
    if(mp.size()>23)
    {
        std::cerr<<"\n\tERROR: Too many parameters in file??"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if(mp.size()<22)
    {
        std::cerr<<"\n\tERROR: Too few parameters in file??"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    //  Put 9 parameters c_? and d_? into qModelParameters instead of mp
    
    for(int i=0;i<9;++i)
    { 
        qModelParameters.push_back(mp[13+i]);
    }
    
    //  Remove the last 10 entries
    
    mp.resize(13);

    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Write the energy spectrum to a file in an appropriate format given
//! the set of command line arguments
//!
////////////////////////////////////////////////////////////////////////////////

void EntanglementEnergyModel::SpectrumToFile(
    boost::program_options::variables_map* optionList,  
                                //!<    Parsed command line arguments
    double* spectrum,           //!<    Sorted list of entanglement eigenvalues
    std::vector<int>& spectrumLabels,
                                //!<    List of labels identifying the RSES branches
    const int spectrumDim,      //!<    Length of the spectrum array
    const int sector)           //!<    Angular momentum sector labelling the
                                //!     current set of data
    const
{
    std::stringstream fileName;
    
    //  Generate file name
    
    fileName.str("");

	fileName<<(*optionList)["out-path"].as<std::string>()<<"/eigenvalues_entanglement_energy_model"
	<<"_n_"
	<<(*optionList)["nbr-n"].as<int>()<<"_na_"<<(*optionList)["nbr-a"].as<int>()
	<<"_ll_"<<(*optionList)["nbr-ll"].as<int>()<<"_sector_"<<sector<<".dat";
	
	std::ofstream f_spectrum;
	
	f_spectrum.open(fileName.str().c_str(),std::ios::out);
	
	if(!f_spectrum.is_open())
	{
	    std::cerr<<"\n\tERROR: Could not create spectrum file: "<<fileName.str()<<std::endl;
	    exit(EXIT_FAILURE);
	}
	
	//  Output format is two columns: one for the eigenvalue, one for its label
	
	for(int i=0;i<spectrumDim;i++)
	{
	    f_spectrum<<spectrum[i]<<"\t"<<spectrumLabels[i]<<"\n";	
	}
	
	f_spectrum.close();

    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Apply a 1st order perturbation to the current energy levels.
//!
//! For non-degenerate levels there is a correction to each energy level
//! given by <occupation| H |occupation>
//!
//! For degenerate levels there is a correction given by diagonalizing in
//! the sub-space of degenerate levels, forming a matrix like
//!
//! <occupation'| H |occupation>
//!
////////////////////////////////////////////////////////////////////////////////
#if 0
void EntanglementEnergyModel::ApplyPerturbation(
    std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,   
                                    //!<    List of LL occupations
    double* energyLevels,           //!<    List of single particle energy levels
    //double* qLevels,              //!<    List of q levels
    double* spectrum,               //!<    Data structure to store the entanglement spectrum
    const int nbrStates,            //!<    Number of orbitals represented
    std::vector<double>& mp)   
                                    //!<    List of model parameters defining
                                    //!     the perturbation Hamiltonian
    const
{
    //////      GET DEGENERACIES        ////////////////////////////////////////

    //  For each state in the spectrum, bin the degeneracy count. For this
    //  process we can exploit that the spectrum is sorted in ascending
    //  order 
    
    std::vector<unsigned int> degeneracyBins;
    
    //  Allow degeneracyBins to be large enough to count each level as 
    //  non-degenerate
    
    degeneracyBins.reserve(occupationsList.size());
    
    {
        double* p_spectrum = spectrum;
        double* p_last = spectrum+occupationsList.size()-1;
        double degenTol = 0.000000000001;

        while(p_spectrum<p_last)
        {
            double* p_forward = p_spectrum+1;
            unsigned int degeneracyCount = 1;
        
            #if _DEBUG_
            
            PRINT("p_forward",*p_forward);
            PRINT("p_spectrum",*p_spectrum);
            
            #endif
            
            while(fabs((*p_forward-*p_spectrum))<degenTol)
            {
                degeneracyCount++;
                p_forward++;
            }
            
            p_spectrum = p_forward;
            
            degeneracyBins.push_back(degeneracyCount);
        }
        
        if(p_spectrum == p_last)
        {
            degeneracyBins.push_back(1);
        }
        
        if(degeneracyBins[0] > occupationsList.size())
        {
            degeneracyBins[0] = occupationsList.size();
        }
        
        #if _DEBUG_
        
        PRINT("degeneracyBins",&degeneracyBins[0],degeneracyBins.size());
        
        PAUSE();
            
        #endif
    }

    //degeneracyBins.shrink_to_fit();
    
    //////      APPLY FIRST ORDER PERTURBATION      ////////////////////////////

    {
        utilities::cout.SecondaryOutput()<<"\n\tV1 value = "<<mp[9]<<std::endl<<std::endl;
        utilities::cout.SecondaryOutput()<<"\n\tV3 value = "<<mp[10]<<std::endl<<std::endl;
        utilities::cout.SecondaryOutput()<<"\n\tV5 value = "<<mp[11]<<std::endl<<std::endl;
        utilities::cout.SecondaryOutput()<<"\n\tV7 value = "<<mp[12]<<std::endl<<std::endl;

        std::vector<double> pseudopotentials;
        std::vector<double> pseudopotentials2LL;

        //utilities::cout.SecondaryOutput()<<"\n\tVCoulomb value = "<<mp[9]<<std::endl<<std::endl;

        pseudopotentials = diagonalization::GenCoulombPseudopotentials(nbrStates,0,mp[9]);

        pseudopotentials[1] += mp[10];

        //pseudopotentials.push_back(0);
        //pseudopotentials.push_back(mp[9]);
        //pseudopotentials.push_back(0);
        //pseudopotentials.push_back(mp[10]);
        //pseudopotentials.push_back(0);
        //pseudopotentials.push_back(mp[11]);
        //pseudopotentials.push_back(0);
        //pseudopotentials.push_back(mp[12]);
        
        pseudopotentials2LL = pseudopotentials;
        
        //pseudopotentials2LL.push_back(0);
        //pseudopotentials2LL.push_back(mp[11]);
        //pseudopotentials2LL.push_back(0);
        //pseudopotentials2LL.push_back(mp[12]);

        int nbrLevels = occupationsList[0].size();
        
        if(1== nbrLevels)
        {
            diagonalization::SpherePseudopotentialHamiltonian hamiltonian(1,nbrStates,pseudopotentials);

            hamiltonian.BuildLookUpTables();

            hamiltonian.SetOccupationEnergies(energyLevels,nbrStates);
            
            //hamiltonian.SetDensityDensity(qLevels,nbrStates);
            
            double* p_spectrum = spectrum;
            
            int counter = 0;

            //  Generate a block diagonal Hamiltonian for each set of degenerate states
            
            for(auto& degen : degeneracyBins)
            {
                uint64_t* fockStateBuffer = new uint64_t[degen];

                //  Convert bitset types to fockState_t
                
                for(int i=0;i<degen;++i,++counter)
                {
                    fockStateBuffer[i] = occupationsList[counter][0].to_ulong();
                }
                
                hamiltonian.SetFockBasis(fockStateBuffer,degen); 

                hamiltonian.BuildHamiltonian();
                
                hamiltonian.SetNbrEigenvaluesToFind(degen);
                
                hamiltonian.Diagonalize();
                
                hamiltonian.GetEigenvalues(p_spectrum,degen);
                
                hamiltonian.ClearHamiltonian();
                
                p_spectrum += degen;
                
                delete[] fockStateBuffer;
            }
        }
        else if(2 == nbrLevels)
        {
            int totalNbrStates = 2*nbrStates-2;
        
            diagonalization::SphereTwoLevelPseudopotentialHamiltonian hamiltonian(1,totalNbrStates,pseudopotentials,pseudopotentials2LL);

            hamiltonian.BuildLookUpTables();

            //  Change the data storage format of the single particle energy
            //  levels for consistency with what SphereTwoLevelPseudopotentialHamiltonian
            //  expects

            double* newEnergies = new double[totalNbrStates];

            for(int i=0;i<nbrStates-2;++i)
            {
                newEnergies[i] = energyLevels[(i+1)*2];
            }
            
            for(int i=nbrStates-2;i<totalNbrStates;++i)
            {
                newEnergies[i] = energyLevels[(i-nbrStates+2)*2+1];
            }

            hamiltonian.SetOccupationEnergies(newEnergies,totalNbrStates);
            
            delete[] newEnergies;
            
            double* p_spectrum = spectrum;
            
            int counter = 0;
            
            for(auto& degen : degeneracyBins)
            {
                uint64_t* fockStateBuffer  = new uint64_t[degen];
                uint64_t* fockStateBuffer2 = new uint64_t[degen];

                //  Convert bitset types to fockState_t
                
                for(int i=0;i<degen;++i,++counter)
                {
                    fockStateBuffer[i] = occupationsList[counter][0].to_ulong();
                    fockStateBuffer2[i] = occupationsList[counter][1].to_ulong();
                    
                    //  For consistency, we need to right shift values in 
                    //  the LLL buffer by 1
                    
                    fockStateBuffer[i] = fockStateBuffer[i] >> 1;
                }

                hamiltonian.SetFockBasis(fockStateBuffer,fockStateBuffer2,degen); 

                hamiltonian.BuildHamiltonian();
                
                hamiltonian.SetNbrEigenvaluesToFind(degen);
                
                hamiltonian.Diagonalize();
                
                hamiltonian.GetEigenvalues(p_spectrum,degen);
                
                hamiltonian.ClearHamiltonian();
                
                p_spectrum += degen;
                
                delete[] fockStateBuffer;
                delete[] fockStateBuffer2;
                
                //getchar();
            }
            
        }
    }
    
    return;
}
#endif

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//! \brief Public interface to the model construction. This function sequentially
//! calls sub-functions to generate the appropriate spectrum
//!
////////////////////////////////////////////////////////////////////////////////

void EntanglementEnergyModel::BuildEntanglementEnergyModel(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzaMin,                                   //!<    Minimum sector calculated
    const int lzaMax)                                   //!<    Maximum sector calculated
    const
{
    utilities::cout.MainOutput()<<"\n\t=========================    CALCUALTING ENTANGLEMENT ENERGY RSES    ========================="<<std::endl;
    
    //  First get model parameters
    
    std::vector<double> mp;
    std::vector<double> qModelParameters;
    
    this->ParametersFromFile(optionList,mp,qModelParameters);
    
    //  Generate the set of energy levels according to the model description
    
    //  Import command line variables
    
    int nbrParticles = (*optionList)["nbr-n"].as<int>();
    int nbrA = (*optionList)["nbr-a"].as<int>();
    
    if(nbrParticles>_MAX_NUMBER_)
    {
        std::cerr<<"ERROR: entanglement energy model not calculable for more than "<<_MAX_NUMBER_<<" particles - the code needs to be modified."<<std::endl;
        return;
    }
    
    int nbrLevels = (*optionList)["nbr-ll"].as<int>();
    double monopoleStrength = (double)(nbrParticles-nbrLevels*nbrLevels)/(2*nbrLevels);
    //  Dimension of the list of pseudoenergies from -q-l to q+l with l=nbrLevels-1
    int nbrStates = (int)(1+2*(monopoleStrength+nbrLevels-1));

    utilities::cout.MainOutput()<<"\n\t\tNBR LEVELS\t\t"<<nbrLevels<<std::endl;
    utilities::cout.MainOutput()<<"\n\t\tMONOPOLE STRENGTH (=FLUX/2)\t\t"<<monopoleStrength<<std::endl;
    
    double* energyLevels  = new double[nbrStates*nbrLevels];
    double* qLevels       = new double[nbrStates*nbrLevels];
    double normalization  = 0.0;
    double qNormalization = 0.0;

    this->GenerateEntanglementEnergies(energyLevels,monopoleStrength,normalization,nbrLevels,nbrStates,mp);
    
    this->GenerateEntanglementEnergies(qLevels,monopoleStrength,qNormalization,nbrLevels,nbrStates,qModelParameters);

    utilities::cout.AdditionalInfo()<<"\tnormalization = "<<normalization<<std::endl<<std::endl;

    for(int i=0;i<nbrStates*nbrLevels;i++)
    {
        utilities::cout.SecondaryOutput()<<"\tSingle particle energy level ["<<i<<"] = "<<energyLevels[i]<<std::endl;
    }
    
    for(int i=0;i<nbrStates*nbrLevels;i++)
    {
        utilities::cout.SecondaryOutput()<<"\tSingle particle q level ["<<i<<"] = "<<qLevels[i]<<std::endl;
    }

    //  Represent occupation patterns of single particle orbital populations by 
    //  a fixed vector of N-bit bitsets. Each bitset represents the occupation
    //  of single particle orbitals within one Landau level. The vector of these
    //  containers is a list of all possible occupations
    
    for(int s=lzaMin;s<=lzaMax;s++)
    {
        utilities::cout.SecondaryOutput()<<"\n\t\tSECTOR "<<s<<std::endl;

        std::vector<std::vector<std::bitset<_MAX_NUMBER_>>> occupationsList;
        std::vector<int> spectrumLabels;
        
        occupationBasis::GenerateLevelOccupations(occupationsList,spectrumLabels,
                                       nbrParticles,nbrA,nbrLevels,s);
        
        utilities::cout.DebuggingInfo()<<"\n\t\tOCCUPATION PATTERNS"<<std::endl;
        
        for(unsigned int i=0;i<occupationsList.size();i++)
        {
            utilities::cout.DebuggingInfo()<<utilities::cout.HyphenLine()<<std::endl;
            
            for(int j=0;j<nbrLevels;j++)
            {
                utilities::cout.DebuggingInfo()<<"\t\t"<<occupationsList[i][j]<<std::endl;
            }
        }

        //  Allocate memory to store the RSES

        double* spectrum = new double[occupationsList.size()];
        
        occupationBasis::AssignEnergies(occupationsList,energyLevels,spectrum,normalization,nbrStates,nbrLevels);
        
        occupationBasis::AssignSquaredEnergies(occupationsList,qLevels,spectrum,qNormalization,nbrStates,nbrLevels);

        for(unsigned int i=0;i<occupationsList.size();i++)
        {
            utilities::cout.AdditionalInfo()<<"\tspectrum["<<i<<"] = "<<spectrum[i]<<std::endl;
        }
        
        #if 0
        if((*optionList)["use-perturbation"].as<bool>())
        {
        
            //  What does this do?
            if(s==lzaMin)
            {
                //  Correct the energy levels for the normalization
                
                double* p_levels = energyLevels;
                
                for(int i=0;i<nbrStates*nbrLevels;++i,++p_levels)
                {
                    *p_levels -= (normalization+2*qNormalization)/(nbrA);
                }
            }

            //  Make sure that the spectrum is sorted in ascending order before proceeding

            double* spectrumCopy = new double[occupationsList.size()];
            memcpy(spectrumCopy,spectrum,occupationsList.size()*sizeof(double));
 
            utilities::QuickSort<double,int,_ASCENDING_ORDER_>(spectrum,spectrumLabels.data(),spectrumLabels.size());
            
            utilities::QuickSort<double,std::vector<std::bitset<_MAX_NUMBER_>>,_ASCENDING_ORDER_>(spectrumCopy,occupationsList.data(),occupationsList.size());

            delete[] spectrumCopy;

            //  Optionally apply a perturbation to the current energy levels

            this->ApplyPerturbation(occupationsList,energyLevels,spectrum,nbrStates,mp);
        }
        #endif
        
        //  Sort the spectrum to be in ascending order
     
        utilities::QuickSort<double,int,_ASCENDING_ORDER_>(spectrum,spectrumLabels.data(),spectrumLabels.size());
        
        for(unsigned int i=0;i<occupationsList.size();i++)
        {
            utilities::cout.AdditionalInfo()<<"\tspectrum["<<i<<"] = "<<spectrum[i]<<std::endl;
        }
 
        this->SpectrumToFile(optionList,spectrum,spectrumLabels,occupationsList.size(),s);
        
        delete[] spectrum;
        
        //getchar();
    }

    delete[] energyLevels;
}
    
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

