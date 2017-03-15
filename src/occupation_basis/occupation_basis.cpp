////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!   	This file contains functions for generating the occupation basis of 
//!     Fractional Quantum Hall Quasihole wave functions for the construction
//!     of the associated real space entanglement spectrum.
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

#include "occupation_basis.hpp"

//////////////////////////////////////////////////////////////////////////////////
//!	\brief This function determines the value a1 or a2 closest to the value b. 
//!
//!	\return b-a1 if abs(b-a1)<abs(b-a2), otherwise b-a2
//////////////////////////////////////////////////////////////////////////////////
int occupationBasis::FindNearestValue(
	const int b,	//!<	An integer
	const int a1,	//!<	Test value 1
	const int a2)	//!<	Test value 2
{
	return abs(b-a1)<fabs(b-a2) ? b-a1 : b-a2;
}

//!
//!	This function calculates the current LzA value associated with a
//!	given signature of occupations for CF LLs
//!
double occupationBasis::FindLzA(
	std::vector<int>& occupations,      //!<    A list of LL occupations       
	int nbrLowest)                      //!<    Total number of spaces in The lowest Landau level
{
	int nbrLLs = occupations.size();
	double LzTot = 0.0;
	for(int l=0; l<nbrLLs; ++l)
	{
		double LzMin = (double)(-nbrLowest+1.0)/2.0 - l;	//	Gives the lowest LzA state 
		//	This calculation assumes that only the lowest LzA states are occupied
		//	i.e there are no excitations
		double LzA = LzMin;
		for(int o=0; o<occupations[l]; ++o, LzA+=1.0)
		{
			LzTot += LzA;
		}
	}
	return LzTot;
}

//!
//!	Generate all the possible assignments of the particle cut to the 
//!	set of LLs for a given sector. Also store their angular momentum relative 
//! to the lowest LzA cut and the allowed 'excitations' for finite-sized
//!	corrected counting.		
//!	
void occupationBasis::GenerateOccupationData(
	const int nbrParticles,	                            //!<    Number of particles in the full system
	const int nbrA,       		                        //!<    Number of particles in the cut	 
	const int lzA,                                      //!<    The angular momentum value
	const int nbrColours, 		                        //!<    The number of U(1) currents in the Hilbert space
	std::vector<std::vector<int> >& deltaNa0List,       //!<    List of minimum (Na-Na0) [to be populated]	
	std::vector<std::vector<int> >& nbrEachLL,          //!<    List of number in each LL [to be populated]
	std::vector<int>& deltaLza,	                        //!<    Difference in Lz compared to Na0 cut,
                                                        //!     [to be populated]
	std::vector<std::vector<int> >& maxExcitationSize)  //!<    Maximum allowed excitation in a given LL, 
                                                        //!	    [to be populated]
{
	switch(nbrColours)
	{
		case 1:		//	Single current case e.g. Laughlin, nu=1 IQHE
		{
			std::vector<int> deltaNa0(nbrColours);
			std::vector<int> occupations(nbrColours);
			std::vector<int> allowedExcitationSize(nbrColours);
			//	Set Na-Na0 to have the minimum absolute difference 
			deltaNa0[0] = FindNearestValue(nbrA, (int)ceil((double)nbrParticles/2.0), 
			                              (int)floor((double)nbrParticles/2.0));
			occupations[0] = nbrA;
			allowedExcitationSize[0] = nbrParticles-nbrA;
			deltaNa0List.push_back(deltaNa0);
			nbrEachLL.push_back(occupations);
			deltaLza.push_back(0);
			maxExcitationSize.push_back(allowedExcitationSize);
			break;
		}
		case 2:		//	2-current case e.g. 2/5 composite fermions, nu=2 IQHE 
		{	
			//	The number in the cut is defined for each CF LL in the construction
			//	shifted by an equal cut in each CF LL

			//	Now generate all possible distributions of N_A amongst the levels
			//	Starting from the one with the lowest L_z^A (corresponding to "sector 0")
			//	and working up to include all those less than the sector we wish
			//	to diagonalize for
			std::vector<int> occupations(nbrColours);
			std::vector<int> baseOccupations(nbrColours);
			double baseLzA, currLzA;
			int nbrLowest = nbrParticles/2-1;		//	states in lowest LL
			baseOccupations[0] = (int)floor((double)nbrA/2.0);
			baseOccupations[1] = (int)ceil((double)nbrA/2.0);
			baseLzA = occupationBasis::FindLzA(baseOccupations,nbrLowest);		
			currLzA = baseLzA;
			occupations = baseOccupations;
			//	Try moving one at a time from the top to the bottom row
			while(floor(currLzA-baseLzA)<=lzA && occupations[0]<=nbrLowest && occupations[1]<=nbrLowest+2)
			{			
				std::vector<int> deltaNa0(nbrColours);
				deltaNa0[0] = occupations[0] - baseOccupations[0];
				deltaNa0[1] = occupations[1] - baseOccupations[1];
				std::vector<int> allowedExcitationSize(nbrColours);
				allowedExcitationSize[0] = nbrLowest - occupations[0];
				allowedExcitationSize[1] = nbrLowest + 2 - occupations[1];
				deltaLza.push_back((int)floor(currLzA-baseLzA));
				nbrEachLL.push_back(occupations);
				deltaNa0List.push_back(deltaNa0);
				maxExcitationSize.push_back(allowedExcitationSize);
				occupations[0]++;
				occupations[1]--;
				currLzA = occupationBasis::FindLzA(occupations,nbrLowest);
			}
			//	Now try moving one at a time from the bottom to the top row
			//	But don't include the baseOccupations case again
			occupations = baseOccupations;
			occupations[0]--;
			occupations[1]++;
			currLzA = occupationBasis::FindLzA(occupations, nbrLowest);
			while(floor(currLzA-baseLzA)<=lzA && occupations[0]<=nbrLowest && occupations[1]<=nbrLowest+2)
			{				
				std::vector<int> deltaNa0(nbrColours);
				deltaNa0[0] = occupations[0] - baseOccupations[0];
				deltaNa0[1] = occupations[1] - baseOccupations[1];
				std::vector<int> allowedExcitationSize(nbrColours);
				allowedExcitationSize[0] = nbrLowest - occupations[0];
				allowedExcitationSize[1] = nbrLowest + 2 - occupations[1];
				deltaLza.push_back((int)floor(currLzA-baseLzA));
				nbrEachLL.push_back(occupations);
				deltaNa0List.push_back(deltaNa0);
				maxExcitationSize.push_back(allowedExcitationSize);
				occupations[0]--;
				occupations[1]++;
				currLzA = occupationBasis::FindLzA(occupations, nbrLowest);
			}
			break;
		}
		case 3:		//	3-current case e.g. 3/7 composite fermions 
        {
		    //	The number in the cut is defined for each CF LL in the construction
			//	shifted by an equal cut in each CF LL
			//	Now generate all possible distributions of N_A amongst the levels
			//	Starting from the one with the lowest L_z^A (corresponding to "sector 0")
			//	and working up to include all those less than the sector we wish
			//	to diagonalize for
			std::vector<int> baseOccupations(nbrColours);
			double baseLzA,currLzA;
			int nbrLowest = nbrParticles/3-2;		//	states in lowest LL
			baseOccupations[0] = (int)floor(((double)nbrA-1.0)/3.0);
			baseOccupations[1] = (int)floor(((double)nbrA-baseOccupations[0])/2.0);
			baseOccupations[2] = nbrA - baseOccupations[0] - baseOccupations[1];
			baseLzA = occupationBasis::FindLzA(baseOccupations, nbrLowest);
            //  Generate all possible occupations of the CF LLs, then
            //  filter out those with floor(currLzA-baseLzA)<=lzA
		    std::vector<std::vector<int> > partitionList;
		    std::vector<std::vector<int> > fullOccupationList;
			utilities::GenerateIntegerPartitions(partitionList, nbrA, nbrColours, nbrA);
			int nbrPartitions = partitionList.size();
			for(int i=0; i<nbrPartitions; ++i)
			{
			    partitionList[i].resize(nbrColours);
			    std::vector<std::vector<int> > uniquePermList;
			    utilities::UniqueObjectPermutations<int>(uniquePermList, partitionList[i]);
			    //  Check to see if the occupations list is consistent with 
			    //  the maximum allowed occupations of the set of CF LLs
			    int nbrPermutations = uniquePermList.size();
			    for(int j=0; j<nbrPermutations; ++j)
			    {
			        std::vector<int> newOccupation = uniquePermList.back();
			        if(newOccupation[0]<=nbrLowest && newOccupation[1]<=nbrLowest+2 
			            && newOccupation[2]<=nbrLowest+4)
			        {
			            fullOccupationList.push_back(newOccupation);
			        }
			        uniquePermList.pop_back();
			    }
			}
			nbrPartitions = fullOccupationList.size();
			for(int i=0; i<nbrPartitions; ++i)
			{
                std::vector<int> occupations(nbrColours);
			    occupations[0] = fullOccupationList[i][0];
			    occupations[1] = fullOccupationList[i][1];
			    occupations[2] = fullOccupationList[i][2];
			    currLzA = occupationBasis::FindLzA(occupations,nbrLowest);
			    if(floor(currLzA-baseLzA)<=lzA)
			    {
			        deltaLza.push_back((int)floor(currLzA-baseLzA));
			        nbrEachLL.push_back(occupations);
				    std::vector<int> deltaNa0(nbrColours);
			        deltaNa0[0] = occupations[0] - baseOccupations[0];
				    deltaNa0[1] = occupations[1] - baseOccupations[1];
				    deltaNa0[2] = occupations[2] - baseOccupations[2];
				    deltaNa0List.push_back(deltaNa0);
				    std::vector<int> allowedExcitationSize(nbrColours);
				    allowedExcitationSize[0] = nbrLowest - occupations[0];
				    allowedExcitationSize[1] = nbrLowest + 2 - occupations[1];
				    allowedExcitationSize[2] = nbrLowest + 4 - occupations[2];
				    maxExcitationSize.push_back(allowedExcitationSize);
				}
			}
		}
		break;
	}
	return;
}

//!
//! Generate a list of the possible coloured integer partitions that 
//! correspond to all possible occupations of the Landau level orbitals
//!
void occupationBasis::GenerateLevelOccupations(
    std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,   
                                //!<    List of LL occupations
    std::vector<int>& spectrumLabels,         
                                //!<    List of spectrum labels
    const int nbr,              //!<    Number of particles
    const int nbrA,             //!<    Number in the A region
    const int nbrLevels,        //!<    Number of Landau levels
    const int sector)           //!<    Given angular momentum sector
{    
    //  Generate all possible coloured integer partitions corresponding to all
    //  possible occupations of the single particle orbitals in this sector
    //  and convert the sets of integer partitions into populations of the 
    //  single particle orbitals. 
    //  Then, convert these integer partitions into Landau level
    //  orbital occupation signatures
    //  We start with the lowest lza state (nbrA bits out of nbr
    //  bits occupied, form the left). The integer partitions 
    //  are interpreted as the set of all possible excitations
    //  from the lowest lza state. For consistency, we bit shift
    //  the leftmost bits by each value in the integer partition
    //  starting from the largest excitations first.
    //  The colour label in the partition corresponds to the Landau
    //  level where the excitation is applied
    //  First generate all possible occupations of the Landau levels
    //  (reusing the Hilbert space data function from the HilbertSpace class)
    std::vector<std::vector<int> > deltaNa0List;		//  List of minimum (Na-Na0)
    std::vector<std::vector<int> > nbrEachLL;	        //  List of number in each LL       
    std::vector<int> deltaLza;						    //	Difference in Lz compared
												        //	to Na0 cut
    std::vector<std::vector<int> > maxExcitationSize;   //  Maximum allowed excitation in 
											            //	a given LL
    //	Note Na0 is defined to be the particle cut corresponding to the minimum 
    //	pseudo-energy and such that (Na-Na0) is minimal in cases where Na0 is 
    //  degenerate.
    occupationBasis::GenerateOccupationData(nbr, nbrA, sector, nbrLevels, deltaNa0List,
                                            nbrEachLL, deltaLza, maxExcitationSize);
    const int nbrBlocks = deltaNa0List.size();
    char dummyColours[nbrLevels];
    dummyColours[0] = 'A';
    if(nbrLevels>1)	
    {
        dummyColours[1] = 'B';
    }
    if(nbrLevels>2)	
    {
        dummyColours[2] = 'C';
    }
    if(nbrLevels>3)
    {
	    std::cerr << "\tERROR: Not implemented for > 3 Landau levels." << std::endl;
	    exit(EXIT_FAILURE);
    }
    for(int b=0; b<nbrBlocks; ++b)
    {
        std::vector<std::vector<int> > partitionList;
        utilities::GenerateColouredIntegerPartitions(partitionList, sector-deltaLza[b], nbrA,
                                                     maxExcitationSize[b], nbrEachLL[b], nbrLevels, 
                                                     dummyColours);
        const int dimList = partitionList.size();
        //  Generate the lowest state based on the nbrEachLL value
        //  Bitshift left for the lower Landau levels e.g.
        //
        //  0 0 0 1 1 1 1 1
        //  0 0 0 1 1 1 1 0
        //  0 0 0 1 1 1 0 0
        std::vector<std::bitset<_MAX_NUMBER_> > lowestOccupations(nbrLevels);
        for(int i=0; i<nbrLevels; ++i)
        {
            lowestOccupations[i] = std::bitset<_MAX_NUMBER_>(0);
            for(int j=nbrLevels-1-i; j<nbrEachLL[b][i]+nbrLevels-1-i; ++j)
            {
                lowestOccupations[i].flip(j);
            }
        }
        int blockLabel = 0;
        for(int l=0; l<nbrLevels; ++l)
        {   
            blockLabel += nbrEachLL[b][l]*(l);
        }
        if(dimList==0)  //  Include only ground state configuration
        {
            occupationsList.push_back(lowestOccupations);
            spectrumLabels.push_back(blockLabel);
        }
        for(int i=0; i<dimList; ++i)
        {
            std::vector<std::bitset<_MAX_NUMBER_> > currOccupations = lowestOccupations;
            const int size = partitionList[i].size()/2;
            std::vector<int> excitationCounter(nbrLevels);
            for(int k=0; k<size; ++k)
            {
                for(int c =0; c<nbrLevels; ++c)
                {
                    if(partitionList[i][k+size] == dummyColours[c])
                    {
                        //  Get position of left most set bit, not inducing the previous excitations
                        //  of the same colour (which is taken into account by the excitationCounter)
                        std::bitset<_MAX_NUMBER_> leftMost = std::bitset<_MAX_NUMBER_>(1) << (nbrEachLL[b][c]+nbrLevels-c-2-excitationCounter[c]);
                        currOccupations[c] ^= leftMost;
                        leftMost <<= partitionList[i][k];
                        currOccupations[c] |= leftMost;
                        ++excitationCounter[c];
                    }
                }
            }
            //  Append the current occupation to the full list
            occupationsList.push_back(currOccupations);  
            spectrumLabels.push_back(blockLabel);
        }
    }
    return;
}

//!
//! Assign the a set of single particle entanglement energy levels 
//! to the set of orbital occupations
//!
void occupationBasis::AssignEnergies(
    std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,   
                              //!<    List of LL occupations
    double* energyLevels,     //!<    Set of single particle energy levels
    double* spectrum,         //!<    Data structure to store the entanglement spectrum
    const double norm,        //!<    A value for the overall normalization
    const int nbrStates,      //!<    Number of angular momentum eigenstates
    const int nbrLevels)      //!<    Number of Landau levels
{
    double* p_spectrum = spectrum;
    for(unsigned int i=0; i<occupationsList.size(); ++i, ++p_spectrum)
    {
        double entanglementEnergy = 0;
        for(int j=0; j<nbrLevels; ++j)
        {
            for(int k=0; k<nbrStates; ++k)
            {
                entanglementEnergy += occupationsList[i][j][k]*energyLevels[k*nbrLevels+j];
            }
        }
        entanglementEnergy -= norm;
        utilities::cout.DebuggingInfo() << "\tTotal entanglement energy " << entanglementEnergy << std::endl;
        *(p_spectrum) = entanglementEnergy;
    }
    return;
}

//!
//! Assign the a set of single particle entanglement energy levels 
//! to the set of orbital occupations, then square the result
//!
void occupationBasis::AssignSquaredEnergies(
    std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,   
                              //!<    List of LL occupations
    double* energyLevels,     //!<    Set of single particle energy levels
    double* spectrum,         //!<    Data structure to store the entanglement spectrum
    const double norm,        //!<    A value for the overall normalization
    const int nbrStates,      //!<    Number of angular momentum eigenstates
    const int nbrLevels)      //!<    Number of Landau levels
{
    double* p_spectrum = spectrum;
    for(unsigned int i=0; i<occupationsList.size(); ++i, ++p_spectrum)
    {
        double entanglementEnergy = 0;
        for(int j=0; j<nbrLevels; ++j)
        {
            for(int k=0; k<nbrStates; ++k)
            {
                entanglementEnergy += occupationsList[i][j][k]*energyLevels[k*nbrLevels+j];
            }
        }
        entanglementEnergy -= 2*norm;
        utilities::cout.DebuggingInfo() << "\tTotal entanglement energy " << entanglementEnergy << std::endl;
        *(p_spectrum) += entanglementEnergy*entanglementEnergy;
    }
    return;
}
