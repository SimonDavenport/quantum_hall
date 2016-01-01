////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 15/02/2015
//!
//!  \file 
//!		This cpp file contains functions to evaluate the Z-k parafermion ground
//!     state wave function in the sphere geometry
//!
//!                    Copyright (C) 2015 Simon C Davenport
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

#include "parafermion.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the Parafermion class. 
//!
//!	 Assigns basic wave function characteristics. 
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Parafermion::Parafermion(
	FQHE::WaveFunctionData *wavefunctionData) 		//!<	Address of a struct of wave function data
    : 
    FQHE::WaveFunction(wavefunctionData),
    m_nCluster(0),
    m_nLastCluster(0),
	m_pfaff(0)
{
	//	Generate basic wave function data 
    
    m_wfData->fillNumerator     = m_wfData->nbrClusters;
    m_wfData->fillDenominator   = m_wfData->nbrClusters*m_wfData->jastrowExponent+2;
    m_wfData->shift             = m_wfData->jastrowExponent+2;
    
    if(m_wfData->nbr%m_wfData->nbrClusters!=0)
    {
        std::cerr<<"\n\tERROR in FQHE::Parafermion. Invalid particle number for cluster size"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    if(m_wfData->nbr>64)
    {
        std::cerr<<"\n\tERROR in FQHE::Parafermion. Cannot use more than 64 particles due to 64-bit binary implementation"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
    //  Construct wave function name
    {
        std::stringstream temp;
        temp.str("");
        temp<<"Z_"<<m_wfData->nbrClusters<<"_Parafermion";
        m_wfData->wfName = temp.str();
    }
    
    m_wfData->statistics = _FERMIONS_;

	m_wfData->fillingFactor =(double)m_wfData->fillNumerator/m_wfData->fillDenominator;
	
	if(m_wfData->geometry==_SPHERE_ || m_wfData->geometry==_DISC_)
	{
		m_wfData->flux=(double)m_wfData->nbr/(m_wfData->fillingFactor)-m_wfData->shift;
		m_wfData->monopoleStrength = (m_wfData->nbr*m_wfData->fillDenominator/m_wfData->fillNumerator-m_wfData->shift)/2.0;
	}
    
    //  Now we can generate a file name
	m_wfData->GenerateFileName();

    //  Generate sphere/disc radius value
    m_wfData->InitRadius();

	//  Print a summary
	m_wfData->CheckAndPrint();
	
	//  Generate all independent permutations of index labels. Only labels swapped between
	//  different clusters will contribute

    m_nCluster = m_wfData->nbr/m_wfData->nbrClusters;

    //  Note that we can exploit a factorization in terms of Pfaffian factors to reduce
	//  the number of clusters by a factor of 2 (for even number of clusters) or
	//  by a factor of 2 excluding the final cluster (for odd number of clusters)
    
    //  For odd number of clusters, the last cluster will remain the same size
    m_nLastCluster = m_wfData->nbrClusters&1 ? m_nCluster : 0;
    
    //  Increase the number of particles in each remaining cluster by a factor of 2
    m_nCluster *= 2;
    
    //  Decrease the number of clusters by the appropriate factor (not including the last one)
    m_wfData->nbrClusters = (m_wfData->nbrClusters)/2;
    
    //  Allocate memory to construct a Pfaffian for each new cluster (not including the last one)
    m_pfaff = new dcmplx[m_nCluster*m_nCluster];

    //  Recursively generate subsets of the cluster permutations

	std::vector<unsigned char> originalList(m_wfData->nbr);

	for(unsigned char i=0;i<m_wfData->nbr;++i)
	{
	    originalList[i] = i;
	}
	
	std::vector<unsigned int> nClusterList;

	for(int i=0;i<m_wfData->nbrClusters;++i)
	{
	    nClusterList.push_back(m_nCluster);
	}
	
	if(m_nLastCluster)  nClusterList.push_back(m_nLastCluster);

	//  Allocate temporary working memory
	std::vector<unsigned char> currPerm(m_wfData->nbr);

	this->BuildSubset(m_wfData->nbr,&nClusterList,&originalList,&currPerm,&m_clusterPermutations,0);
	
	std::cout<<"\tTotal number of cluster permutations = "<<m_clusterPermutations.size()<<std::endl<<std::endl;
	
	#if _DEBUG_
	
	//  Print values
	
	for(auto& it : m_clusterPermutations)
	{
	    std::cout<<"\t";
	
	    for(auto& it_it : it)
	    {
	        std::cout<<it_it<<" ";
	    }
	    
	    std::cout<<std::endl;
	}
	
	#endif
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Recursively generate unique subsets of the particle indices in 
//! different clusters
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::Parafermion::BuildSubset(
    const unsigned int nbrRemain,       //!<    Number of indices remaining to be assigned to clusters
    std::vector<unsigned int>* nClusterList,
                                        //!<    List of the number of particles in each cluster
    std::vector<unsigned char>* remainList, 
                                        //!<    List of indexes to be assigned to clusters
    std::vector<unsigned char>* currPerm,     
                                        //!<    Current permutation being constructed (of size nbr)
    std::vector<std::vector<unsigned char> >* clusterPermutations,
                                        //!<    Full list of cluster permutations
    const unsigned int level)           //!<    Recursion level (= which
                                        //!     cluster we're assigning)
    const
{
    //  Test with full set of permutations
    #if 0
    do
    {
        clusterPermutations->push_back(*remainList);
    }
    while(std::next_permutation(remainList->begin(),remainList->end()));
    
    for(auto& it : *clusterPermutations)
    {
        std::cout<<"\tADDED ";
	
        for(auto& it_it : it)
        {
            std::cout<<it_it<<" ";
        }
        
        std::cout<<std::endl;
    }
    
    #endif
   
    #if 1

    unsigned int nCluster         = (*nClusterList)[level];
    uint64_t mask                 = utilities::binary::FirstBinaryHammingNumber64(nCluster);
    const unsigned int nbrSubsets = utilities::BinomialFromTable(nbrRemain,nCluster);
    unsigned int offset           = 0;
    
    for(unsigned int i=0;i<level;++i)
    {
        offset += (*nClusterList)[i];
    }

    //PRINT("mask",mask);
    //PRINT("nbrSubsets",nbrSubsets);

    //  Generate subsets of the original list by masking it with a list of binary Hamming
	//  numbers in lexicographic order

    for(unsigned int i=0;i<nbrSubsets;++i)
    {
        uint64_t temp = mask;
           
        //PRINT("mask",std::bitset<64>(mask));
        //PRINT("offset",offset);

        auto it = currPerm->begin()+offset;

        while(temp!=0)
        {
            //PRINT("temp",std::bitset<64>(temp));
        
            //  Isolate the right-most bit using a 2s complement, then 
            //  use its position to mask the original list

            uint64_t occupied = temp & -temp;
            *it = (*remainList)[utilities::binary::HammingWeight64(occupied-1)];
            
            //PRINT("*it = ",(*remainList)[utilities::binary::HammingWeight64(occupied-1)]);

            //  Remove right-most bit and proceed to the next iteration

            temp = temp ^ occupied;
            ++it;
        }

        //  Put remaining unused indices in a separate list

        temp = (~mask) & utilities::binary::FirstBinaryHammingNumber64(nbrRemain);

        std::vector<unsigned char> subRemainList(nbrRemain-nCluster);
        
        it = subRemainList.begin();
        
        while(temp!=0)
        {
            //PRINT("temp",std::bitset<64>(temp));
            
            uint64_t occupied = temp & -temp;
            *it = (*remainList)[utilities::binary::HammingWeight64(occupied-1)];
            
            temp = temp ^ occupied;
            ++it;
        }
        
        //  Append to main list
        if(level == nClusterList->size()-1)
        {
            clusterPermutations->push_back(*currPerm);
            
            //std::cout<<"\tADDED ";
	
	        //for(auto& it : *currPerm)
	        //{
	        //    std::cout<<it<<" ";
	        //}
	        
	        //std::cout<<std::endl;
            
            //getchar();
        }
        else
        {
            this->BuildSubset(nbrRemain-nCluster,nClusterList,&subRemainList,currPerm,clusterPermutations,level+1);
        }
        
        //  Generate the next permutation
	    mask = utilities::binary::NextHammingNumber64(mask);
	    
	    //PRINT("next mask",std::bitset<64>(mask));
    }
    
    #endif
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the Parafermion class. 
//!
//!	 Deallocates memory as required
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Parafermion::~Parafermion()
{
    if(m_pfaff!=0)  delete[] m_pfaff;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
#if _ENABLE_SPHERE_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Parafermion wave function algorithm in the sphere geometry 
//!
//!	For more information see 
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::Parafermion::EvaluateWfSphere(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* u,	    //!<	The memory address of the start of an array of spinor u coordinates
	dcmplx* v)	    //!<	The memory address of the start of an array of spinor v coordinates
	const
{
    dcmplx wfValue = 0.0;
    bool isFirstPerm = true;

    //  Note that we can exploit a factorization in terms of Pfaffian factors to reduce
	//  the number of clusters by a factor of roughly 2.

    for(auto& it_cluster : m_clusterPermutations)
    {
        dcmplx accum = 0.0;
    
        //  Build a Pfaffian factor and Jastrow factor for each cluster permutation
        
        for(int a=0;a<m_wfData->nbrClusters;++a)
        {
            dcmplx* p_row = m_pfaff;

            for(unsigned int i=0;i<m_nCluster;++i,p_row+=m_nCluster)
            {
                for(unsigned int j=i+1;j<m_nCluster;++j)
                {
	                *(p_row+j)=dcmplx(1,0)/(*(u+it_cluster[a*m_nCluster+i]) * *(v+it_cluster[a*m_nCluster+j]) - *(v+it_cluster[a*m_nCluster+i]) * *(u+it_cluster[a*m_nCluster+j]));
                }
            }
            
            accum += utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,m_nCluster);
            
            for(unsigned int i=0;i<m_nCluster;++i)
            {
                for(unsigned int j=i+1;j<m_nCluster;++j)
                {
                    accum += log(*(u+it_cluster[a*m_nCluster+i]) * *(v+it_cluster[a*m_nCluster+j]) - *(v+it_cluster[a*m_nCluster+i]) * *(u+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //  Include the last cluster in terms of the standard pairing formulation (k odd only)
        
        if(m_nLastCluster>0)
        {
            unsigned int a = m_wfData->nbrClusters;
        
            for(unsigned int i=0;i<m_nLastCluster;++i)
            {
                for(unsigned int j=i+1;j<m_nLastCluster;++j)
                {
                    accum += 2.0*log(*(u+it_cluster[a*m_nCluster+i]) * *(v+it_cluster[a*m_nCluster+j]) - *(v+it_cluster[a*m_nCluster+i]) * *(u+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //	Calculate running total using logarithmic addition

        if(isFirstPerm)
        {
            isFirstPerm = false;
	        wfValue     = accum;
        }
        else
        {
	        wfValue = log(exp(wfValue-accum)+1.0)+accum;
        }
    }

    //  Full method
    #if 0
    for(auto& it_cluster : m_clusterPermutations)
    { 
        dcmplx accum = 0.0;
        
        //  Inter-cluster terms
        
        for(int a=0;a<m_wfData->nbrClusters;++a)
        {
            for(int i=0;i<m_nCluster;++i)
            {
                for(int j=i+1;j<m_nCluster;++j)
                {
                    accum += 2.0*log(*(u+it_cluster[a*m_nCluster+i]) * *(v+it_cluster[a*m_nCluster+j]) - *(v+it_cluster[a*m_nCluster+i]) * *(u+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //	Calculate running total using logarithmic addition

	    if(isFirstPerm)
	    {
            isFirstPerm = false;
		    wfValue     = accum;
	    }
	    else
	    {
		    wfValue = log(exp(wfValue-accum)+1.0)+accum;
	    }
    }
    #endif
    
    //  Include Jastrow factor in front
    
    for(int i=0;i<n;++i)
    {
        for(int j=i+1;j<n;++j)
        {
            wfValue += (double)m_wfData->jastrowExponent*log(*(u+i) * *(v+j) - *(v+i) * *(u+j));
        }
    }
    
    //PRINT("wfValue",wfValue);
    //getchar();
    
    return wfValue;
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////
#if _ENABLE_DISC_GEOMETRY_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Parafermion wave function algorithm in the disc geometry 
//!
//!	For more information see 
//!
////////////////////////////////////////////////////////////////////////////////

dcmplx FQHE::Parafermion::EvaluateWfDisc(
	const int n,	//!<	The number of particles to perform the evaluation for
	dcmplx* z)	    //!<	The memory address of the start of an array of disc coordinates
	const
{
    dcmplx wfValue = 0.0;
    bool isFirstPerm = true;
    
    //  Note that we can exploit a factorization in terms of Pfaffian factors to reduce
	//  the number of clusters by a factor of roughly 2.

    for(auto& it_cluster : m_clusterPermutations)
    {
        dcmplx accum = 0.0;
    
        //  Build a Pfaffian factor and Jastrow factor for each cluster permutation
        
        for(int a=0;a<m_wfData->nbrClusters;++a)
        {
            dcmplx* p_row = m_pfaff;

            for(unsigned int i=0;i<m_nCluster;++i,p_row+=m_nCluster)
            {
                for(unsigned int j=i+1;j<m_nCluster;++j)
                {
	                *(p_row+j)=dcmplx(1,0)/(*(z+it_cluster[a*m_nCluster+i]) - *(z+it_cluster[a*m_nCluster+j]));
                }
            }
            
            accum += utilities::linearAlgebra::LogPfaffian<dcmplx>(m_pfaff,m_nCluster);
            
            for(unsigned int i=0;i<m_nCluster;++i)
            {
                for(unsigned int j=i+1;j<m_nCluster;++j)
                {
                    accum += log(*(z+it_cluster[a*m_nCluster+i]) - *(z+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //  Include the last cluster in terms of the standard pairing formulation (k odd only)
        
        if(m_nLastCluster>0)
        {
            unsigned int a = m_wfData->nbrClusters;
        
            for(unsigned int i=0;i<m_nLastCluster;++i)
            {
                for(unsigned int j=i+1;j<m_nLastCluster;++j)
                {
                    accum += 2.0*log(*(z+it_cluster[a*m_nCluster+i]) - *(z+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //	Calculate running total using logarithmic addition

        if(isFirstPerm)
        {
            isFirstPerm = false;
	        wfValue     = accum;
        }
        else
        {
	        wfValue = log(exp(wfValue-accum)+1.0)+accum;
        }
    }
    
    //  Full method
    #if 0
    for(auto& it_cluster : m_clusterPermutations)
    {
        dcmplx accum = 0.0;
    
        //  Inter-cluster terms
        
        for(int a=0;a<m_wfData->nbrClusters;++a)
        {
            for(int i=0;i<m_nCluster;++i)
            {
                for(int j=i+1;j<m_nCluster;++j)
                {
                    accum += 2.0*log(*(z+it_cluster[a*m_nCluster+i]) - *(z+it_cluster[a*m_nCluster+j]));
                }
            }
        }
        
        //	Calculate running total using logarithmic addition

		if(isFirstPerm)
		{
            isFirstPerm = false;
			wfValue     = accum;
		}
		else
		{
			wfValue=log(exp(wfValue-accum)+1.0)+accum;
		}
    }
    #endif

    //  Include Jastrow factor in front
    
    for(int i=0;i<n;++i)
    {
        for(int j=i+1;j<n;++j)
        {
            wfValue+= (double)m_wfData->jastrowExponent*log(*(z+i) - *(z+j));
        }
    }
    
	for(int i=0;i<n;++i)
	{
		//	Implement Gaussian exponential part
		double tmp = abs(*(z+i));
		wfValue -= tmp*tmp/4.0;	
	}
	
	return wfValue;	
}

#endif
//////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
