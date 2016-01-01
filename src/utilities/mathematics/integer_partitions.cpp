////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 01/12/2014
//!
//!  \file 
//!  	This file implements integer partition functions.
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

#include "integer_partitions.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function generates a full set of integer partitions and then
//! filters the list down to (k,r) admissible partitions. 
//!
//! (k,r) admissible partitions satisfy:
//!
//! \[ \lambda_i - \lambda_{i+k} \ge r  \]
//!
//! Note that this algorithm is very inefficient because it first generates
//! the much larger set of integer patitions.
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::Generate_k_r_AdmissiblePartitions(
    partition_t& krAdmissible,
                        //!<    A list of k-r admissible partitions to be populated
    const int n,        //!<    Number to be partitioned
    const int p,        //!<    Size of partition
    const int k,        //!<    k index in admissibility condition
    const int r,        //!<    r index in admissibility condition
    const bool shift)   //!<    the shift flag tells us to allow for partitions at n=1
{
    //std::cout<<"\n\t ============ GENERATING (k,r) ADMISSIBLE PARTITIONS ============\n"<<std::endl;

    //  First generate a list of all integer patitions

    partition_t partitionList;
    
    utilities::GenerateIntegerPartitions(partitionList,n,p,n);
    
    //  Treat the paritions as all of length p e.g. for n=p=3:
    //
    //  3 0 0
    //  2 1 0
    //  1 1 1

    //  Iterate over the list and check k-r admissibility of each term
    
    int counter = 0;
    int j=0;
    
    for(auto it_partition=partitionList.begin();it_partition<partitionList.end();it_partition++,j++)
    {
        const int dimPartition = it_partition->size();
        
        bool admissible = true;
        
        //std::cout<<"\tTEST";  
        //for(int q=0;q<dimPartition;q++)
        //{
        //    std::cout<<partition[q]<<" ";
        //}
        //std::cout<<std::endl;
        
        for(int i=0;i<dimPartition;i++)
        {
            if(i+k<dimPartition) //  Check condition within the list
            {
                admissible = admissible && ((*it_partition)[i+k] - (*it_partition)[i] >= r);

                //if(!admissible) std::cout<<"FAILED CHECK 1"<<std::endl;
            }
            
            if(!admissible) break;
            
            if(k-i<=p-dimPartition)  //  Check the condition with trailing 0s
            {
                admissible = admissible && ((*it_partition)[i] >= r);
                
                //if(!admissible) std::cout<<"FAILED CHECK 2"<<std::endl;
            }
            else
            {
                admissible = false;
                
                //std::cout<<"FAILED CHECK 3"<<std::endl;
            }
               
            if(!admissible) break;
        }
        
        //  If admissible, then append to the list the be returned
        
        if(admissible)
        {
            //std::cout<<"\t ADDED No. "<<counter<<"\t";
            //for(int q=0;q<dimPartition;q++)
            //{
            //    std::cout<<it_partition[q]<<" ";
            //}
            //std::cout<<std::endl;
            
            krAdmissible.push_back(*it_partition);
            counter++;
        }
    }
    
    //  Allow a partition at n=1 if the shift flag is set
    
    if(shift && 1 == n && p >= 1)
    {
        std::vector<int> dummy;
    
        dummy.push_back(1);
    
        krAdmissible.push_back(dummy);
        counter++;
    }
    
    //std::cout<<"\n\t ============ DONE ============\n"<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function returns a list of integer partitions of the number n into
//!	no more than p parts and nbrColours colours.
//!	e.g. for n=3,p=2 and 1 colour it would return {3,J},{2,1,J,J}
//!
//!	This function acts as a wrapper for the ColouredPartitionIterator function. 
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::GenerateColouredIntegerPartitions(
    partition_t& partitions,        //!<    A list of integer partitions to be populated
	const int n,		            //!<	Number to be partitioned
	const int p,	                //!<	Maximum number of parts in the partition
	std::vector<int>& maxPartSize,	//!<	A list of maximum partition sizes
	                                //!     allowed for each colour
	std::vector<int>& maxColNbr,    //!<	A list of the maximum number of each 
	                                //!     colour allowable in each parititon
	const int nbrColours,			//!<	Number of colours.
	const char*  colours)			//!<	Character array containing colour 
	                                //!     labels. Must be of length nbrColours.
{
	int temp[n];

	ColouredPartitionIterator(n,p,maxPartSize,maxColNbr,temp,0,nbrColours,partitions,colours);
	
	//  Sort each of the integer partitions into descending order

	for(auto& part : partitions)
	{
	    unsigned int size = part.size()/2;
	    utilities::QuickSort<int,int,_DESCENDING_ORDER_>(&part[0],&part[size],size);
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function returns a list of integer partitions of the number n into
//!	no more than p parts. e.g. for n=3,p=2 it would return {3},{2,1}
//!
//!	This function acts as a wrapper for the more general PartitionIterator function.
//!
//! This version is an overload to make using a single coloured integer partition
//! function call simpler.
//!
//!	\return	A list of partitions
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::GenerateIntegerPartitions(
    partition_t& partitions,    //!<    A list of integer partitions to be populated
	const int n,	            //!<	Number to be partitioned
	const int p,                //!<	Maximum number of parts in the partition
	const int maxPartSize)      //!<    Maximum size of a number in a given partition
{
	int temp[n];

	PartitionIterator(n,p,maxPartSize,temp,0,partitions);
	
	//  Sort each of the integer partitions into descending order

	for(auto& part : partitions)
	{
	    std::sort(part.begin(),part.end());
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function returns a list of integer partitions of the number n 
//!	into exactly maxPartSize.size() parts, including all possible permutations. 
//! Limit the maximum partition size in each part to no more than maxPartSize.
//!	e.g. for n=4 and maxPartSize={3,4} it would return {3,1},{2,2},{1,3},{0,4}
//!
//!	This function acts as a wrapper for the Partition_Iterator function.
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::GenerateIntegerPartitionsWithPermutations(
    partition_t& partitions,        //!<    A list of integer partitions to be populated
    const int n,                    //!<	Number to be partitioned
    std::vector<int>& maxPartSize)  //!<	Maximum number allowed in each place 
                                    //!     in the partition 
{
	int temp[n];

    if(0==n)
    {
        //  Return an empty partition
        
        std::vector<int> dummy(maxPartSize.size(),0);
        
        partitions.push_back(dummy);
        
        return;
    }

    partition_t originalPartitions;

	PartitionIterator(n,maxPartSize.size(),n,temp,0,originalPartitions);
	
	//  Extend the list of partitions by including all unique permutations
    
	for(auto& part : originalPartitions)
	{
	    //  First pad the partition list with zeros
	
	    part.resize(maxPartSize.size(),0);
	
	    //  Sort each of the integer partitions into descending order
	
	    std::sort(part.begin(),part.end());

	    //  Then append all permutations of partitions to the list

        do
        {
            //  Test for maxPartSize constraint before appending
            //  to the list
            
            bool allowed = true;
            auto it_maxPart = maxPartSize.begin();
            
            for(auto it_part=part.begin();it_part<part.end();it_part++,it_maxPart++)
            {
                if(*it_part > *it_maxPart)  allowed = false;
            }
            
            if(allowed) partitions.push_back(part);
        } 
        while (std::next_permutation(part.begin(),part.end()));
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Iteratively calling this function generates a list of integer partitions
//!	of n into no more than p parts. The list is returned as a vector of 
//!	integer vectors
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::PartitionIterator(
	const int n,					//!<	Number to be partitioned
	const int p,					//!<	Maximum number of parts in the partition
	const int maxPartSize,          //!<    Maximum allowed number in the partition
	int* v,					        //!<	A vector of length n to provide working 
	                                //!     memory space 
	int level,						//!<	Iteration level (initial call with level=0)
	partition_t& partitionList)	    //!<	A vector of vectors to be populated with 
	                                //!     a list of all partitions generated here
{
    if(n<1) return;

	v[level]=n;
	
	//	Add the current partition(s) to the overall list
	
    std::vector<int> tempArr((level+1));

    for(int i=0;i<=level;i++)
    {
	    tempArr[i]=v[i];
    }
    
    //	Impose the maxPartSize constraint
    
    bool isAllowed=true;
    
    for(int i=0;i<=level;i++)
	{
	    isAllowed = isAllowed && tempArr[i]<=maxPartSize;
    }
    
    //  Append to the overall list of partitions
    
    if(isAllowed)
	{
        partitionList.push_back(tempArr);
    }
    
    //	Go to next level of iteration
	
	if(level<p-1)
	{
		int first=(level==0) ? 1 : v[level-1];

		for(int i=first;i<=n/2;i++)
		{
			v[level]=i;

			//	Call the next level of iteration
			PartitionIterator(n-i,p,maxPartSize,v,level+1,partitionList);
		}
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Iteratively calling this function generates a list of integer partitions
//!	of n into no more than p parts and nbrColours colours. The list is returned
//!	as a vector of a integer vectors e.g. a list containing {3;J},{3;K}
//!	{2,1;J,J},{2,1;J,K},{2,1;K,J},{2,1;K,K}, {1,1,1;J,J,J},{1,1,1;J,J,K},  
//!	{1,1,1;J,K,K},{1,1,1;K,K,K}
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::ColouredPartitionIterator(
	const int n,					//!<	Number to be partitioned
	const int p,					//!<	Maximum number of parts in the partition
	std::vector<int>& maxPartSize,  //!<	A list of maximum partition sizes  
	                                //!     allowed for each colour 
	std::vector<int>& maxColNbr,    //!<	A list of the maximum number of each 
	                                //!     colour allowable in each parititon
	int* v,				    	    //!<	A vector of length n to provide working 
	                                //!     memory space 
	int level,						//!<	Iteration level (initial call with level=0)
	const int nbrColours,			//!<	Number of colours [N.B it's expected 
	                                //!     that maxPartSize and maxColNbr are of this size]
	partition_t& partitionList,	    //!<	A vector of vectors to be populated with 
	                                //!     a list of all partitions generated here
	const char* colours)			//!<	Labels to assign colours of this 
	                                //!     partition. Used to distinguish in 
	                                //!     cases of  coloured partitions
{
	if(n<1) return;

	v[level]=n;

    //	Add the current partition(s) to the overall list

    std::vector<int> tempArr(2*(level+1));

    for(int i=0;i<=level;i++)
    {
	    tempArr[i]=v[i];
    }

    //  Generate all possible associated colour assignment for 
    //  the length (level+1) of the current integer partition

    std::vector<std::vector<char> > colourList;

    char tempChar[(level+1)];

    ColourIterator(nbrColours,(level+1),maxColNbr,0,0,tempChar,colourList,colours);

    //  For each possible colour assignment, generate a list
    //  of unique permitations of that colour where 'unqiue'
    //  means that that e.g. {1,1,1;J,J,K} = {1,1,1;J,K,J}
    //  and e.g. {2,1;J,K} != {2,1;K,J}

    UnqiueColourPerms(v,colourList,(level+1));

    //  Include all of the different colour combinations in the overall list

    for(unsigned int c=0;c<colourList.size();c++)
    {
	    for(int i=level+1;i<=2*level+1;i++)
	    {
		    tempArr[i]=colourList[c][i-level-1];
	    }
	
	    //	Impose the maxPartSize constraint
	
	    bool isAllowed=true;
	
	    for(int i=0;i<=level;i++)
	    {
            int col=0;
        
		    for(int j=0;j<nbrColours;j++)
		    {
			    if(tempArr[i+level+1]==colours[j])
			    {
				    col = j;
			    }
		    }
		
		    isAllowed = isAllowed && tempArr[i]<=maxPartSize[col];
	    }

	    if(isAllowed)
	    {
		    partitionList.push_back(tempArr); 
	    }
    }
    
	//	Go to next level of iteration
	
	if(level<p-1)
	{
		int first=(level==0) ? 1 : v[level-1];

		for(int i=first;i<=n/2;i++)
		{
			v[level]=i;

			//	Call the next level of iteration
			ColouredPartitionIterator(n-i,p,maxPartSize,maxColNbr,v,level+1,nbrColours,partitionList,colours);
		}
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Iteratively calling this function generates a list of colour 
//!	association of nColour different colours into nSlots labels. The list is
//!	returned as a vector of chars e.g. a list containing AAA, AAB, ABB,BBB
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::ColourIterator(
	int nbrColours,			    //!<	Number of colours to be used 
	const int nSlots,		    //!<	Number of slots to label
	std::vector<int>& maxColNbr,//!<	Maximum no. slots that can be filled 
	                            //!     with a given colour
	int firstSlot,				//!<	The first slot to fillingFactor at each 
	                            //!     level of iteration
	int level,					//!<	Iteration level (initial call with level=0)
	char*  v,					//!<	A vector of length n to provide working 
	                            //!     memory space (of size nSlots)
	std::vector<std::vector<char> >& colourList,
	                            //!<	A vector of chars to be populated with a 
	                            //!     list of all partitions generated here
	const char*  colours)		//!<	Labels to assign colours of this partition.
	                            //!     IMPORTANT: length of colours array MUST equal 
	                            //!     the highest value of nbrColours used.
{
	//	Set all the remaining colours to be the last colour in the list
	
	for(int i=firstSlot;i<nSlots;i++)
	{
		v[i] = colours[level];
	}
	
	//  Count the number of each colour in the list
	
	int totNbrCol = maxColNbr.size();
	
	std::vector<int > counter(totNbrCol);
	
	for(int j=0;j<totNbrCol;j++)
	{
	    counter[j]=0;
	}
	
	for(int i=0;i<nSlots;i++)
	{
	    for(int j=0;j<totNbrCol;j++)
	    {
	        if(v[i]==colours[j])
	        {
	            counter[j]++;
	        }
	    }
	}
	
	bool appendToList=true;
	
	for(int j=0;j<totNbrCol;j++)
	{
	    //	Impose a condition to have a maximum number of partitions in a given colour
	
	    appendToList = appendToList && (counter[j]<=maxColNbr[j]);
	}

    if(appendToList)
    {
	    //	Append current partition to the list of all partitions
	
	    std::vector<char> tempArr(nSlots);
	
	    for(int i=0;i<nSlots;i++)
	    {
		    tempArr[i] = v[i];
	    }
	    
	    //std::cout<<"Level "<<level<<" Added to list: "<<appendToList<<std::endl;
        //for(i=0;i<nSlots;i++)
        //{
        //    std::cout<<tempArr[i]<<" ";
        //}
        //std::cout<<std::endl;
        //getchar();

	    colourList.push_back(tempArr);
	}
	
	//	Go to next level of iteration
	
	if(nbrColours>1)
	{
		for(int i=firstSlot;i<nSlots;i++)
		{
			for(int j=firstSlot;j<i;j++)
			{
				v[j] = colours[level];
			}
			
			ColourIterator(nbrColours-1,nSlots,maxColNbr,i,level+1,v,colourList,colours);
		}
	}
	
	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	This function generates all unique permutations of a set of colours, 
//!  where the uniquness of the permutation is decided by association to a
//!  partition of number labels. e.g. {1,1,1;J,J,K} = {1,1,1;J,K,J}
//!  and e.g. {2,1;J,K} != {2,1;K,J} 
//!
//!  The list of colours 'colourList' will be appended with any additional 
//!  unique permutations.
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::UnqiueColourPerms(
	int* partition,			    //!<	An integer partition into nSlots parts	
	std::vector<std::vector<char> >& colourList,
	                            //!<	A list of possible colour assignments
	int nSlots)					//!<	Number of parts in the partition
{
	//  Determine the number of colour associations to deal with e.g. JJJ, JJK, JKK, KKK
	
	int nAssociations = colourList.size();

	//  For each associataion, calculate all permutations of the colours
	//  If the permutation is "unique" then add to the list colourList
	
	for(int j=0;j<nAssociations;j++)
	{
	    //  Declare local variables
	
	    std::vector<ColouredPartition> uniquePerms;    //  A list to store all unique
	                                                //  permutations generated
	
	    std::vector<char> currColPerm;              //  Container for the current 
	                                                //  colour permutation

        ColouredPartition currPerm(nSlots);            //  Container for the current 
                                                    //  coloured integer partition

	    currColPerm = colourList[j];
	    
	    //  Put the first unique permutation in the list
	
	    //  Make a vector containing the current permutation as a list of 
        //  e.g. 1J, 2J, 3K...
        
        for(int k=0;k<nSlots;k++)
        {
            currPerm.Set(k,partition[k]);
            currPerm.Set(k,currColPerm[k]);
        }
        
        //  Sort the vector

        std::sort(currPerm.BeginAddress(),currPerm.EndAddress(),ColourPermsComp);

        uniquePerms.push_back(currPerm);
        
        //std::cout<<"unique perms contains:"<<std::endl;
        //for(int q=0;q<uniquePerms.size();q++)
        //{
        //    uniquePerms[q].Print();
        //}

        //  Generate colour permutations

        //  To ensure that we start at the lexicographically first permutation
        //  before we start we need to sort the colour assignment
        
        std::sort(currColPerm.begin(),currColPerm.end());

        //  Return if we're already at the last lexicographic permutation
        //  note that this step also skips to the first permutation, since 
        //  the identity permutation is already taken into account

        do 
        {
            //  Test if this is a "unique permutation" or not
            
            //  Make a vector containing the current permutation as a list of 
            //  e.g. {1J, 2J, 3K...}
            
            for(int k=0;k<nSlots;k++)
            {
                currPerm.Set(k,currColPerm[k]);
            }   
            
            //  Sort the vector
            
            //std::cout<<"unsorted:"<<std::endl;
            //currPerm.Print();
            
            std::sort(currPerm.BeginAddress(),currPerm.EndAddress(),ColourPermsComp);
            
            //std::cout<<"sorted:"<<std::endl;
            
            //std::cout<<"colour assignmnet:"<<std::endl;
            //currPerm.Print();
            
            //  Check to see if this vector already appears in the list
            //  of unique permutaitons.

            //  if(std::binary_search(uniquePerms.BeginAddress(),uniquePerms.EndAddress(),currPerm))
            //
            //  NOTE AT THIS POINT WE WANT TO USE BINARY SEARCH (SINCE THIS PROVIDES
            //  A HIGHLY EFFICIENT SEARCH ALGORITHM), HOWEVER, SINCE IT'S NOT ALWAYS
            //  POSSIBLE TO DEFINE A DEFINITE ORDERING OF INTEGER PARTITIONS
            //  (C.F HASSE DIAGRAMS) WE MUST DIRECTLY COMPARE TO EVERY VALUE IN THE
            //  LIST

            if(currPerm.PermutationInList(uniquePerms))
            {
                 //  If it was not there then we have generated a new unique permutation
            
                 //std::cout<<"\tOUTCOME: was there"<<std::endl;
            }
            else
            {
                //std::cout<<"\tOUTCOME: was not there"<<std::endl;
            
                //  Include this colour permutation in the list of unique
                //  Colour permuations to use

                colourList.push_back(currColPerm);    
            
                //  And add this new permuation to the list of unique perms
                
                uniquePerms.push_back(currPerm);

                //std::cout<<"unique perms contains:"<<std::endl;
                //for(int q=0;q<uniquePerms.size();q++)
                //{
                //    uniquePerms[q].Print();
                //}
            }
        } 
        while ( std::next_permutation(currColPerm.begin(),currColPerm.end()));
        //getchar();
	}
	
	return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief  This function compares two terms of the form e.g. {1J, 2J, 3K...} 
//!  and returns a bool to determine whether these are in ascending lexicographic
//!  order or not. The order is defined by the integers e.g. 1<2, but if integers
//!  are equal, then the order is defined by the letters e.g. J<K (according to
//!  the numerical values associated to those characters)
//!	
//!  \return Return true if the terms are in descending lexicographic
//!  order, and false otherwise
//!
//////////////////////////////////////////////////////////////////////////////////

bool utilities::ColourPermsComp(
	const PartitionData &lhs,		//!<	One of the terms of compare 
	const PartitionData& rhs)		//!<	The other of the terms of compare
{
    if(lhs.integer<rhs.integer)
    {
        return true;
    }
    else if(lhs.integer==rhs.integer)
    {
        if(lhs.colour<rhs.colour)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Constructor for an object of type ColouredPartition
//!	
////////////////////////////////////////////////////////////////////////////////

utilities::ColouredPartition::ColouredPartition(
	int dim)	//!<	Number of parts in the partition
{
	m_dim=dim;

	//std::cout<<"CALLED CTOR"<<std::endl;

	m_list = new PartitionData[m_dim];
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Destructor for an object of type ColouredPartition
//!	
////////////////////////////////////////////////////////////////////////////////

utilities::ColouredPartition::~ColouredPartition()
{
	delete[] m_list;
	
	//std::cout<<"CALLED DTOR"<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Copy constructor for an object of type ColouredPartition
//!	
////////////////////////////////////////////////////////////////////////////////

utilities::ColouredPartition::ColouredPartition(
	const ColouredPartition& other)	//!<	Address of the object to be coppied
{
	m_dim = other.m_dim;

	m_list = new PartitionData[m_dim];
	
	for(int i=0;i<m_dim;i++)
	{
		m_list[i]=other.m_list[i];
	}
	
	//std::cout<<"CALLED COPY CTOR"<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Prints out the current value of the coloured integer partition stored
//!	
////////////////////////////////////////////////////////////////////////////////

void utilities::ColouredPartition::Print()
{
	for(int i=0;i<m_dim;i++)
	{
		utilities::cout.AdditionalInfo()<<m_list[i].integer<<m_list[i].colour<<" ";
	}
	utilities::cout.AdditionalInfo()<<std::endl;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Compares all terms in the search list and see if any are equal
//!	to any coloured integer parition stored in this object. 
//!	\return true if the provided list does contain the value, false otherwise
//!	
////////////////////////////////////////////////////////////////////////////////

bool utilities::ColouredPartition::PermutationInList (
	std::vector<ColouredPartition> searchList)	//!<	List of coloured integer 
	                                            //!     partitions to check against
{
	for(unsigned int j=0;j<searchList.size();j++)
	{

		//std::cout<<"Compare: ";
		//Print();
		//std::cout<<"With: ";
		//searchList[j].Print();

		bool areEqual = true;
	
		for(int i=0;i<m_dim;i++)
		{
			areEqual = areEqual && m_list[i].integer==searchList[j].m_list[i].integer
			             && m_list[i].colour==searchList[j].m_list[i].colour;
		
			//std::cout<<"So far "<<(m_list[i].integer==searchList[j].m_list[i].integer 
			//  && m_list[i].colour==searchList[j].m_list[i].colour)<<std::endl;
		}
		
		//std::cout<<"Conclude "<<areEqual<<std::endl;getchar();
		
		if(areEqual)
		{
			return true;
		}
	}
	
	return false;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Set a certian integer value in the coloured integer partition array
//!	
////////////////////////////////////////////////////////////////////////////////

void utilities::ColouredPartition::Set(
	int index,		//!<	Index of value to set
	int integer)	//!<	Integer to set it to 
{
	m_list[index].integer = integer;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Set a certain colour value in the coloured integer partition array
//!	
////////////////////////////////////////////////////////////////////////////////

void utilities::ColouredPartition::Set(
	int index,		//!<	Index of value to set
	char colour)			//!<	Colour to set it to 
{
	m_list[index].colour=colour;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Return the first memory address in the list
//!	\return First memory address in the list
////////////////////////////////////////////////////////////////////////////////

utilities::PartitionData* utilities::ColouredPartition::BeginAddress()
{
	return m_list;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Return the last memory address in the list
//!	\return Last memory address in the list
////////////////////////////////////////////////////////////////////////////////
		
utilities::PartitionData* utilities::ColouredPartition::EndAddress()
{
	return m_list+m_dim;
}

////////////////////////////////////////////////////////////////////////////////
#if 0

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Iteratively calling this function generates a list of the number of
//!	integer partitions of the number n into no more than p parts. This function
//!	populates a list of the length of each partition.
//!	
//////////////////////////////////////////////////////////////////////////////////

void utilities::CountIterator(
	 int n,			        //!<	Number to be partitioned
	 const int p,		    //!<	Maximum number of parts in the partition
	 int* v,			    //!<	A vector of length n to provide working memory
	 int level,		        //!<	Iteration level (initial call with level=0)
	 const int nbrColours,  //!<	Number of colours   
	 std::vector<int>& nbrTerms)    
	                        //!<	A list which will be populated with the 
	                        //!     number of terms in each integer partition 
{
	if(n<1) return;

	v[level]=n;
	
	int colCount=0;
	
	//  Generate all possible associated colour assignment associated with 
	//  the length (level+1) of the current integer partition

	ColourCountIterator(nbrColours,(level+1),0,0,&colCount);
	
	//  For each possible colour assignment, generate a list
	//  of unique permitations of that colour where 'unqiue'
	//  means that that e.g. {1,1,1;J,J,K} = {1,1,1;J,K,J}
	//  and e.g. {2,1;J,K} != {2,1;K,J}
	//  Increase colCount accordingly
	
	//Unqiue_Colour_Perms(v,colourList,(level+1));
	
	std::cout<<"Counted "<<colCount<<std::endl;

	//	Add the length of the current partition(s) to the list of all lengths
	//  Record duplicate lengths to take into account all colour combinations
	
	for(int c=0;c<colCount;c++)
	{
	    nbrTerms.push_back((level+1));
	}

	//	Go to next level of iteration

	if(level<p-1)
	{
		int first=(level==0) ? 1 : v[level-1];
		
		for(int i=first;i<=n/2;i++)
		{
			v[level]=i;
			//	Call the next level of iteration
			CountIterator(n-i,p,v,level+1,nbrColours,nbrTerms);
		}
	}
	
	return;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if 0

//////////////////////////////////////////////////////////////////////////////////
//!	\brief This function returns a list of the lengths of the set of integer
//!	partitions of the number n into no more than p parts and nbrColours colours.
//!	e.g. for n=3,p=2 and 1 colour it would return {1,2} (counting the length of
//!	partitions {3} and {2,1}).
//!
//!	This function acts as a wrapper for the CountIterator function.
//!
//!	\return	A list of the length of each partition
//!	
//////////////////////////////////////////////////////////////////////////////////

std::vector<int> utilities::CountIntegerPartitions(const int n,	//!<	Number to be partitioned	
	const int p,		//!<	Maximum number of parts in the partition		
	const int nbrColours)//!<	Number of colours
{
	int temp[n];
	std::vector<int> nbrTerms;

	CountIterator(n,p,temp,0,nbrColours,nbrTerms);

	//std::cout<<"\n\tNUMBER OF PARTITIONS "<<nbrTerms.size()<<std::endl;
	//std::cout<<"\n\n\tMEMORY REQUIRED TO STORE PARTITIONS "
	<<(sizeof(std::vector<std::vector<int> >)*2*totNbrTerms)/(1024.0*1024.0*1024.0)
	<<" GB"<<std::endl;getchar();

	return nbrTerms;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
#if 0

//////////////////////////////////////////////////////////////////////////////////
//!	\brief	Iteratively calling this function generates a the number of terms in
//!	a list of colour associations of nColour different colours into nSlots labels.
//!
//////////////////////////////////////////////////////////////////////////////////

void utilities::ColourCountIterator(
	int nbrColours,	//!<	Number of colours to be used
	const int nSlots,	//!<	Number of slots to label
	int firstSlot,		//!<	The first slot to fillingFactor at each level of iteration
	int level,			//!<	Iteration level (initial call with level=0)
	int* nCalls)		//!<	A memory address to store the call count
{	
	int first;
	int i,j;

	nCalls[0]++;

	//	Go to next level of iteration

	if(nbrColours>1)
	{
		for(i=firstSlot;i<nSlots;i++)
		{
			ColourCountIterator(nbrColours-1,nSlots,i,level+1,nCalls);
		}
	}
	
	return;
}

#endif
////////////////////////////////////////////////////////////////////////////////

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

