////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!  	Header file for implementations of regular and coloured integer 
//!     partitions. 
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

#ifndef _INTEGER_PARTITIONS_HPP_INCLUDED_
#define _INTEGER_PARTITIONS_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../algorithms/quick_sort.hpp"
#include "../general/cout_tools.hpp"
#include <vector>
#include <algorithm>
#if _DEBUG_
#include "../general/debug.hpp"
#endif

typedef std::vector<std::vector<int> > partition_t;

namespace utilities
{
	void Generate_k_r_AdmissiblePartitions(partition_t& krAdmissible,
        const int n, const int p, const int k, const int r, const bool shift);	
	void GenerateColouredIntegerPartitions(partition_t& partitions,
		const int n, const int p, std::vector<int>& maxPartSize,
		std::vector<int>& maxNbr, const int nbrColours, const char* colours);	
	void GenerateIntegerPartitions(partition_t& partitions,
		const int n, const int p, const int maxPartSize);    
    void GenerateIntegerPartitionsWithPermutations(partition_t& partitions,
		const int n, std::vector<int>& maxPartSize);    
	void PartitionIterator(const int n, const int p, const int maxPartSize, int* v,
	    int level, partition_t& partitionList);        
    void ColouredPartitionIterator(const int n, const int p, std::vector<int>& maxPartSize,
        std::vector<int>& maxNbr, int* v,int level, const int nbrColours, 
        partition_t& partitionList, const char* colours);
	void ColourIterator(int nbrColours, const int nSlots, std::vector<int>& maxNbr,
        int firstSlot, int level, char* v, std::vector<std::vector<char> >& colourList,
	    const char* colours);	
	void UnqiueColourPerms(int* partition, std::vector<std::vector<char> >& colourList, int nSlots);
    
    //!
    //!	Data struct to specify data used for sorting coloured partitions.
    //!	Contains an integer and its associated colour label.
    //!    
	struct PartitionData
	{
	    int integer;
	    char colour;
	};
	bool ColourPermsComp(const PartitionData &lhs, const PartitionData& rhs);
	
	//!
	//!	A class to define an object which is simply a list of PartitionData
	//!	type data structures (to represent a coloured integer partition)
	//!	By representing coloured integer partition data within this class we 
	//!	are able to perform sorting functions on those partitions
	//!	
	class ColouredPartition
	{
        private:
        PartitionData *m_list;
        int m_dim;
        public:
        ColouredPartition(int dimension);
        ~ColouredPartition();
        ColouredPartition(const ColouredPartition& other);
        void Print();
        bool PermutationInList (std::vector<ColouredPartition> searchList);
        void Set(int index,int integer);
        void Set(int index,char colour);
        PartitionData* BeginAddress();
        PartitionData* EndAddress();
	};

    #if 0
	void CountIterator(int n, const int p, int* v, int level, 
	    const int nbrColours, std::vector<int>& nbrTerms);
	std::vector<int> Count_Integer_Partitions(const int n,
		const int p, const int nbrColours);
    void ColourCountIterator(int nbrColours, const int nSlots,
		int firstSlot, int level, int *nCalls);
	#endif
}  //   End namespace utilities
#endif
