////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 11/04/2014
//!
//!  \file 
//!   	Library file for the HilbertSpace class. Used to generate and diagonalize
//!	 	a matrix representation of an entanglement Hamiltonian written in terms
//!	 	of U(1) current operators and Majorana operators.
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

#include "hilbert.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief Default constructor for an object of type HilbertSpace. 
//!
//! Initialize class flags and pointers
//!
//////////////////////////////////////////////////////////////////////////////////

HilbertSpace::HilbertSpace()
    :
    m_matrixRepresentation(0),
    m_blockRepresentation(0),
    m_eigenvalues(0),
    m_eigenvalueLabels(0),
    m_hilbertSpaceBuilt(false),
    m_matrixGenerated(false),
    m_matrixDiagonalized(false)
{
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Destructor for an object of type HilbertSpace. 
//! 
//! Deallocates memory to store matrix representation of the Hamiltonian.
//!	
////////////////////////////////////////////////////////////////////////////////

HilbertSpace::~HilbertSpace()
{
	//	Deallocate memory if allocated

    if(m_eigenvalues!=0)            delete[] m_eigenvalues;
    if(m_eigenvalueLabels!=0)       delete[] m_eigenvalueLabels;
	if(m_matrixRepresentation!=0)	delete[] m_matrixRepresentation;
	if(m_blockRepresentation!=0)
	{
		for(unsigned int i=0;i<m_blockDims.size();i++)
		{
			delete[] m_blockRepresentation[i];
		}
		
		delete[] m_blockRepresentation;
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief A function to construct the Hilbert space
//!
//! In this constructor we generate the full Hilbert space for a given number
//!	of U(1) currents, a given system size and for a given number of orbital
//!	angular momentum sectors. The constructor also allocates memory to store
//!	the matrix representation of the Hamiltonian.
//!	
//////////////////////////////////////////////////////////////////////////////////

void HilbertSpace::BuildHilbertSpace(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzA)                                      //!<    Angular momentum
{
    if(!m_hilbertSpaceBuilt)
    {
	    utilities::cout.MainOutput()<<"\t- GENERATING HILBERT SPACE"<<std::endl;
 
        //  Extract command line variables
        const int nbrColours = (*optionList)["nbr-ll"].as<int>();
        const int nbrA = (*optionList)["nbr-a"].as<int>();
        const int nbr  = (*optionList)["nbr-n"].as<int>();
        m_useBlocks = (*optionList)["blocks"].as<bool>();
             
        //  Generate a list of colour labels for the neutral modes
        
	    char colours[nbrColours];
	    colours[0]=_U1_LABEL_1_;

	    if(nbrColours>1)	colours[1]=_U1_LABEL_2_;

	    if(nbrColours>2)	colours[2]=_U1_LABEL_3_;

	    if(nbrColours>3)
	    {
		    std::cerr<<"\tERROR: HilbertSpace note implemented for > 3 currents."<<std::endl;
		    exit(EXIT_FAILURE);
	    }

        //  Time the construction
	    double timer=clock();

        //  If we are to include majorana modes also then we need to bi-partition
        //  the excitations and consider all possible combinations of charge
        //  and current modes adding up to a given total angular momentum, 
        //  for a given number of particles.
        
        //  Note that in the entanglement wave function picture, the Pfaffian
        //  factor is separate from the Jastrow factors. This means that the 
        //  number of neutral and majorana modes is independent (and both
        //  are restricted by the number nbr-a of particles in the 
        //  particle space cut)
        
        std::vector<int> neutralModePartition;
        const bool useMajoranaModes = (*optionList)["majorana-modes"].as<bool>();
        
        if(useMajoranaModes && nbrColours > 0)
        {
            //  Generate a list of bi-partitions of the angular momentum
            //  sector into majorana and neutral modes

            //  Always include a zero-dimensional partition
            //  (not otherwise generate in the list of partitions)
            
            //  Complication: When Na is even, there is no charge
            //  mode at lzA = 1 (it is there, however, for na odd)
                    
	        //  Explicitly exclude the case of lzA = 1 for the majorana
	        //  modes since this will give no admissible partitions
	        //  The way the code works below, if there are no majorana
	        //  modes then it's assumed that the neutral modes are 
	        //  still possible (which is fine for the majorana mode
	        //  vacuum, but not fine for the lzA = 1 state)
            
            neutralModePartition.reserve(std::max(1,lzA));
            
            for(int aPartition=0;aPartition<=lzA;aPartition++)
            {
                if( (!((nbrA & 1) == 0) && (lzA-1 != aPartition)))
                {
                    neutralModePartition.push_back(aPartition);
                }
            }
        }
        else if(useMajoranaModes && 0 == nbrColours)
        {
            //  If we don't include neutral modes then the list of 
            //  neutral partitions contains only a single term
            
            neutralModePartition.push_back(0);
        }
        else
        {
            //  If we don't include majorana modes then the list of 
            //  neutral partitions contains only a single term
            
            neutralModePartition.push_back(lzA);
        }

	    //  For each bi-partition of angular momentum, generate all the possible
	    //  assignments of the particle cut to the set of LLs for a given sector
	    //	Also store their angular momentum relative to the lowest LzA cut

	    for(unsigned int p=0;p<neutralModePartition.size();p++)
	    {
	        //  First generate the majorana mode partitions, if requested
	        
	        std::vector<std::vector<int> > majoranaPartitionList;
	        
	        if(useMajoranaModes)
	        {
	            //  Generate a complete list of majorana mode partitions
                
                //  For the Pfaffian state this is given by the list 
                //  of (k,r) admissible integer partitions (for k=r=2)
                
                utilities::Generate_k_r_AdmissiblePartitions(majoranaPartitionList,
		        (lzA-neutralModePartition[p]),nbrA,2,2,(nbrA & 1));	
		        
		        //  These partitions do not have colour labels   
		    }
		    
		    const int dimMajorana = majoranaPartitionList.size();
        
	        //  If we have ONLY majorana modes, then we need a separate
	        //  loop to construct the charge mode only Hilbert space
	        
	        if(0 == nbrColours)
	        {
	            if(m_useBlocks)
                {  
                    m_blockDims.push_back(dimMajorana);
                }   

		        for(int j=0;j<dimMajorana;j++)
	            {
	                const int size = majoranaPartitionList[j].size();
	            
	                std::vector<short int> tempInt(size);
	                std::vector<char> tempChar(size);
	            
                    for(int k=0;k<size;k++)
	                {
		                tempInt[k] = (-majoranaPartitionList[j][k]);
	                }

	                for(int k=0;k<size;k++)
	                {
		                tempChar[k] = _MAJORANA_CURRENT_;
	                }

		            //	Calculate normalization factor
		            
		            std::vector<int> dummyNbr(1);
		            std::vector<int> dummNa0(1);
		            
		            dummyNbr[0] = nbr;
		            dummNa0[0]  = 0;

		            Term state(1.0,tempInt,tempChar,dummyNbr,dummNa0);

                    state.SetCoefficient(1.0/sqrt(InnerProduct(state,state)));

                    //  Include this term in the Hilbert space

		            m_hilbertSpace.push_back(state);
		        }
		        
	            //  Skip the construction of the neutral modes   
	            
                break;
	        }

	        //  For the coloured partition case as applied to composite fermions we 
	        //  need to take into account the sets of counting associated with 
	        //  each configuration of the particles in the cut distributed amongst 
	        //  the CF LLs

	        //  This amounts to the construction of a set of identical Hilbert spaces, 
	        //  each starting from a  different angular momentum sector and each being 
	        //  associated with a different colour configuration of the particle cut. 
	
	        //	If possible, we can separate the matrix into a block structure.
	        //	This occurs in cases where there are sectors with different occupation
	        //	numbers. In this case we need to specify a list of the dimensions of 
	        //	each block.
	
	        //	For the 1-current case, there will only be one type of cut (the one 
	        //	specified by nbrA out  of nbrParticles)
	    
	        std::vector<std::vector<int> > deltaNa0List;		//  List of minimum (Na-Na0)
	        std::vector<std::vector<int> > nbrEachLL;	        //  List of number in each LL       
	        std::vector<int> deltaLza;						    //	Difference in Lz compared
														        //	to Na0 cut
	        std::vector<std::vector<int> > maxExcitationSize;   //  Maximum allowed excitation in 
													            //	a given LL
													            
	        //	Note Na0 is defined to be the particle cut corresponding to the minimum 
	        //	pseudo-energy and such that (Na-Na0) is minimal in cases where Na0 is 
	        //  degenerate.
	
	        occupationBasis::GenerateOccupationData(nbr,nbrA,neutralModePartition[p],nbrColours,
	                                 deltaNa0List,nbrEachLL,deltaLza,maxExcitationSize);

            //  We can bloch diagonalize the Hamiltoian at the very least
            //  in terms of the NA0 sectors (since these sectors are 
            //  orthogonal by construction)

            const int nbrBlocks = deltaNa0List.size();
	
	        for(int b=0;b<nbrBlocks;b++)
	        {
		        //utilities::cout.DebuggingInfo()<<"\n\lzA-deltaLza[b] = "<<lzA-deltaLza[b]<<std::endl;
		
		        //utilities::cout.DebuggingInfo()<<"nbr in LL 0: "<<nbrEachLL[b][0]<<std::endl;
                //utilities::cout.DebuggingInfo()<<"nbr in LL 1: "<<nbrEachLL[b][1]<<std::endl;
                //utilities::cout.DebuggingInfo()<<"nbr in LL 2: "<<nbrEachLL[b][2]<<std::endl;
		
		        if(neutralModePartition[p]-deltaLza[b]>0)
		        {
			        //	Generate integer partitions
					
			        std::vector<std::vector<int> > neutralPartitionList;

			        utilities::GenerateColouredIntegerPartitions(neutralPartitionList,neutralModePartition[p]-deltaLza[b],
			        nbrA,maxExcitationSize[b],nbrEachLL[b],nbrColours,colours);

                    //  These partitions have colour labels

			        const int dimNeutral = neutralPartitionList.size();

			        if(m_useBlocks && dimNeutral != 0)
                    {  
                        if(dimMajorana == 0)
                        {
                            m_blockDims.push_back(dimNeutral);
                        }
                        else
                        {
                            m_blockDims.push_back(dimNeutral*dimMajorana);
                        }
                    }   
			
			        //  Add these integer partitions to the Hilbert space
			
			        for(int i=0;i<dimNeutral;i++)
			        {
			            const int size = neutralPartitionList[i].size()/2;
			        
			            std::vector<short int> tempInt(size);
			            std::vector<char> tempChar(size);

			            for(int k=0;k<size;k++)
			            {
				            tempInt[k]=-neutralPartitionList[i][k];
			            }

			            for(int k=size;k<size*2;k++)
			            {
				            tempChar[k-size]=(char)neutralPartitionList[i][k];
			            }
			        
			            for(int j=0;j<dimMajorana;j++)
			            {
			                std::vector<short int> tempInt2(tempInt);
			                std::vector<char> tempChar2(tempChar);
			            
				            const int size2 = majoranaPartitionList[j].size();
	
	                        for(int k=0;k<size2;k++)
			                {
				                tempInt2.push_back(-majoranaPartitionList[j][k]);
			                }

			                for(int k=0;k<size2;k++)
			                {
				                tempChar2.push_back(_MAJORANA_CURRENT_);
			                }
	
				            //	Calculate normalization factor

				            Term state(1.0,tempInt2,tempChar2,nbrEachLL[b],deltaNa0List[b]);

                            state.SetCoefficient(1.0/sqrt(InnerProduct(state,state)));

                            //  Include this term in the Hilbert space

				            m_hilbertSpace.push_back(state);

				            //utilities::cout.DebuggingInfo()<<"\n\nNORM IS "<<matrixElement.FullCommute()<<std::endl;getchar();
				        }
				        
				        //  Or, if there are no majorana modes...
				        
				        if(dimMajorana == 0)
				        {
				            //	Calculate normalization factor

				            Term state(1.0,tempInt,tempChar,nbrEachLL[b],deltaNa0List[b]);

                            state.SetCoefficient(1.0/sqrt(InnerProduct(state,state)));

                            //  Include this term in the Hilbert space

				            m_hilbertSpace.push_back(state);
				        }
			        }
		        }
		        else if(dimMajorana > 0)
		        {
		            //  If there are only majorana modes...

                    if(m_useBlocks)
                    {  
                        m_blockDims.push_back(dimMajorana);
                    }   

			        for(int j=0;j<dimMajorana;j++)
		            {
		                const int size = majoranaPartitionList[j].size();
		            
		                std::vector<short int> tempInt(size);
		                std::vector<char> tempChar(size);
		            
                        for(int k=0;k<size;k++)
		                {

			                tempInt[k] = (-majoranaPartitionList[j][k]);
		                }

		                for(int k=0;k<size;k++)
		                {
			                tempChar[k] = _MAJORANA_CURRENT_;
		                }

			            //	Calculate normalization factor

			            Term state(1.0,tempInt,tempChar,nbrEachLL[b],deltaNa0List[b]);

                        state.SetCoefficient(1.0/sqrt(InnerProduct(state,state)));

                        //  Include this term in the Hilbert space

			            m_hilbertSpace.push_back(state);
			        }  
		        }
		        else if(0 == neutralModePartition[p]-deltaLza[b])
		        {
		            //  If there are no majorana or neutral modes then
			        //	include only the vacuum (an empty Term) in the Hilbert space
			        
			        m_hilbertSpace.push_back(Term(1.0,nbrEachLL[b],deltaNa0List[b]));
			
			        if(m_useBlocks)
                    { 
                        m_blockDims.push_back(1);
                        //utilities::cout.DebuggingInfo()<<"Block Dimension "<<b<<" = "<<m_blockDims[b] <<std::endl;
                    } 
		        }
		        
		        //  Generate a label that uniquely identifies each type of CF branch. 
	            //  For this use the effective cyclotron energy, such that 
	            //  the LCFLL has energy 0 per electron, the 2nd LL has energy 1 per
	            //  electron, etc. 
	            
	            if(m_useBlocks)
                {
	                int blockLabel = 0;

                    for(int c=0;c<nbrColours;c++)
                    {   
                        blockLabel += nbrEachLL[b][c]*(c);
                    }
		            
		            m_blockLabels.push_back(blockLabel);
		        }
	        }
	    }

        //////  PRINT OUT HILBERT SPACE INFORMATION     ////////////////////////
    
	    utilities::cout.SecondaryOutput()<<"\n\t\tHILBERT SPACE CONTENTS: "<<std::endl;

	    for(unsigned int j=0;j<m_hilbertSpace.size();j++)
	    {
		    utilities::cout.SecondaryOutput()<<"\n\t\t";
		    m_hilbertSpace[j].PrintTerm();
		    utilities::cout.SecondaryOutput()<<"|";		
		    m_hilbertSpace[j].PrintDeltaNa0();
		    utilities::cout.SecondaryOutput()<<">"<<std::endl;
	    }
	
	    m_dimension = m_hilbertSpace.size();
		
	    utilities::cout.MainOutput()<<"\n\t\tTIME TAKEN TO GENERATE HILBERT SPACE:\t"
	             <<(clock()-timer)/CLOCKS_PER_SEC<<" seconds.\n"<<std::endl;
	
	    utilities::cout.SecondaryOutput()<<"\t\tHILBERT SPACE DIMENSION:\t"<<m_dimension<<std::endl;

	    int memRequired=0;

	    if(m_useBlocks)
	    {	
		    for(unsigned int i=0;i<m_blockDims.size();i++)
		    {
			    memRequired += sizeof(double)*m_blockDims[i]*m_blockDims[i];
		    }
	    }
	    else
	    {
		    memRequired = sizeof(double)*m_dimension*m_dimension;
	    }

	    utilities::cout.SecondaryOutput()<<"\n\t\tMEMORY REQUIRED TO STORE HAMITLONIAN (ESIMATE)\t"
	             <<(memRequired)/(1024.0*1024.0)<<" MB"<<std::endl;
	
	    //////      ALLOCATE MEMORY TO STORE HAMILTONIAN        ////////////////

	    if(m_useBlocks)
	    {
	        //  If we use blocks then we need to assign the right amount of
            //  memory to store each block

		    m_blockRepresentation = new (std::nothrow) double*[m_blockDims.size()];
		
		    for(unsigned int i=0;i<m_blockDims.size();i++)
		    {
			    m_blockRepresentation[i] = new (std::nothrow) double[m_blockDims[i]*m_blockDims[i]];
		    }
	    }
	    else
	    {
		    m_matrixRepresentation = new (std::nothrow) double[m_dimension*m_dimension];
	    }

	    utilities::cout.SecondaryOutput()<<std::endl;
	
	    m_hilbertSpaceBuilt = true;
	}
	else
	{
	    std::cerr<<"\n\t\tERROR IN BUILD HILBERT SPACE: IT'S ALREADY BUILT!"<<std::endl;
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function converts an operator representation of an entanglement
//!	Hamiltonian into a numerical matrix. 
//!
//! This is done by evaluating matrix
//!	elements of the Hamiltonian with all states in the Hilbert space.
//!	
////////////////////////////////////////////////////////////////////////////////

void HilbertSpace::GenerateMatrixElements(
	ListOfTerms hamiltonian) 	//!<    Operator representation of entanglement	Hamiltonian 
							    //!	 	as a vector of Term objects							
{
    if(m_hilbertSpaceBuilt)
    {
	    //	Calculate matrix elements of the matrix to be diagonalized.

	    utilities::cout.MainOutput()<<"\n\t- POPULATING HAMILTONIAN MATRIX...\n"<<std::endl;

	    double timer = clock();
	
	    if(m_useBlocks)
	    {
	        int totalDimension=0;
	        
	        //	Make a vector H |ket>

	        std::vector<ListOfTerms> ketVector(hamiltonian.OperateOnKet(m_hilbertSpace));

	        //	Populate the full matrix by closing the Hilbert space
	        //	on the ketVector this created. This operation can be made efficient
	        //	by using that zero terms occur when the ketVector is zero. The remaining
	        //	zero and non-zero elements are also evaluated at this stage (simply as
	        //	products of bra and ket vectors)
	        
	        for(unsigned int b=0;b<m_blockDims.size();b++)
	        {
	            //  Populate each block matrix
	            
	            double* p_matrix;
	            double* p_block;

	            int blockDim = m_blockDims[b];
	            bool printMatrix = false;
	            int i;

	            if(blockDim<20)
	            {
		            printMatrix = true;
	            }
	            
	            utilities::cout.MainOutput()<<"\n\t\tBLOCK MATRIX "<<b<<std::endl<<std::endl;
	            
	            p_block = m_blockRepresentation[b];

	            for(p_matrix=p_block,i=0;i<blockDim;i++)
		        {  
			        //Term rowTerm = m_hilbertSpace[i+totalDimension];		//	Used the default copy constructor here
			
			        //rowTerm.Conjugate();
			
			        if(printMatrix)		utilities::cout.SecondaryOutput()<<"\t\t";

			        for(int j=0;j<i;j++,p_matrix++)
			        {
				        *(p_matrix) = p_block[j*blockDim+i];

				        if(printMatrix)		utilities::cout.SecondaryOutput()<<*(p_matrix)<<" ";
			        }
			
			        //	Evaluate matrix elements for upper triangular part only

			        for(int j=i;j<blockDim;j++,p_matrix++)
			        {
				        //ListOfTerms matrixElement(rowTerm,hamiltonian,m_hilbertSpace[j+totalDimension]);

				        //*(p_matrix) = matrixElement.FullCommute();
		
		                *(p_matrix) = InnerProduct(m_hilbertSpace[i+totalDimension],ketVector[j+totalDimension]);
		
				        if(printMatrix)		utilities::cout.SecondaryOutput()<<*(p_matrix)<<" ";

			        }
			
			        if(printMatrix)		utilities::cout.SecondaryOutput()<<std::endl;

			        if(!printMatrix)
			        {
				        utilities::cout.SecondaryOutput()<<std::fixed<<"\tCALCULATED "<<100.0*(i+1)/blockDim<<"%\r";

				        fflush(stdout);
			        }
		        }
		        
		        totalDimension += blockDim;
		        
		        if(!printMatrix)
		        {
		            utilities::cout.SecondaryOutput()<<std::endl;
		        }
	        }   
	    }
	    else
	    {
		    double* p_matrix;
		    int i;
		    bool printMatrix=false;

	        if(m_dimension<25)
	        {
		        printMatrix=true;
	        }
		
		    //	Make a vector H |ket>

		    std::vector<ListOfTerms> ketVector(hamiltonian.OperateOnKet(m_hilbertSpace));

		    //	Populate the full matrix by closing the Hilbert space
		    //	on the ketVector this created. This operation can be made efficient
		    //	by using that zero terms occur when the ketVector is zero. The remaining
		    //	zero and non-zero elements are also evaluated at this stage (simply as
		    //	products of bra and ket vectors)

            //utilities::cout.DebuggingInfo()<<"Ket vector"<<std::endl;

            //for(int i=0;i<ketVector.size();i++)
            //{
            //    ketVector[i].PrintMatrixElement(true);
            //}

		    for(p_matrix=m_matrixRepresentation,i=0;i<m_dimension;i++)
		    {
			    if(printMatrix)		utilities::cout.SecondaryOutput()<<"\t\t";
	
			    for(int j=0;j<i;j++,p_matrix++)
			    {
				    *(p_matrix) = m_matrixRepresentation[j*m_dimension+i];

				    if(printMatrix)		utilities::cout.SecondaryOutput()<<*(p_matrix)<<" ";
			    }
			
			    //	Evaluate matrix elements for upper triangular part only

			    for(int j=i;j<m_dimension;j++,p_matrix++)
			    {
				    //m_hilbertSpace[i].PrintTerm();
				    //utilities::cout.SecondaryOutput()<<" ";
				    //ketVector[j].PrintMatrixElement(false);

				    *(p_matrix) = InnerProduct(m_hilbertSpace[i],ketVector[j]);

				    if(printMatrix)		utilities::cout.SecondaryOutput()<<*(p_matrix)<<" ";

			    }
			
			    if(printMatrix)		utilities::cout.SecondaryOutput()<<std::endl;

			    if(!printMatrix)
			    {
				    utilities::cout.SecondaryOutput()<<std::fixed<<"\tCALCULATED "<<100.0*(i+1)/m_dimension<<"%\r";
				    fflush(stdout);
			    }
		    }
	    }
        utilities::cout.SecondaryOutput()<<std::endl;

	    utilities::cout.MainOutput()<<"\t\tTIME TAKEN TO CONSTRUCT HAMILTONIAN: "<<(clock()-timer)/CLOCKS_PER_SEC<<" seconds.\n"<<std::endl;
        
        m_matrixGenerated = true;
    }
	else
	{
	    std::cerr<<"\n\t\tERROR IN GENERATE MATRIX ELEMENTS: NO HILBERT SPACE BUILT!"<<std::endl;
	}	
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief Calculate the eigenvalues of the entanglement Hamiltonian matrix 
//!	stored in this object
//!
//!	The eigenvalues are then immediately written to a file by the DataToFile 
//!	function
//!
////////////////////////////////////////////////////////////////////////////////

void HilbertSpace::Diagonalize()
{
    if(m_matrixGenerated && !m_matrixDiagonalized)
    {
        //  Allocate memory to store eigenvalues and labels
        
	    m_eigenvalues    = new (std::nothrow) double[m_dimension];
	    m_eigenvalueLabels = new (std::nothrow) int[m_dimension];
	
	    utilities::cout.MainOutput()<<"\t- GENERATING EIGENVALUES..."<<std::endl;
	
	    double timer = clock();
	
	    if(m_useBlocks)
	    {
	        double* p_eigenvalues = m_eigenvalues;
	        int totDim=0;
	
		    for(unsigned int i=0;i<m_blockDims.size();i++)
		    {
			    utilities::linearAlgebra::DiagonalizeSymmetricMatrix<double>(m_blockRepresentation[i],p_eigenvalues,m_blockDims[i],'U');
			
			    p_eigenvalues += m_blockDims[i];
			
			    for(int k=totDim;k<totDim+m_blockDims[i];k++)
	            {
	                m_eigenvalueLabels[k] = m_blockLabels[i];
	                
	                //utilities::cout.DebuggingInfo()<<blockLabelList[k]<<std::endl;
	            }
	            
	            totDim += m_blockDims[i];
		    }
		
		    //  Re-sort eigenvalues and sort the label list to match
		    //  This is important because we need to have the eigenvalues
		    //  in size order

	        utilities::QuickSort<double,int,_ASCENDING_ORDER_>(m_eigenvalues,m_eigenvalueLabels,m_dimension);
	        
	    }
	    else
	    {
		    utilities::linearAlgebra::DiagonalizeSymmetricMatrix<double>(m_matrixRepresentation,m_eigenvalues,m_dimension,'U');
	        
	        for(int k=0;k<m_dimension;k++)
	        {
	            m_eigenvalueLabels[k] = 0;
	        }
	    }

	    //	Print eigenvalues

	    utilities::cout.SecondaryOutput()<<"\n\t\tEIGENVALUES\tBLOCK ID \n"<<std::endl;

	    for(int i=0;i<m_dimension;i++)
	    {
		    utilities::cout.SecondaryOutput()<<std::setw(10)<<"\t"<<m_eigenvalues[i]<<std::setw(2)<<"\t"<<m_eigenvalueLabels[i]<<std::endl;
	    }

	    utilities::cout.MainOutput()<<"\n\t\tTIME TAKEN TO DIAGONALIZE: "<<(clock()-timer)/CLOCKS_PER_SEC<<" seconds.\n"<<std::endl;

        m_matrixDiagonalized = true;
    }
	else
	{
	    std::cerr<<"\n\t\tERROR IN DIAGONALIZE: MATRIX NOT CONSTRUCTED!"<<std::endl;
	}	
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function generates a file name string to identify output files
//!
////////////////////////////////////////////////////////////////////////////////

std::string HilbertSpace::GenFileName(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzA)                                      //!<    Angular momentum
    const
{
	std::stringstream fileName;	//	string stream to build file name

	fileName.str("");

	fileName<<(*optionList)["out-path"].as<std::string>()<<"/eigenvalues_"<<(*optionList)["hamiltonian"].as<std::string>()
	<<"_model_n_"<<(*optionList)["nbr-n"].as<int>()<<"_na_"<<(*optionList)["nbr-a"].as<int>()
	<<"_ll_"<<(*optionList)["nbr-ll"].as<int>()<<"_sector_"<<lzA<<".dat";

    return fileName.str();
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief This function outputs eigenvalue data to a text file
//!
////////////////////////////////////////////////////////////////////////////////

void HilbertSpace::EigenvaluesToFile(
    boost::program_options::variables_map* optionList,  //!<    Parsed command line arguments
    const int lzA)                                      //!<    Angular momentum sector
    const
{
    if(m_matrixDiagonalized)
    {
        std::string fileName = this->GenFileName(optionList,lzA);
        
        std::ofstream f_out;

	    f_out.open(fileName.c_str());

        if(!f_out.is_open())
	    {
		    std::cerr<<"\tFile: "<<fileName<<" not found !!! "<<std::endl<<std::endl;
	    }

	    for(int i=0;i<m_dimension;i++)
	    {
		    f_out<<m_eigenvalues[i]<<"\t";
		    f_out<<m_eigenvalueLabels[i]<<"\n";
	    }

	    f_out.close();
	
	    utilities::cout.MainOutput()<<"\n\t- DATA WRITTEN TO FILE: "<<fileName<<std::endl<<std::endl;
	}
	else
	{
	    std::cerr<<"\n\t\tERROR IN EIGENVLAUES TO FILE: MATRIX NOT DIAGONALIZED!"<<std::endl;
	}

    return;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

