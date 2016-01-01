#include "algebra.hpp"				

#define  half 0.5
#define  one  1.0
#define  fpf  5.5
#define  zero 0.0 

typedef std::complex<double> dcmplx;
using namespace std;

//extern "C"{
	//	Define external lapack routine to get the singular value decomposition
//	int zgesvd_(char*, char*, int*, int*, dcmplx*,int*, double*,dcmplx*, int*, dcmplx*, int*,dcmplx*, int*, double*, int*);
  //      int zgetrf_(int*, int*, double*, int*, int*, int* );
    //    int zgeev_(char*, char*, int*, dcmplx* , int* , dcmplx* , dcmplx*, int*, dcmplx*, int*, dcmplx*, int*, double*, int* );
      //    }



/*
 void Det(int dim, dcmplx **A, dcmplx* Det1)
      {
      dcmplx PhaseLU;
      double AbsLU, aux0, PhaseLU1,*AT,mantisa;
      int *IPVT, LDA, INFO, i, j,expo;


     
      LDA=dim;


      IPVT = new int[dim];
      AT = new double[2*(dim+1)*(dim+1)];

//       IPVT = (int*)calloc(dim, sizeof(int));
//       AT = (double*)calloc(2*dim*dim, sizeof(double));    
        


      // to call a Fortran routine from C we have to transform the matrix 

       for (i=0; i<dim; i++)		
            for(j=0; j<dim; j++)              
               {
                AT[2*(j+dim*i)]=real(A[j][i]);
                AT[2*(j+dim*i)+1]=imag(A[j][i]);
                }
      

             		      
      zgetrf_(&LDA, &LDA, AT, &LDA, IPVT, &INFO );
     
      aux0=0;
      for (j=0;j<dim;j++)
          {
           if (IPVT[j] != j+1)  
           aux0=aux0+1;
           }
  
      if (fmod(aux0,2.0)==0)   
       aux0=-1;
       else 
       aux0=1;


        *Det1=1.0;
        for(i=0;i<dim;i++)
           *Det1=(*Det1)*(dcmplx(AT[2*(i+dim*i)],AT[2*(i+dim*i)+1]));

        *Det1=(*Det1)*aux0;


    delete [] IPVT;
    delete [] AT;

         }
    


*/


double factorial(unsigned int n) 
{
    if (n == 0)
       return 1.0;
    return (double)(n) * factorial(n - 1);
}



void binomial (double **binom, double topbinom)

       {  
        int i, j;
 
        binom[0][0]=1;
        binom[1][0]=1;
        binom[1][1]=1;

        for (j=1;j<=topbinom;j++)
        binom[j][0]=1;
                     
             for (i=1;i<=topbinom;i++)
                 {
                  for (j=1;j<=i;j++)
                       binom[i][j]=binom[i-1][j-1] + binom[i-1][j];
                  }      
       }



void rand_init(char *argv[1])
     {
 
      int seed;

//      MersenneTwister mt;

     //////////////////////////////////////////////////////////////////
	
	//	Simon:Initialize mt random number generator
	
	//	Seed the standard random number generator with cluster job ID or system time
	
	int rank = 0;
          
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	if(getenv("PBS_JOBID")==NULL)	
          //seed=time(0);

            seed=time(0) + 10*rank;
	else
            seed=atoi(getenv("PBS_JOBID")) + 10*rank;

     

	srand(seed);
	 
	//	Now seed the Mersenne Twister random number generator using
	//	an array of these standard random numbers
	{
		int INIT_ARRAY_SIZE=200;
		
		long unsigned int initArray[INIT_ARRAY_SIZE];

		for(int i=0;i<INIT_ARRAY_SIZE;i++)
		{
			initArray[i]=rand();
		}

		mt.init_by_array(initArray,INIT_ARRAY_SIZE);
		//	Now the mt random number generator is initialized
		//	generate random doubles in [0,1) with mt.random()
		//	generate random integers with mt.genrand_int31();
	}

	
//////////////////////////////////////////////////////////////////

 
      }      






double nfmod(double a, double b)
{
	//	A description of what this function does?
	return a - b *floor(a/b);
}



double Beta_I(double X, double a, double b)
      {
       return BETAI(a,b,X)*Gamma(a)*Gamma(b)/Gamma(a+b);
       }








void combinationUtil(int arr[], int data[], int start, int end, int index, int r, vector <int> & aux)
{
   
    // Current combination is ready to be printed, print it
    if (index == r)
    {
        for (int j=0; j<r; j++)
            {
         //    printf("%d ", data[j]);
             aux.push_back(data[j]);          
             }
        printf("\n");
        return;
    }
 
    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
    for (int i=start; i<=end && end-i+1 >= r-index; i++)        
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r, aux);
    }
}

 
/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Staring and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */



void Combinations(int n, int r, int** comb, int**comb_comp, int data[])
     {
      int i,j, topbinom, count=0, arr[n], aux1, **arr_aux, count_aux;
      vector <int> aux; 
      double ** binom;

    for (i=0;i<n;i++)
        {arr[i]=i; arr_aux=0;}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// BINOMIAL //////////////////////////////////////////////////

 topbinom=39;

 binom=new double* [int(topbinom+2)];
   for (i=0;i<topbinom+1;i++)
 binom[i]=new double[int(topbinom+2)]();

 binomial (binom, topbinom);

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// INITIALIZE RANDOM NUMBER MERSENNE GENERATOR /////////////////////////////


      combinationUtil(arr, data, 0, n-1, 0, r, aux);

    
      for (i=0; i<binom[n][r]; ++i)
          for (j=0;j<r;++j)
              {
               comb[i][j]=aux[count];
               count+=1;
              }

     arr_aux=new int* [(int)binom[n][r]+2];
       for (i=0;i<binom[n][r]+1;i++)
           arr_aux[i]=new int[n]();



   for (i=0; i<binom[n][r]; ++i)
       for (j=0;j<r;++j)
           arr_aux[i][comb[i][j]]=1;


   for (i=0; i<binom[n][r]; ++i)
       {
        count_aux=0;
        for (j=0;j<n;++j)
            {
             if (arr_aux[i][j]==0)
                {
                 comb_comp[i][count_aux]=j;
                 count_aux+=1;
                 }
            }
        }



      for (i=0;i<binom[n][r]+1;i++)
        delete [] arr_aux[i];
    delete [] arr_aux;


    for (i=0;i<topbinom+1;i++)
         delete [] binom[i];
    delete [] binom;


          
  }
           



