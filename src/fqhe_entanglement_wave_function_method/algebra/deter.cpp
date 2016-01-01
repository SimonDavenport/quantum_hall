#include <complex>
#include <iostream>

using namespace std;

typedef std::complex<double> dcmplx;

extern "C"{
           int zgetrf_ ( int* , int* , dcmplx *, int* , int *, int * );
          }



    dcmplx LogDeterminant(dcmplx** M, int dim)    
    {        
        int *num_piv, info_LU, i, sign;    
        dcmplx det, *A;             

 
        //    First perform an LU decomposition of the matrix A



        A=new dcmplx [dim*dim]; 
        num_piv = new int [(int) dim];

        for (int i=0; i<dim; i++)		
            for(int j=0; j<dim; j++)              
                A[dim*j+i]=M[i][j];


       zgetrf_ ( &dim, &dim, A, &dim, num_piv, &info_LU);
 
      



       det=1.0;
       for (i=0;i<dim;++i)
           det*=A[i*dim+i];


  

       sign=1.0;
       for (i=1;i<=dim;++i)
           {
           if(num_piv[i-1]!=i)
             sign=-sign;      
           }


       det*=sign;

       
        delete [] A;
        delete [] num_piv;


        return log(det);

    }
 


