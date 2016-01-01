#include <complex>
typedef std::complex<double> dcmplx;

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
 
    //////////////////////////////////////////////////////////////////////////////////
    //!    \brief This function performs a LU decomposition of the matrix A using the
    //!    Crout algorithm. This template uses type U variables. 
    //!
    //!    Algorithm modified from 
    //!    http://www.mymathlib.com/matrices/linearsystems/crout.html
    //!     
    //!    Note that the input matrix A will be overwritten with the LU output.
    //!     
    //!    Template argument U corresponds to the arbitrary precision type.
    //!
    //!    \return The number of pivoting operators that took place
    //!     
    //////////////////////////////////////////////////////////////////////////////////
 
    template <typename U>
    int CroutLuDecomposition(
        U *A,            //!<    Type U template variable
        const int dim,    //!<    Dimension of matrix         
        int* pivot)        //!<    Pivot vector: stores a list of the row    permutations
                        //!        when pivoting occurs: i-th element is pivot row 
                        //!        interchanged with row i     
    {
        //    Declare local variables     
        int i,j,k;                            //    Generic loop variables
        U *p_k,*p_row,*p_col;                //    Pointers to a row k, or a given row or 
                                            //    column of the LU matrix     
        U Max;                                //    Max absolute value of a LU element - use
                                            //    to check for pivoting
        double Abs,Abs1;                    //    Absolute values of elements of LU rows -
                                            //    use to check for pivoting
        int pivotCount=0;                    //    Counts the number of pivot permutations
                                            //    to keep track of the sign of LU     
                                         
        //    For each row and column, k = 0, ..., dim-1,
 
        for (k = 0, p_k = A; k < dim; p_k += dim, k++) {
 
            //    Find the pivot row
 
            pivot[k] = k;Max=*(p_k + k);
            for (j = k + 1, p_row = p_k + dim; j < dim; j++, p_row += dim)
            {
                Abs=abs(Max);
                Abs1=abs(*(p_row+k));
 
                if(Abs<Abs1)
                {
                    Max=*(p_row + k);pivot[k] = j;p_col = p_row;
                }
            }
 
            //    And if the pivot row differs from the current row, then interchange the two rows.
 
            if (pivot[k] != k){pivotCount++;
                for (j = 0; j < dim; j++) {
                    Max=*(p_k + j);
                    *(p_k + j)=*(p_col + j);
                    *(p_col + j)=Max;
                }
            }
 
            //    Find the upper triangular matrix elements for row k.
 
            for (j = k+1; j < dim; j++) {*(p_k + j)/=*(p_k + k);}
 
            //    Update remaining matrix
 
            for (i = k+1, p_row = p_k + dim; i < dim; p_row += dim, i++)
            {
                for (j = k+1; j < dim; j++)
                {
                    *(p_row + j)-=*(p_row + k)**(p_k + j);
                }
            }
        }
 
        return pivotCount;
    }
 
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//


//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
 
    //////////////////////////////////////////////////////////////////////////////////
    //!    \brief Calculates the log of the determinant of a matrix.
    //!     
    //!    This is done using LU decomposition then taking the trace.
    //!    This is an overload of the LogDeterminant function to U types.
    //!
    //!    Note that the input matrix A will be overwritten with the LU output.
    //!
    //!    \return The (complex) value of the log of determinant
    //!
    //////////////////////////////////////////////////////////////////////////////////
 
    template <typename U>
    dcmplx LogDeterminant(U** M, //!<    Memory address of the stored matrix
        const int dim)    //!<    Dimension of the matrix
    {
        //    Declare local variables     
        int pivotCount;     //  Count the number of pivot permutations of the matrix.
                            //  Each permutation multiplies the determinant by -1.
        int pivot[dim];        //    Pivoting vector     
        U det, *A;                //    To store the value of the determinant
        dcmplx returnVal;    //    Value to return 
        const long double PI=3.1415926535897932384626433832795028841972;
 
        //    First perform an LU decomposition of the matrix A

         A=new U [dim*dim+1]; 

        for (int i=0; i<dim; i++)		
            for(int j=0; j<dim; j++)              
                A[i+dim*j]=M[j][i];


 
        pivotCount=CroutLuDecomposition<U>(A,dim,pivot);
 
        //    Calculate the determinant from the trace
 
        det = 1.0;
 
        for(int i=0;i<dim;i++)
        {
            det*=*(A+i*dim+i);
        }
 

        delete [] A;
        //    Take the log
 
//        det=log(det);
 
        returnVal=det;
 
        //    Include contribution from pivoting steps
 
        returnVal+=dcmplx(0,PI*pivotCount);
//        returnVal+=dcmplx(0,PI*pivotCount);
 
        return returnVal;
    }
 
//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

