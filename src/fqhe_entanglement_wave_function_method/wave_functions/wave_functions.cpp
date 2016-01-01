#include <iostream>
#include <cstdlib>
#include <math.h>
#include <memory.h>
#include <complex>

#include "../../utilities/mathematics/mt.hpp"
#include "wave_functions.hpp"
#include "../algebra/algebra.hpp"
//#include "../Algebra/determinant.hpp"
#include "../algebra/deter.hpp"


//extern "C"{
	//	Define external lapack routine to get the singular value decomposition
//	int zgesvd_(char*, char*, int*, int*, dcmplx*,int*, double*,dcmplx*, int*, dcmplx*, int*,dcmplx*, int*, double*, int*);
//        int zgetrf_(int*, int*, double*, int*, int*, int* );
//        int zgeev_(char*, char*, int*, dcmplx* , int* , dcmplx* , dcmplx*, int*, dcmplx*, int*, dcmplx*, int*, double*, int* );
//          }


using namespace std;
typedef std::complex<double> dcmplx;


Wave_Functions::Wave_Functions(int ne, int npart, int nu, int LL, char fluxes)
{

int i;

ne_c=ne;
npart_c=npart;
nu_c=nu;
LL_c=LL;
flux=fluxes;

Lz_EWF_aux = new double [3*ne+LL]();

Matrix=new dcmplx* [npart+1];
for (i=0;i<npart+1;i++)
    Matrix[i]=new dcmplx[npart+1]();

IPVT = new int [npart+1]();
AT = new double[2*(npart+1)*(npart+1)]();


int topbinom=39;


 binom=new double* [int(topbinom+2)];
   for (i=0;i<topbinom+1;i++)
        binom[i]=new double[int(topbinom+2)]();

 binomial (binom, topbinom);


}





Wave_Functions::~Wave_Functions()
{

int i, topbinom=39;

//cout<<"destruct"<<endl;
//cin.get();

delete [] Lz_EWF_aux;
delete [] IPVT;
delete [] AT;


//cout<<"ivan"<<endl;
//cin.get();

  for (i=0;i<topbinom+1;i++)
        delete [] binom[i];
   delete [] binom;


for (i=0; i<npart_c; ++i)
     delete [] Matrix[i];
delete [] Matrix;

}




 void Wave_Functions::w_func (dcmplx &norm, double *Lz_EWF, int ne, dcmplx *zu, dcmplx *zv, double **binom, int nu, int LL)
      {
       dcmplx state1, state2, state3;
       double coeff;
       int np_eff, i, j;
       const long double PI=3.1415926535897932384626433832795028841972;


 

  if (flux=='P')            
    {
      if (LL==1 && nu==1)
         {
          CF_EWF_Laugh(state1, Lz_EWF, ne, LL, zu, zv);                              
          norm=log(exp(state1)*Lz_EWF[(int) (LL+Lz_EWF[0])]);  /// REMEMBER THAT SIGN=LZ_EWF[LL+LZ_EWF[0]]    
           }
      else if (LL==1 && nu==2)
              {               

 
               Lz_EWF_aux[0]=Lz_EWF[0];              
               for (j=0;j<Lz_EWF[0];++j)
                    Lz_EWF_aux[j+LL]=Lz_EWF[j+LL];       

               CF_EWF_Laugh(state1, Lz_EWF_aux, ne, LL, zu, zv);
           
               for (j=0;j<Lz_EWF[0];++j)
                    Lz_EWF_aux[j+LL]=Lz_EWF[(int)(j+Lz_EWF[0]+LL)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+nu*Lz_EWF[0])];
          
               norm=(state1+state2) + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
         
               }
      else if (LL==1 && nu==3)
              {               
               Lz_EWF_aux[0]=Lz_EWF[0];              
               for (j=0;j<Lz_EWF[0];++j)
                    Lz_EWF_aux[j+LL]=Lz_EWF[j+LL];
               CF_EWF_Laugh(state1, Lz_EWF_aux, ne, LL, zu, zv);

               for (j=0;j<Lz_EWF[0];++j)
                    Lz_EWF_aux[j+LL]=Lz_EWF[(int)(j+Lz_EWF[0]+LL)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               for (j=0;j<Lz_EWF[0];++j)
                    Lz_EWF_aux[j+LL]=Lz_EWF[(int)(j+2*Lz_EWF[0]+LL)];
               CF_EWF_Laugh(state3, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+nu*Lz_EWF[0])];

               norm=state1 + state2 + state3 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
               }

      else if (LL==2 && nu==0)
              {

               CF_EWF_Jain_UnProj_2_LL(state1, Lz_EWF, ne, LL, zu, zv, binom);
       
               coeff=Lz_EWF[(int) (LL+Lz_EWF[0]+Lz_EWF[1])];
               norm=state1 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));                

               }
      else if (LL==2 && nu==1)
              {       
               Lz_EWF_aux[0]=Lz_EWF[0];              
               Lz_EWF_aux[1]=Lz_EWF[1];
               np_eff=Lz_EWF[0]+Lz_EWF[1];

               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+LL]=Lz_EWF[i+LL];

               CF_EWF_Jain_Proj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);
//               CF_EWF_Jain_UnProj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);

               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+np_eff)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+(nu+1)*np_eff)];

               norm = state1 + state2 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
           
 
               }
      else if (LL==2 && nu==2)
              {
               Lz_EWF_aux[0]=Lz_EWF[0];              
               Lz_EWF_aux[1]=Lz_EWF[1];
               np_eff=Lz_EWF[0]+Lz_EWF[1];


               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+LL]=Lz_EWF[i+LL];

               CF_EWF_Jain_Proj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);
//               CF_EWF_Jain_UnProj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);

               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+np_eff)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+2*np_eff)];
               CF_EWF_Laugh(state3, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+(nu+1)*np_eff)];

               norm = state1 + state2 + state3 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
 
               }
      else if (LL==3 && nu==0)
              {
               //CF_EWF_Jain_Proj_3_LL(state1, Lz_EWF, ne, LL, zu, zv, binom);
               CF_EWF_Jain_UnProj_3_LL(state1, Lz_EWF, ne, LL, zu, zv, binom);

               coeff=Lz_EWF[(int) (LL+Lz_EWF[0]+Lz_EWF[1]+Lz_EWF[2])];
               norm = state1 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
               }
      else if (LL==3 && nu==1)
              {

               
               Lz_EWF_aux[0]=Lz_EWF[0];              
               Lz_EWF_aux[1]=Lz_EWF[1];
               Lz_EWF_aux[2]=Lz_EWF[2];
               np_eff=Lz_EWF[0]+Lz_EWF[1]+Lz_EWF[2];

               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+LL]=Lz_EWF[i+LL];
                   CF_EWF_Jain_Proj_3_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);


               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+np_eff)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+(nu+1)*np_eff)];               
               norm= state1 + state2 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         


               }

      else if (LL==3 && nu==2)
              {
               Lz_EWF_aux[0]=Lz_EWF[0];              
               Lz_EWF_aux[1]=Lz_EWF[1];
               Lz_EWF_aux[2]=Lz_EWF[2];
               np_eff=Lz_EWF[0]+Lz_EWF[1]+Lz_EWF[2];


               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+LL]=Lz_EWF[i+LL];
                   CF_EWF_Jain_Proj_3_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);


               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+np_eff)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+2*np_eff)];
               CF_EWF_Laugh(state3, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+(nu+1)*np_eff)];
               norm=state1+state2+state3+dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
 
               }

       }

  else if (flux=='N')
         {

           if (LL==2 && nu==1)
              {
               Lz_EWF_aux[0]=Lz_EWF[0];              
               Lz_EWF_aux[1]=Lz_EWF[1];
               np_eff=Lz_EWF[0]+Lz_EWF[1];

               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+LL]=Lz_EWF[i+LL];
                   CF_EWF_Jain_negative_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);
//                   CF_EWF_Jain_Proj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);


//               CF_EWF_Jain_UnProj_2_LL(state1, Lz_EWF_aux, ne, LL, zu, zv, binom);

               Lz_EWF_aux[0]=np_eff;
               for (i=0;i<np_eff;++i)
                   Lz_EWF_aux[i+1]=Lz_EWF[(int) (i+LL+np_eff)];
               CF_EWF_Laugh(state2, Lz_EWF_aux, ne, LL, zu, zv);

               coeff=Lz_EWF[(int) (LL+(nu+1)*np_eff)];

               norm = state1 + state2 + dcmplx(0.0, PI/2.0 -(PI/2.0)*coeff/abs(coeff)) + log(abs(coeff));         
           
 
               }

          }



       }




 





/////////////////////////// ROUTINES COMPUTING DIFERENT FQH STATES ///////////////////////////////////////////////





 void Wave_Functions::CF_EWF_Jain_Proj_3_LL (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv, double **binom)
      {

       int i, j, k, dim_M1, n1, n2, n3, LLd;
       double QJain, nee, pp=1.0;
       dcmplx state, FJ[5][5][60],PJ[5][5][60], *zu_eps, *zv_eps;


       nee=ne;
       LLd=LL;
       QJain=(nee/LLd-3.0)/2.0;
       n1=Lz_EWF[0];
       n2=Lz_EWF[1];
       n3=Lz_EWF[2];

       dim_M1=n1+n2+n3;


    


     ///////////////////  LLL PROJECTIONS FROM JAIN'S BOOK //////////////////////      

      for (i=0;i<n1+n2+n3;i++)
          for (k=0;k<n1+n2+n3;k++)
              {
               if (k!=i)
                  {
                   FJ[1][0][i]+= (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k]));
                   FJ[0][1][i]+= (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])); 

                   FJ[2][0][i]+= pow((zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0)*pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
                   FJ[1][1][i]+= pow((zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0)*pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);
                   FJ[0][2][i]+= pow((zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0)*pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0);

                   }
               }
          

      for (k=0; k<n1+n2+n3; k++)
          { 
           PJ[1][0][k]=FJ[1][0][k];
           PJ[0][1][k]=FJ[0][1][k];


           PJ[1][1][k]=pow(pp,2.0)*FJ[0][1][k]*FJ[1][0][k] - pp*FJ[1][1][k];
           PJ[2][0][k]=pow(pp,2.0)*pow(FJ[1][0][k],2) - pp*FJ[2][0][k];
           PJ[0][2][k]=pow(pp,2.0)*pow(FJ[0][1][k],2) - pp*FJ[0][2][k];
           }


   //////////////////////////////////////////////////////////////////////////////




       for (i=0;i<dim_M1;++i)
           {

           for(j=0;j<n1;++j)              
               Matrix[i][j]= pow(zu[i],(int)(QJain+Lz_EWF[j+LL]))*pow(zv[i],(int)(QJain-Lz_EWF[j+LL]));
               //               cout<<Matrix[i][j]<<endl;}
           
           for(j=n1;j<n1+n2;++j)              
               Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*
                             ( binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL]+1)]*zv[i]*PJ[0][1][i] -
                               binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL])]*zu[i]*PJ[1][0][i] );
                            //  cout<<Matrix[i][j]<<endl;}

           for(j=n1+n2;j<n1+n2+n3;++j)              
              Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*
                            ( binom[(int)(2*QJain+2)][(int)(QJain+2-Lz_EWF[j+LL])]*pow(zv[i],2.0)*PJ[0][2][i] -
                          2.0*binom[(int)(2*QJain+2)][(int)(QJain+1-Lz_EWF[j+LL])]*zu[i]*zv[i]*PJ[1][1][i] +
                              binom[(int)(2*QJain+2)][(int)(QJain-Lz_EWF[j+LL])]*pow(zu[i],2.0)*PJ[2][0][i] );
              // cout<<Matrix[i][j]<<endl;}

        
           }



//cin.get();

        norm=LogDeterminant(Matrix,dim_M1);
//      norm=LogDeterminant<dcmplx>(Matrix,dim_M1);
   //        Det1(dim_M1, Matrix, &state);
     //      norm=state;




       }







 void Wave_Functions::CF_EWF_Jain_UnProj_3_LL (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv, double **binom)
      {

       int i, j, k, dim_M1, n1, n2, n3;
       double QJain, nee, pp=1.0, LLd;
       dcmplx state, FJ[5][5][50],PJ[5][5][50], *zu_eps, *zv_eps;


       nee=ne;
       LLd=LL;
       QJain=(nee/LLd-3.0)/2.0;
       n1=Lz_EWF[0];
       n2=Lz_EWF[1];
       n3=Lz_EWF[2];

       dim_M1=n1+n2+n3;


  
       for (i=0;i<dim_M1;++i)
           {

           for(j=0;j<n1;++j)              
               Matrix[i][j]= pow(zu[i],(int)(QJain+Lz_EWF[j+LL]))*pow(zv[i],(int)(QJain-Lz_EWF[j+LL]));
                
           
           for(j=n1;j<n1+n2;++j)              
               Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*
                             ( binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL]+1)]*pow(abs(zv[i]),2.0)-
                               binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL])]*pow(abs(zu[i]),2.0) );
             

           for(j=n1+n2;j<n1+n2+n3;++j)              
               Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*
                            ( binom[(int)(2*QJain+2)][(int)(QJain+2-Lz_EWF[j+LL])]*pow(abs(zv[i]),4.0) -
                          2.0*binom[(int)(2*QJain+2)][(int)(QJain+1-Lz_EWF[j+LL])]*pow(abs(zv[i]),2.0)*pow(abs(zu[i]),2.0) +
                              binom[(int)(2*QJain+2)][(int)(QJain-Lz_EWF[j+LL])]*pow(abs(zu[i]),4.0) );





        
           }


             norm=LogDeterminant(Matrix,dim_M1);
            // norm=LogDeterminant<dcmplx>(Matrix,dim_M1);

   //        Det1(dim_M1, Matrix, &state);
   //        norm=state;


       }


















 void Wave_Functions::CF_EWF_Jain_Proj_2_LL (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv, double **binom)
      {
       int i, j, k, dim_M1, n1, n2;
       double QJain, nee;
       dcmplx  state, FJ[5][5][60],PJ[5][5][60], *zu_eps, *zv_eps;
    
       nee=ne;
       QJain=(nee/2.0-2.0)/2.0;
       n1=Lz_EWF[0];
       n2=Lz_EWF[1];

       dim_M1=n1+n2;



     ///////////////////  LLL PROJECTIONS FROM JAIN'S BOOK //////////////////////      

      for (i=0;i<n1+n2;i++)
          for (k=0;k<n1+n2;k++)
              {
               if (k!=i)
                  {
                   FJ[1][0][i] += (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k]));
                   FJ[0][1][i] += (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])); 
                   }
               }
          


      for (k=0; k<n1+n2;k++)
          { 
           PJ[1][0][k]=FJ[1][0][k];
           PJ[0][1][k]=FJ[0][1][k];
           }


   //////////////////////////////////////////////////////////////////////////////

 
       for (i=0;i<dim_M1;++i)
           {
           for(j=0;j<n1;++j)              
               Matrix[i][j]= pow(zu[i],(int)(QJain+Lz_EWF[j+LL]))*pow(zv[i],(int)(QJain-Lz_EWF[j+LL]));
               
               
           
           for(j=n1;j<n1+n2;++j)              
               Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*
                             ( binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL]+1)]*zv[i]*PJ[0][1][i] -
                               binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL])]*zu[i]*PJ[1][0][i] );

        
           }


               norm=LogDeterminant(Matrix,dim_M1);


       //      norm=LogDeterminant<dcmplx>(Matrix,dim_M1);
        
     //      Det1(dim_M1, Matrix, &state);
     //      norm=state;

           
        }













 void Wave_Functions::CF_EWF_Jain_UnProj_2_LL (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv, double **binom)
      {
       int i, j, k, dim_M1, n1, n2;
       double QJain, nee;
       dcmplx  state;         
    
       nee=ne;
       QJain=(nee/2.0-2.0)/2.0;
       n1=Lz_EWF[0];
       n2=Lz_EWF[1];

       dim_M1=n1+n2;

  

  
       for (i=0;i<dim_M1;++i)
           {
           for(j=0;j<n1;++j)              
               Matrix[i][j]= pow(zu[i],(int)(QJain+Lz_EWF[j+LL]))*pow(zv[i],(int)(QJain-Lz_EWF[j+LL]));
              

           
           for(j=n1;j<n1+n2;++j)                            
               Matrix[i][j]= pow(zu[i],(int) (QJain+Lz_EWF[j+LL]))*pow(zv[i],(int) (QJain-Lz_EWF[j+LL]))*              
                             ( binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL]+1)]*pow((double)abs(zv[i]),(int)2) -
                               binom[(int)(2*QJain+1)][(int) (QJain-Lz_EWF[j+LL])]*pow((double)abs(zu[i]),(int)2) );
  
        
               
           }


/*          cout<<"{";
      for (i=0; i<dim_M1; ++i)
          {
          cout<<"{";
          for (j=0; j<n1+n2-1;++j)
             cout<<real(Matrix[i][j])<<"+I*"<<imag(Matrix[i][j])<<",";
          for (j=n1+n2-1; j<n1+n2;++j)
             cout<<real(Matrix[i][j])<<"+I*"<<imag(Matrix[i][j]);
          cout<<"},"<<endl;
          }
          cout<<"}";
      cin.get(); */
 
//             norm=LogDeterminant<dcmplx>(Matrix,dim_M1);
//       cout<<"det1:  "<<norm<<endl;

           norm=LogDeterminant(Matrix,dim_M1);
 //      cout<<"det2:  "<<norm<<endl;
   //    cout<<"Log(det2):  "<<log(norm)<<endl;
  //    cin.get();

//           Det1(dim_M1, Matrix, &state); 
//           norm=state;     
           
        }















void Wave_Functions::CF_EWF_Jain_negative_2_LL (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv, double **binom)
     {

     int i,j,k, dim_M1, n1, n2;
     dcmplx FJ[20][20][20],PJ[20][20][20],laugh,expo,wfuncCF,wfuncCF1,wfuncCF2;
     double p,QJain,QLaugh;


 p=1.0;

 n1=Lz_EWF[0];
 n2=Lz_EWF[1];
 dim_M1=n1+n2;




// FJ[20][20][20]={0};
// PJ[20][20][20]={0};

      for (i=0;i<20;i++)
          for (j=0;j<20;j++)
              for (k=0;k<20;k++)
                  {FJ[i][j][k]=dcmplx(0.0,0.0);PJ[i][j][k]=dcmplx(0.0,0.0);}






 for (i=0;i<dim_M1;i++)
     for (k=0;k<dim_M1;k++)
         {
          if (k != i) 
             {
              FJ[1][0][i]= FJ[1][0][i] + pow( (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
              FJ[0][1][i]= FJ[0][1][i] + pow( (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);


              FJ[2][0][i]= FJ[2][0][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
              FJ[1][1][i]= FJ[1][1][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);
              FJ[0][2][i]= FJ[0][2][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0);


              FJ[3][0][i]= FJ[3][0][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
              FJ[2][1][i]= FJ[2][1][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);
              FJ[1][2][i]= FJ[1][2][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0);
              FJ[0][3][i]= FJ[0][3][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0);


              FJ[4][0][i]= FJ[4][0][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),4.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
              FJ[3][1][i]= FJ[3][1][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);
              FJ[2][2][i]= FJ[2][2][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0);
              FJ[1][3][i]= FJ[1][3][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0);
              FJ[0][4][i]= FJ[0][4][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),4.0);


              FJ[5][0][i]= FJ[5][0][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),5.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0);
              FJ[4][1][i]= FJ[4][1][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),4.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0);
              FJ[3][2][i]= FJ[3][2][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0);
              FJ[2][3][i]= FJ[2][3][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),2.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),3.0);
              FJ[1][4][i]= FJ[1][4][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),1.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),4.0);
              FJ[0][5][i]= FJ[0][5][i] + pow(  (zv[k]/(zu[i]*zv[k]-zv[i]*zu[k])),0.0) *  pow( (-zu[k]/(zu[i]*zv[k]-zv[i]*zu[k])),5.0); 
              }
         }






 
  for (j=0;j<dim_M1;j++)
      {

   PJ[1][0][j]=p*FJ[1][0][j];
   PJ[0][1][j]=p*FJ[0][1][j];

   PJ[1][1][j]=pow(p,2)*FJ[0][1][j]*FJ[1][0][j] - p*FJ[1][1][j];
   PJ[2][0][j]=pow(p,2)*pow(FJ[1][0][j],2) - p*FJ[2][0][j];
   PJ[0][2][j]=pow(p,2)*pow(FJ[0][1][j],2) - p*FJ[0][2][j];

   PJ[3][0][j]=pow(p,3)*pow(FJ[1][0][j],3) - 3*pow(p,2)*FJ[1][0][j]*FJ[2][0][j] + 
               2*p*FJ[3][0][j];
   PJ[0][3][j]=pow(p,3)*pow(FJ[0][1][j],3) - 3*pow(p,2)*FJ[0][1][j]*FJ[0][2][j] + 
               2*p*FJ[0][3][j];
   PJ[2][1][j]=pow(p,3)*FJ[0][1][j]*pow(FJ[1][0][j],2) - 2*pow(p,2)*FJ[1][0][j]*FJ[1][1][j] - 
               pow(p,2)*FJ[0][1][j]*FJ[2][0][j] + 2*p*FJ[2][1][j];
   PJ[1][2][j]=pow(p,3)*pow(FJ[0][1][j],2)*FJ[1][0][j] - pow(p,2)*FJ[0][2][j]*FJ[1][0][j] - 
               2*pow(p,2)*FJ[0][1][j]*FJ[1][1][j] + 2*p*FJ[1][2][j];
  
   PJ[4][0][j]= pow(p,4)*pow(FJ[1][0][j],4) - 6*pow(p,3)*pow(FJ[1][0][j],2)*FJ[2][0][j] + 
                3*pow(p,2)*pow(FJ[2][0][j],2) + 8*pow(p,2)*FJ[1][0][j]*FJ[3][0][j] - 6*p*FJ[4][0][j];
   PJ[0][4][j]= pow(p,4)*pow(FJ[0][1][j],4) - 6*pow(p,3)*pow(FJ[0][1][j],2)*FJ[0][2][j] + 
                3*pow(p,2)*pow(FJ[0][2][j],2) + 8*pow(p,2)*FJ[0][1][j]*FJ[0][3][j] - 6*p*FJ[0][4][j];
   PJ[1][3][j]= pow(p,4)*pow(FJ[0][1][j],3)*FJ[1][0][j] - 3*pow(p,3)*FJ[0][1][j]*FJ[0][2][j]*
                FJ[1][0][j] + 2*pow(p,2)*FJ[0][3][j]*FJ[1][0][j] - 3*pow(p,3)*pow(FJ[0][1][j],2)*
                FJ[1][1][j] + 3*pow(p,2)*FJ[0][2][j]*FJ[1][1][j] + 6*pow(p,2)*FJ[0][1][j]*
                FJ[1][2][j] - 6*p*FJ[1][3][j];
   PJ[3][1][j]= pow(p,4)*FJ[0][1][j]*pow(FJ[1][0][j],3) - 3*pow(p,3)*pow(FJ[1][0][j],2)*FJ[1][1][j] - 
                3*pow(p,3)*FJ[0][1][j]*FJ[1][0][j]*FJ[2][0][j] + 3*pow(p,2)*FJ[1][1][j]*FJ[2][0][j] + 
                6*pow(p,2)*FJ[1][0][j]*FJ[2][1][j] + 2*pow(p,2)*FJ[0][1][j]*FJ[3][0][j] - 6*p*FJ[3][1][j];
   PJ[2][2][j]= pow(p,4)*pow(FJ[0][1][j],2)*pow(FJ[1][0][j],2) - pow(p,3)*FJ[0][2][j]*pow(FJ[1][0][j],2) - 
                4*pow(p,3)*FJ[0][1][j]*FJ[1][0][j]*FJ[1][1][j] + 2*pow(p,2)*pow(FJ[1][1][j],2) + 
                4*pow(p,2)*FJ[1][0][j]*FJ[1][2][j] - pow(p,3)*pow(FJ[0][1][j],2)*FJ[2][0][j] + 
                pow(p,2)*FJ[0][2][j]*FJ[2][0][j] + 4*pow(p,2)*FJ[0][1][j]*FJ[2][1][j] - 6*p*FJ[2][2][j];

   PJ[5][0][j]= pow(p,5)*pow(FJ[1][0][j],5) - 10*pow(p,4)*pow(FJ[1][0][j],3)*FJ[2][0][j] + 
                15*pow(p,3)*FJ[1][0][j]*pow(FJ[2][0][j],2) + 20*pow(p,3)*pow(FJ[1][0][j],2)*FJ[3][0][j] - 
                20*pow(p,2)*FJ[2][0][j]*FJ[3][0][j] - 30*pow(p,2)*FJ[1][0][j]*FJ[4][0][j] + 24*p*FJ[5][0][j];
   PJ[0][5][j]= pow(p,5)*pow(FJ[0][1][j],5) - 10*pow(p,4)*pow(FJ[0][1][j],3)*FJ[0][2][j] + 
                15*pow(p,3)*FJ[0][1][j]*pow(FJ[0][2][j],2) + 20*pow(p,3)*pow(FJ[0][1][j],2)*FJ[0][3][j] - 
                20*pow(p,2)*FJ[0][2][j]*FJ[0][3][j] - 30*pow(p,2)*FJ[0][1][j]*FJ[0][4][j] + 24*p*FJ[0][5][j];
   PJ[4][1][j]= pow(p,5)*FJ[0][1][j]*pow(FJ[1][0][j],4) - 4*pow(p,4)*pow(FJ[1][0][j],3)*FJ[1][1][j] - 
                6*pow(p,4)*FJ[0][1][j]*pow(FJ[1][0][j],2)*FJ[2][0][j] + 12*pow(p,3)*FJ[1][0][j]*
                FJ[1][1][j]*FJ[2][0][j] + 3*pow(p,3)*FJ[0][1][j]*pow(FJ[2][0][j],2) + 
                12*pow(p,3)*pow(FJ[1][0][j],2)*FJ[2][1][j] - 12*pow(p,2)*FJ[2][0][j]*FJ[2][1][j] + 
                8*pow(p,3)*FJ[0][1][j]*FJ[1][0][j]*FJ[3][0][j] - 8*pow(p,2)*FJ[1][1][j]*FJ[3][0][j] - 
                24*pow(p,2)*FJ[1][0][j]*FJ[3][1][j] - 6*pow(p,2)*FJ[0][1][j]*FJ[4][0][j] + 24*p*FJ[4][1][j];
   PJ[1][4][j]= pow(p,5)*pow(FJ[0][1][j],4)*FJ[1][0][j] - 6*pow(p,4)*pow(FJ[0][1][j],2)*FJ[0][2][j]*
                FJ[1][0][j] + 3*pow(p,3)*pow(FJ[0][2][j],2)*FJ[1][0][j] + 8*pow(p,3)*FJ[0][1][j]*
                FJ[0][3][j]*FJ[1][0][j] - 6*pow(p,2)*FJ[0][4][j]*FJ[1][0][j] - 
                4*pow(p,4)*pow(FJ[0][1][j],3)*FJ[1][1][j] + 12*pow(p,3)*FJ[0][1][j]*FJ[0][2][j]*FJ[1][1][j] - 
                8*pow(p,2)*FJ[0][3][j]*FJ[1][1][j] + 12*pow(p,3)*pow(FJ[0][1][j],2)*FJ[1][2][j] - 
                12*pow(p,2)*FJ[0][2][j]*FJ[1][2][j] - 24*pow(p,2)*FJ[0][1][j]*FJ[1][3][j] + 24*p*FJ[1][4][j];
   PJ[2][3][j]= pow(p,5)*pow(FJ[0][1][j],3)*pow(FJ[1][0][j],2) - 3*pow(p,4)*FJ[0][1][j]*FJ[0][2][j]*
                pow(FJ[1][0][j],2) + 2*pow(p,3)*FJ[0][3][j]*pow(FJ[1][0][j],2) - 
                6*pow(p,4)*pow(FJ[0][1][j],2)*FJ[1][0][j]*FJ[1][1][j] + 6*pow(p,3)*FJ[0][2][j]*
                FJ[1][0][j]*FJ[1][1][j] + 6*pow(p,3)*FJ[0][1][j]*pow(FJ[1][1][j],2) + 
                12*pow(p,3)*FJ[0][1][j]*FJ[1][0][j]*FJ[1][2][j] - 12*pow(p,2)*FJ[1][1][j]*FJ[1][2][j] - 
                12*pow(p,2)*FJ[1][0][j]*FJ[1][3][j] - pow(p,4)*pow(FJ[0][1][j],3)*FJ[2][0][j] + 
                3*pow(p,3)*FJ[0][1][j]*FJ[0][2][j]*FJ[2][0][j] - 2*pow(p,2)*FJ[0][3][j]*FJ[2][0][j] + 
                6*pow(p,3)*pow(FJ[0][1][j],2)*FJ[2][1][j] - 6*pow(p,2)*FJ[0][2][j]*FJ[2][1][j] - 
                18*pow(p,2)*FJ[0][1][j]*FJ[2][2][j] + 24*p*FJ[2][3][j];
   PJ[3][2][j]= pow(p,5)*pow(FJ[0][1][j],2)*pow(FJ[1][0][j],3) - pow(p,4)*FJ[0][2][j]*pow(FJ[1][0][j],3) - 
                6*pow(p,4)*FJ[0][1][j]*pow(FJ[1][0][j],2)*FJ[1][1][j] + 6*pow(p,3)*FJ[1][0][j]*
                pow(FJ[1][1][j],2) + 6*pow(p,3)*pow(FJ[1][0][j],2)*FJ[1][2][j] - 
                3*pow(p,4)*pow(FJ[0][1][j],2)*FJ[1][0][j]*FJ[2][0][j] + 3*pow(p,3)*FJ[0][2][j]*
                FJ[1][0][j]*FJ[2][0][j] + 6*pow(p,3)*FJ[0][1][j]*FJ[1][1][j]*FJ[2][0][j] - 
                6*pow(p,2)*FJ[1][2][j]*FJ[2][0][j] + 12*pow(p,3)*FJ[0][1][j]*FJ[1][0][j]*FJ[2][1][j] - 
                12*pow(p,2)*FJ[1][1][j]*FJ[2][1][j] - 18*pow(p,2)*FJ[1][0][j]*FJ[2][2][j] + 
                2*pow(p,3)*pow(FJ[0][1][j],2)*FJ[3][0][j] - 2*pow(p,2)*FJ[0][2][j]*FJ[3][0][j] - 
                12*pow(p,2)*FJ[0][1][j]*FJ[3][1][j] + 24*p*FJ[3][2][j];     

     }







       QJain=(ne/2.0 - 2.0)/2.0;




        for (i=0;i<n1;i++)
            for (j=0;j<dim_M1;j++)
                 Matrix[i][j]= PJ[(int) (QJain+Lz_EWF[i+2])][(int) (QJain-Lz_EWF[i+2])][j];
   
 

                    
        for (i=n1;i<dim_M1;i++)
            for (j=0;j<dim_M1;j++)
                 {
                  Matrix[i][j]= binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF[i+2])]*zv[j]*
                                PJ[(int) (QJain+Lz_EWF[i+2])][(int) (QJain-Lz_EWF[i+2]+1)][j] - 
                                    binom[(int) (2.0*QJain+1.0)][(int) (QJain-Lz_EWF[i+2])]*zu[j]*
                                PJ[(int) (QJain+Lz_EWF[i+2]+1)][(int) (QJain-Lz_EWF[i+2])][j];    

    //             cout<<(int) (QJain+Lz_EWF[i+2]) <<" "<<(int) (QJain-Lz_EWF[i+2]+1)<<endl;
      //           cout<<binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF[i+2])]<<endl;
             //    cin.get();
     
                 } 






               norm=LogDeterminant(Matrix,dim_M1);
         //    norm=LogDeterminant<dcmplx>(Matrix,dim_M1);




      }























































































 void Wave_Functions::CF_EWF_Laugh (dcmplx &norm, double *Lz_EWF, int ne, int LL, dcmplx *zu, dcmplx *zv)
      {
       int i, j, k, dim_M1;
       double QLaugh;
       dcmplx state;         
    
       QLaugh=(ne-1.0)/2.0;

       dim_M1=Lz_EWF[0];
 



       
       for (i=0;i<dim_M1;++i)
           for(j=0;j<dim_M1;++j)
               Matrix[i][j]=pow(zu[j],(int)(QLaugh+Lz_EWF[i+1]))*pow(zv[j],(int)(QLaugh-Lz_EWF[i+1]));
              //  cout<<"www: "<<Matrix[i][j]<<endl;}
             
           

             norm=LogDeterminant(Matrix,dim_M1);
//      norm=LogDeterminant<dcmplx>(Matrix,dim_M1);
 

//      cout<<"llll: "<<norm<<endl;
//      cin.get();

//           Det1(dim_M1, Matrix, &state);
 
  //         norm=state;

           
        }














//////////////////////////////// METROPOLIS WAVE FUNCTIONS ////////////////////////////////






void Wave_Functions::Laugh_Test(dcmplx &norm, int ne, int LL, dcmplx *zu, dcmplx *zv, int nu, int npart, char AorB)
     {
      int i,j;
      double QLaugh;
      dcmplx expo, shift;

      QLaugh=(ne-1.0)/2.0;

      norm=1.0;

   
      for (i=0;i<npart;++i)
          for (j=i+1;j<npart;++j)
              norm=norm*(zu[i]/zv[i]-zu[j]/zv[j]);


      expo=1.0;


      for (i=0;i<npart;++i)
          expo=expo*zv[i];

      expo=pow(expo, 2.0*QLaugh);

               
      if (AorB=='A')
          norm=pow(norm*expo,1.0*nu);
      else if (AorB=='B')
              {
               shift=1.0;
               for (i=0;i<npart;++i)
                    shift=shift*(zu[i]/zv[i]);   
               shift=pow(shift,(double)(ne-npart));           
               norm=pow(norm*expo*shift,1.0*nu);
               }



     }




void Wave_Functions::Laugh_Test_A(dcmplx &norm, int ne, int LL, dcmplx *zu, dcmplx *zv, int nu, int npart)
     {
      int i,j;
      double QLaugh;
      dcmplx expo;

      QLaugh=(ne-1.0)/2.0;

   
      norm=0.0;
      for (i=0;i<npart;++i)
          for (j=i+1;j<npart;++j)
              norm+=log(zu[i]/zv[i]-zu[j]/zv[j]);


      expo=0.0;
      for (i=0;i<npart;++i)
          expo+=2.0*QLaugh*log(zv[i]);

            
      norm=double(nu)*(norm+expo);

//      norm=exp(norm);
  
//      cout<<"norma: "<<norm<<endl; 

     }
 

void Wave_Functions::Laugh_Test_B(dcmplx &norm, int ne, int LL, dcmplx *zu, dcmplx *zv, int nu, int npart)
     {
      int i,j;
      double QLaugh;
      dcmplx expo, shift;

      QLaugh=(ne-1.0)/2.0;

      norm=0.0;

      for (i=0;i<npart;++i)
          for (j=i+1;j<npart;++j)
              norm+= log(zu[i]/zv[i]-zu[j]/zv[j]);

      expo=0.0;

      for (i=0;i<npart;++i)
          expo+= 2.0*QLaugh*log(zv[i]);


     shift=0.0;
     for (i=0;i<npart;++i)
         shift+=double(ne-npart)*log(zu[i]/zv[i]);
       
     norm=double(nu)*(norm + expo + shift);
   
  //   norm=exp(norm);

  
     }







void Wave_Functions::Laugh_Sphere(dcmplx &norm, int ne, int LL, dcmplx *zu, dcmplx *zv, int nu)
{

      int i,j;
      double QLaugh;
      dcmplx expo;

      QLaugh=(ne-1.0)/2.0;

   
      norm=0.0;
      for (i=0;i<ne-1;++i)
          for (j=i+1;j<ne;++j)
              norm+=log(zu[i]/zv[i]-zu[j]/zv[j]);


      expo=0.0;
      for (i=0;i<ne;++i)
          expo+=2.0*QLaugh*log(zv[i]);

            
      norm=double(nu)*(norm+expo);

}
















