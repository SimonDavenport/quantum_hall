#include <iostream>
#include <cstdlib>
#include <math.h>
#include <memory.h>
#include <complex>
#include <vector>
#include "hilbert_ewf.hpp"

using namespace std;


void EWF_Space(double Lz, int ne, int &dim, int npart, int LL, vector <vector <double> > &Fspace, 
               vector <vector <double> > &Fspace_comp, int nu)
     {
      double nee=ne, Lzz,L1,L2, QJain, QLaugh, Lz_Jain, Lz_Laugh;
      int i,j,k,l,dim1stLL, dim2ndLL, N_1, N_2, N_3, Lz_1, Lz_2, *sgn_bin1, *sgn_bin2, dim_aux, aux, dim1, dim2, N1_Jain, N2_Jain, 
          dim_eff_Fspace, Lz_3, dim3, m, dim3rdLL, N3_Jain;
      vector <vector <double> > Fspace1, Fspace2, Fspace1_inv, Fspace2_inv, Fspace_comp1, Fspace_comp2, Fspace3, Fspace_comp3, Fspace3_inv;
      vector< vector <int> > bin1, bin2, bin3, bin_comp1, bin_comp2, bin_comp3;
      vector <double> extraLz;


if (LL==1 && nu==1)
   { 
    if (npart<=16)
        partition(Lz+ ((nee-1.0)/2.0)*npart +npart, npart, ne, Fspace, dim1stLL);
    else 
        {
         int npart_eff=15 ;
         double Lz_extra=(npart-npart_eff)*(npart-npart_eff-1)/2.0 + (npart-npart_eff) + (npart-npart_eff)*npart_eff;

         for (i=0;i<npart-npart_eff;++i)
              extraLz.push_back(i+1);

         partition(Lz+ ((nee-1.0)/2.0)*npart +npart - Lz_extra, npart_eff, ne, Fspace, dim1stLL);  

     
         for (i=0;i<dim1stLL;++i)
             {
              for (j=0;j<npart_eff;++j)
                   Fspace[i][j]=Fspace[i][j]+npart-npart_eff;

              Fspace[i].insert(Fspace[i].begin(),extraLz.begin(),extraLz.end());
              }
         }




    for (i=0;i<dim1stLL;++i)
        {
         Fspace[i].resize(npart+LL+1);
         Fspace[i].insert(Fspace[i].begin()+1,Fspace[i].begin(), Fspace[i].end());
         Fspace[i][0]=npart;
         for (j=0+LL;j<npart+LL;++j)
             Fspace[i][j]=Fspace[i][j]-(nee-1)/2.0-1.0;

        }
   
 
    dim=dim1stLL;
 

   
    //////////////// GENERATE SIGN OF PERMUTATION ///////////////////


    if (dim!=-1)
       {
        binary(Fspace, bin1, ne, dim, LL);
        sgn_bin1=new int [dim+1];
   
        sign_binary(bin1, sgn_bin1, ne, dim, LL);

        for (i=0;i<dim;++i)
            {
             Fspace[i][Fspace[i][0]+LL]=1.0;
             Fspace[i][Fspace[i][0]+LL]*=sgn_bin1[i];
             }

        delete [] sgn_bin1;
 
   //////////////// GENERATE COMPLEMENT //////////////////////////////////////
   
        bin_complement(bin1, bin_comp1, ne, dim, LL);     
        inv_binary(Fspace_comp, bin_comp1, ne, dim, LL);

        for (i=0; i<dim; ++i)
            Fspace_comp[i][ne-npart+LL]=1;

        }
   }





else if (LL==1 && nu==2)
        {


        dim=0;      
        QLaugh=(nee-1.0)/2.0;
        Lz_Laugh=(-npart*(npart-1.0)/2.0 + QLaugh*npart);
        dim_eff_Fspace=0;



        for (Lz_1=0;Lz_1<2*Lz_Laugh+1;++Lz_1)
            for (Lz_2=0;Lz_2<2*Lz_Laugh+1;++Lz_2)   
                {


                 if ( (Lz_1+Lz_2-2*Lz_Laugh)==Lz )
                    {

                     EWF_Space(Lz_1-Lz_Laugh, ne, dim1, npart, 1, Fspace1, Fspace_comp1, 1);                    
                     EWF_Space(Lz_2-Lz_Laugh, ne, dim2, npart, 1, Fspace2, Fspace_comp2, 1);

                                  
                     if (dim1==-1 || dim2==-1)
                         continue;
                     else{
                          dim_eff_Fspace+=(dim1+1)*(dim2+1);                   

                          Fspace.resize(dim_eff_Fspace);
                          for (i=0;i<dim_eff_Fspace;++i)
                               Fspace[i].resize(nu*npart+LL+1);

                          Fspace_comp.resize(dim_eff_Fspace);
                          for (i=0;i<dim_eff_Fspace;++i)
                              Fspace_comp[i].resize(nu*(ne-npart)+LL+1);


            
                          for (i=0;i<dim1;++i)
                              for (j=0;j<dim2;++j)
                                  {

                                   Fspace[dim][0]=Fspace1[i][0];
                                   for (m=0;m<npart;++m)
                                       {                                   
                                       Fspace[dim][m+LL]=Fspace1[i][m+LL];
                                       Fspace[dim][m+LL+npart]=Fspace2[j][m+LL];
                                        }                                


                                   Fspace[dim][LL+nu*npart]=Fspace1[i][npart+LL]*Fspace2[j][npart+LL];

                                   dim=dim+1;
                    
                                   }
                           }
                       }

             }



          Fspace1.resize(dim+1);
          Fspace2.resize(dim+1);
          
          for (i=0;i<dim;++i)
              {
               Fspace1[i].resize(ne+LL);
               Fspace2[i].resize(ne+LL);          
               }

         for (i=0;i<dim;++i)
             for (j=0;j<ne+LL;++j)
                  {
                   Fspace1[i][j]=0.0;
                   Fspace2[i][j]=0.0;
                   }



        
          for (i=0;i<dim;++i)
              {
               Fspace1[i][0]=Fspace[i][0];
               Fspace2[i][0]=Fspace[i][0];
               for (j=0;j<npart;++j)
                   {
                    Fspace1[i][j+LL]=Fspace[i][j+LL];
                    Fspace2[i][j+LL]=Fspace[i][j+LL+npart];
                    }            
              }



          ///////// ELLIMINATE EQUALS CONFIGURATIONS AND COMPUTE THE CORRESPONDING MULTIPLICITIES ////////////


	  //          for (i=0;i<dim;++i)
            //            Fspace[i][3*Fspace[i][0]+LL]=1;
       
                    for (i=0;i<dim-1;++i)      
                        {
                         aux=1;
                         for (j=i+1;j<dim;++j)
                             {   
                              if ((Fspace1[i]==Fspace2[j]) && (Fspace2[i]==Fspace1[j]))
                                 {      
                                  aux=aux+1;                                 
                                  Fspace.erase(Fspace.begin()+j);
                                  Fspace1.erase(Fspace1.begin()+j);
                                  Fspace2.erase(Fspace2.begin()+j);
                                  dim=dim-1;
                                  --j;
                                 }                             
                              }
                         Fspace[i][nu*Fspace[i][0]+LL]*=aux;
                         }




        //////////// COMPUTE COMPLEMENT STATES, i.e. STATES ON B SUBSYSTEM /////////////


         binary(Fspace1, bin1, ne, dim, LL);
         binary(Fspace2, bin2, ne, dim, LL);


         bin_complement(bin1, bin_comp1, ne, dim, LL);
         bin_complement(bin2, bin_comp2, ne, dim, LL);
  

         inv_binary(Fspace1_inv, bin_comp1, ne, dim, LL);
         inv_binary(Fspace2_inv, bin_comp2, ne, dim, LL);






        Fspace_comp.resize(dim+1);
        for (i=0;i<dim;++i)
            Fspace_comp[i].resize(nu*ne+LL+1);


     
        for (i=0;i<dim;++i)
            {
             Fspace_comp[i][0]=nee-Fspace[i][0];
             for (j=0;j<Fspace_comp[i][0];++j)
                 {
                  Fspace_comp[i][j+LL]=Fspace1_inv[i][j+LL];
                  Fspace_comp[i][j+LL+Fspace_comp[i][0]]=Fspace2_inv[i][j+LL];
                  }
             Fspace_comp[i][2*(ne-npart)+LL]=1;
             }






         }






else if (LL==1 && nu==3)
        {


        dim=0;      
        QLaugh=(nee-1.0)/2.0;
        Lz_Laugh=(-npart*(npart-1.0)/2.0 + QLaugh*npart);
        dim_eff_Fspace=0;



        for (Lz_1=0;Lz_1<2*Lz_Laugh+1;++Lz_1)
            for (Lz_2=0;Lz_2<2*Lz_Laugh+1;++Lz_2)   
                 for (Lz_3=0;Lz_3<2*Lz_Laugh+1;++Lz_3)   
                    {
                    if ( (Lz_1+Lz_2+Lz_3-3*Lz_Laugh)==Lz )
                      {

                    
                     EWF_Space(Lz_1-Lz_Laugh, ne, dim1, npart, 1, Fspace1, Fspace_comp1, 1);                    
                     EWF_Space(Lz_2-Lz_Laugh, ne, dim2, npart, 1, Fspace2, Fspace_comp2, 1);
                     EWF_Space(Lz_3-Lz_Laugh, ne, dim3, npart, 1, Fspace3, Fspace_comp3, 1);

                     if (dim1==-1 || dim2==-1 || dim3==-1)
                         continue;
                     else
                         {
                          dim_eff_Fspace+=(dim1+1)*(dim2+1)*(dim3+1);                   

                          Fspace.resize(dim_eff_Fspace);
                          for (i=0;i<dim_eff_Fspace;++i)
                               Fspace[i].resize(3*npart+LL+1);

                          Fspace_comp.resize(dim_eff_Fspace);
                          for (i=0;i<dim_eff_Fspace;++i)
                               Fspace_comp[i].resize(3*(ne-npart)+LL+1);

            
                          for (i=0;i<dim1;++i)
                              for (j=0;j<dim2;++j)
                                  for (k=0;k<dim3;++k)
                                      {
                                       Fspace[dim][0]=Fspace1[i][0];
                                       for (m=0;m<npart;++m)
                                           {                                   
                                            Fspace[dim][m+LL]=Fspace1[i][m+LL];
                                            Fspace[dim][m+LL+npart]=Fspace2[j][m+LL];
                                            Fspace[dim][m+LL+2*npart]=Fspace3[k][m+LL];
                                            }                                
                                        Fspace[dim][LL+3*npart]=Fspace1[i][npart+LL]*Fspace2[j][npart+1]*Fspace3[k][npart+1];     
                                        dim=dim+1;                              
                                       }
                            }
                       }

             }





          Fspace1.resize(dim+1);
          Fspace2.resize(dim+1);
          Fspace3.resize(dim+1);
          for (i=0;i<dim;++i)
              {
               Fspace1[i].resize(ne+LL);
               Fspace2[i].resize(ne+LL);
               Fspace3[i].resize(ne+LL);
               }


         for (i=0;i<dim;++i)
             for (j=0;j<ne+LL;++j)
                  {
                   Fspace1[i][j]=0.0;
                   Fspace2[i][j]=0.0;
                   Fspace3[i][j]=0.0;
                   }


        
          for (i=0;i<dim;++i)
              {
               Fspace1[i][0]=Fspace[i][0];
               Fspace2[i][0]=Fspace[i][0];
               Fspace3[i][0]=Fspace[i][0];
               for (j=0;j<npart;++j)
                   {
                    Fspace1[i][j+LL]=Fspace[i][j+LL];
                    Fspace2[i][j+LL]=Fspace[i][j+LL+npart];
                    Fspace3[i][j+LL]=Fspace[i][j+LL+2*npart];
                    }            
              }




          ///////// ELLIMINATE EQUALS CONFIGURATIONS AND COMPUTE THE CORRESPONDING MULTIPLICITIES ////////////

       
                    for (i=0;i<dim-1;++i)      
                        {
                         aux=1;
                         for (j=i+1;j<dim;++j)
                             {   
                              if ((Fspace1[i]==Fspace2[j] && Fspace2[i]==Fspace1[j] && Fspace3[i]==Fspace3[j]) ||
                                  (Fspace1[i]==Fspace2[j] && Fspace2[i]==Fspace3[j] && Fspace3[i]==Fspace1[j]) ||
                                  (Fspace1[i]==Fspace3[j] && Fspace2[i]==Fspace1[j] && Fspace3[i]==Fspace2[j]) ||
                                  (Fspace1[i]==Fspace3[j] && Fspace2[i]==Fspace2[j] && Fspace3[i]==Fspace1[j]) ||
                                  (Fspace1[i]==Fspace1[j] && Fspace2[i]==Fspace3[j] && Fspace3[i]==Fspace2[j]))
                                 {
                                  aux=aux+1;                                 
                                  Fspace.erase(Fspace.begin()+j);
                                  Fspace1.erase(Fspace1.begin()+j);
                                  Fspace2.erase(Fspace2.begin()+j);
                                  Fspace3.erase(Fspace3.begin()+j);
                                  dim=dim-1;
                                  --j;
                                 }                             
                              }
                         Fspace[i][3*Fspace[i][0]+LL]*=aux;
                         }




        //////////// COMPUTE COMPLEMENT STATES, i.e. STATES ON B SUBSYSTEM /////////////


         binary(Fspace1, bin1, ne, dim, LL);
         binary(Fspace2, bin2, ne, dim, LL);
         binary(Fspace3, bin3, ne, dim, LL);


         bin_complement(bin1, bin_comp1, ne, dim, LL);
         bin_complement(bin2, bin_comp2, ne, dim, LL);
         bin_complement(bin3, bin_comp3, ne, dim, LL);



         inv_binary(Fspace1_inv, bin_comp1, ne, dim, LL);
         inv_binary(Fspace2_inv, bin_comp2, ne, dim, LL);
         inv_binary(Fspace3_inv, bin_comp3, ne, dim, LL);




        Fspace_comp.resize(dim+1);
        for (i=0;i<dim;++i)
            Fspace_comp[i].resize(nu*ne+LL);


     
        for (i=0;i<dim;++i)
            {
             Fspace_comp[i][0]=nee-Fspace[i][0];
             for (j=0;j<Fspace_comp[i][0];++j)
                 {
                  Fspace_comp[i][j+LL]=Fspace1_inv[i][j+LL];
                  Fspace_comp[i][j+LL+Fspace_comp[i][0]]=Fspace2_inv[i][j+LL];
                  Fspace_comp[i][j+LL+2*Fspace_comp[i][0]]=Fspace3_inv[i][j+LL];
                  }
             Fspace_comp[i][3*(ne-npart)+LL]=1;
             }

             


         }







else if (LL==2 && nu==0)
        {
    
         dim=0;
         nee=ne;
         double QJain_0th=(nee/2.0-2.0)/2.0;
         double QJain_1st=nee/4.0;

         for (N_1=0; N_1<=npart ; ++N_1)
             for (N_2=0; N_2<=npart ; ++N_2)
                 {  
                  if ((N_1+N_2)==npart && N_1<=nee/2.0-1 && N_2<=nee/2.0+1)
                     {                 

                      double Lz_N1= -N_1*(N_1-1.0)/2.0 + N_1*QJain_0th;
                      double Lz_N2= -N_2*(N_2-1.0)/2.0 + N_2*QJain_1st;

                     for (Lz_1=0; Lz_1<=2.0*Lz_N1; ++Lz_1)   ///angular mom. in LLL
                         for (Lz_2=0; Lz_2<=2.0*Lz_N2; ++Lz_2)         ///angular mom. in 1st LL
                             if ((Lz_1+Lz_2-Lz_N1-Lz_N2)==Lz ) 
                                {             

                                   EWF_Space(Lz_1-Lz_N1, nee/2.0-1.0, dim1stLL, N_1, 1, Fspace1, Fspace_comp1, 1);
                                   EWF_Space(Lz_2-Lz_N2, nee/2.0+1.0, dim2ndLL, N_2, 1, Fspace2, Fspace_comp2, 1);
                                   
            
                                   if (dim1stLL==-1 && N_1!=0 || dim2ndLL==-1 && N_2!=0)
                                      continue;
                                   else
                                       {
                                        if (dim1stLL==-1)
                                            dim1stLL=1;
                                        if (dim2ndLL==-1)
                                            dim2ndLL=1;

                                        for (i=0; i<dim1stLL; ++i)
                                             for (j=0; j<dim2ndLL; ++j)
                                                 {                                   
                                                  Fspace.resize(dim+1);
                                                  Fspace[dim].resize(npart+LL+1);
                                                  Fspace[dim][0]=N_1;
                                                  Fspace[dim][1]=N_2;
                                                  for (k=0;k<N_1;++k)     
                                                       Fspace[dim][k+LL]=Fspace1[i][k+1];
                                                  for (k=0;k<N_2;++k)     
                                                       Fspace[dim][k+N_1+LL]=Fspace2[j][k+1];
                                                  Fspace[dim][npart+LL]=1.0;
                                                  dim=dim+1;  
                                                  }
                                        } 







                                }
                     }
               }     

  



       //////////////// GENERATE SIGN OF PERMUTATION ///////////////////



        binary(Fspace, bin1, ne, dim, LL);

        sgn_bin1=new int [dim+1];
        sign_binary(bin1, sgn_bin1, ne, dim, LL);


        for (i=0;i<dim;++i)
            {
             Fspace[i][Fspace[i][0]+Fspace[i][1]+LL]=sgn_bin1[i];
             }

        delete [] sgn_bin1;

      //////////////// GENERATE COMPLEMENT //////////////////////////////////////
    

        bin_complement(bin1, bin_comp1, ne, dim, LL);

    
        inv_binary(Fspace_comp, bin_comp1, ne, dim, LL);    

        for (i=0; i<dim; ++i)
            Fspace_comp[i][ne-npart+LL]=1.0;    
       
        }






else if (LL==2 && (nu==2||nu==1))
        {
        dim=0;      

        QJain=(nee/2.0-2.0)/2.0;
        QLaugh=(nee-1.0)/2.0;

        N1_Jain=int((2.0*QJain+1.0)/2.0);
        N2_Jain=npart-N1_Jain;
        Lz_Jain= -N1_Jain*(N1_Jain-1.0)/2.0 + QJain*N1_Jain - N2_Jain*(N2_Jain-1.0)/2.0 + (QJain+1.0)*N2_Jain;
        Lz_Laugh=( -npart*(npart-1.0)/2.0 + QLaugh*npart)*nu;
        dim_eff_Fspace=0;



      if (npart==nee)
        {
         Lz_Jain=0.0;
         Lz_Laugh=0.0;
         }

         
        for (Lz_1=0;Lz_1<2*Lz_Jain+1;++Lz_1)
            for (Lz_2=0;Lz_2<2*Lz_Laugh+1;++Lz_2)   
                {


                 if ( (Lz_1+Lz_2-Lz_Jain-Lz_Laugh)==Lz )
                    {        


                     EWF_Space(Lz_1-Lz_Jain, ne, dim1, npart, 2, Fspace1, Fspace_comp1, 0);                                         
                     EWF_Space(Lz_2-Lz_Laugh, ne, dim2, npart, 1, Fspace2, Fspace_comp2, nu);

                     
                    dim_eff_Fspace+=(dim1+1)*(dim2+1);                   

                     Fspace.resize(dim_eff_Fspace);
                     for (i=0;i<dim_eff_Fspace;++i)
                         Fspace[i].resize((nu+1)*npart+LL+1);

                     Fspace_comp.resize(dim_eff_Fspace);
                     for (i=0;i<dim_eff_Fspace;++i)
                         Fspace_comp[i].resize((nu+1)*(ne-npart)+LL+1);

            
                     for (i=0;i<dim1;++i)
                         for (j=0;j<dim2;++j)
                             {

                              Fspace[dim][0]=Fspace1[i][0];
                              Fspace[dim][1]=Fspace1[i][1];
                              for (k=0;k<npart;++k)                                   
                                   Fspace[dim][k+LL]=Fspace1[i][k+LL];                                                                   
                              for (k=0;k<nu*npart;++k)                                   
                                   Fspace[dim][k+LL+npart]=Fspace2[j][k+1];

                              Fspace[dim][LL+(nu+1)*npart]=Fspace1[i][npart+LL]*Fspace2[j][nu*npart+1];     
                    
                   
                              Fspace_comp[dim][0]=Fspace_comp1[i][0];
                              Fspace_comp[dim][1]=Fspace_comp1[i][1];
                              for (k=0;k<ne-npart;++k)                                   
                                   Fspace_comp[dim][k+LL]=Fspace_comp1[i][k+LL];
                              for (k=0;k<nu*(ne-npart);++k)                                   
                                   Fspace_comp[dim][k+LL+ne-npart]=Fspace_comp2[j][k+1];

                              Fspace_comp[dim][LL+(nu+1)*(ne-npart)]=1;     

                              dim=dim+1;
                              }

                     }
                  
              

                 }




         }

else if (LL==3 && nu==0)
        {
    
         dim=0;
         nee=ne;
         double LLd=LL;
         double QJain_0th=(nee/LLd-LLd)/2.0;
         double QJain_1st=(nee/LLd-LLd+2)/2.0;
         double QJain_2nd=(nee/LLd-LLd+4)/2.0;

         for (N_1=0; N_1<=npart ; ++N_1)
             for (N_2=0; N_2<=npart ; ++N_2)
                 for (N_3=0; N_3<=npart ; ++N_3)
                 {  
                  if ((N_1+N_2+N_3)==npart && N_1<=nee/LLd-2 && N_2<=nee/LLd && N_3<=nee/LLd+2)
                     {                 

                      double Lz_N1= -N_1*(N_1-1.0)/2.0 + N_1*QJain_0th;
                      double Lz_N2= -N_2*(N_2-1.0)/2.0 + N_2*QJain_1st;
                      double Lz_N3= -N_3*(N_3-1.0)/2.0 + N_3*QJain_2nd;

                     for (Lz_1=0; Lz_1<=2.0*Lz_N1; ++Lz_1)   ///angular mom. in LLL
                         for (Lz_2=0; Lz_2<=2.0*Lz_N2; ++Lz_2)         ///angular mom. in 1st LL
                             for (Lz_3=0; Lz_3<=2.0*Lz_N3; ++Lz_3)         ///angular mom. in 2nd LL
                                 if ((Lz_1+Lz_2+Lz_3-Lz_N1-Lz_N2-Lz_N3)==Lz) 
                                    {             
                                     EWF_Space(Lz_1-Lz_N1, nee/LLd-2.0, dim1stLL, N_1, 1, Fspace1, Fspace_comp1, 1);
                                     EWF_Space(Lz_2-Lz_N2, nee/LLd, dim2ndLL, N_2, 1, Fspace2, Fspace_comp2, 1);
                                     EWF_Space(Lz_3-Lz_N3, nee/LLd+2.0, dim3rdLL, N_3, 1, Fspace3, Fspace_comp3, 1);
                                   

            
                                   if (dim1stLL==-1 && N_1!=0 || dim2ndLL==-1 && N_2!=0 || dim3rdLL==-1 && N_3!=0)
                                      continue;
                                   else
                                       {                                     

                                        if (dim1stLL==-1)
                                            dim1stLL=1;
                                        if (dim2ndLL==-1)
                                            dim2ndLL=1;
                                        if (dim3rdLL==-1)
                                            dim3rdLL=1;

                                       
                                        for (i=0; i<dim1stLL; ++i)
                                             for (j=0; j<dim2ndLL; ++j)
                                                 for (l=0; l<dim3rdLL; ++l)
                                                     {                                                                        
                                                      Fspace.resize(dim+1);
                                                      Fspace[dim].resize(npart+LL+1);
                                                      Fspace[dim][0]=N_1;
                                                      Fspace[dim][1]=N_2;
                                                      Fspace[dim][2]=N_3;
                                                      for (k=0;k<N_1;++k)     
                                                           Fspace[dim][k+LL]=Fspace1[i][k+1];
                                                      for (k=0;k<N_2;++k)     
                                                           Fspace[dim][k+N_1+LL]=Fspace2[j][k+1];
                                                      for (k=0;k<N_3;++k)     
                                                           Fspace[dim][k+N_1+N_2+LL]=Fspace3[l][k+1];

                                                  Fspace[dim][npart+LL]=1.0;
                                                  dim=dim+1;  
                                                  }
                                        } 



                                }
                     }
               }     






       //////////////// GENERATE SIGN OF PERMUTATION ///////////////////



        binary(Fspace, bin1, ne, dim, LL);

        sgn_bin1=new int [dim+1];
        sign_binary(bin1, sgn_bin1, ne, dim, LL);


        for (i=0;i<dim;++i)
            {
             Fspace[i][Fspace[i][0]+Fspace[i][1]+Fspace[i][2]+LL]=sgn_bin1[i];
             }

        delete [] sgn_bin1;

      //////////////// GENERATE COMPLEMENT //////////////////////////////////////
    

        bin_complement(bin1, bin_comp1, ne, dim, LL);
    
        inv_binary(Fspace_comp, bin_comp1, ne, dim, LL);    

        for (i=0; i<dim; ++i)
            Fspace_comp[i][ne-npart+LL]=1.0;    

               
        }






  
else if (LL==3 && (nu==2||nu==1))
        {
        dim=0;      

        QJain=(nee/3.0-3.0)/2.0;
        QLaugh=(nee-1.0)/2.0;




        N1_Jain=int((2.0*QJain+1.0)/2.0);
        N2_Jain=int((2.0*QJain+3.0)/2.0);
        N3_Jain=npart-N1_Jain-N2_Jain;

        Lz_Jain= -N1_Jain*(N1_Jain-1.0)/2.0 + QJain*N1_Jain - N2_Jain*(N2_Jain-1.0)/2.0 + (QJain+1.0)*N2_Jain 
                - N3_Jain*(N3_Jain-1.0)/2.0 + (QJain+2.0)*N3_Jain;;

        Lz_Laugh=( -npart*(npart-1.0)/2.0 + QLaugh*npart)*nu;
        dim_eff_Fspace=0;



      if (npart==nee)
        {
         Lz_Jain=0.0;
         Lz_Laugh=0.0;
         }


         
        for (Lz_1=0;Lz_1<2*Lz_Jain+1;++Lz_1)
            for (Lz_2=0;Lz_2<2*Lz_Laugh+1;++Lz_2)   
                {
                 if ( (Lz_1+Lz_2-Lz_Jain-Lz_Laugh)==Lz )
                    {        

                     EWF_Space(Lz_1-Lz_Jain, ne, dim1, npart, 3, Fspace1, Fspace_comp1, 0);                                         
                     EWF_Space(Lz_2-Lz_Laugh, ne, dim2, npart, 1, Fspace2, Fspace_comp2, nu);

                     
                    dim_eff_Fspace+=(dim1+1)*(dim2+1);                   

                     Fspace.resize(dim_eff_Fspace);
                     for (i=0;i<dim_eff_Fspace;++i)
                         Fspace[i].resize((nu+1)*npart+LL+1);

                     Fspace_comp.resize(dim_eff_Fspace);
                     for (i=0;i<dim_eff_Fspace;++i)
                         Fspace_comp[i].resize((nu+1)*(ne-npart)+LL+1);

            
                     for (i=0;i<dim1;++i)
                         for (j=0;j<dim2;++j)
                             {

                              Fspace[dim][0]=Fspace1[i][0];
                              Fspace[dim][1]=Fspace1[i][1];
                              Fspace[dim][2]=Fspace1[i][2];
                              for (k=0;k<npart;++k)                                   
                                   Fspace[dim][k+LL]=Fspace1[i][k+LL];                                                                   
                              for (k=0;k<nu*npart;++k)                                   
                                   Fspace[dim][k+LL+npart]=Fspace2[j][k+1];

                              Fspace[dim][LL+(nu+1)*npart]=Fspace1[i][npart+LL]*Fspace2[j][nu*npart+1];     
                    
                   
                              Fspace_comp[dim][0]=Fspace_comp1[i][0];
                              Fspace_comp[dim][1]=Fspace_comp1[i][1];
                              Fspace_comp[dim][2]=Fspace_comp1[i][2];
                              for (k=0;k<ne-npart;++k)                                   
                                   Fspace_comp[dim][k+LL]=Fspace_comp1[i][k+LL];
                              for (k=0;k<nu*(ne-npart);++k)                                   
                                   Fspace_comp[dim][k+LL+ne-npart]=Fspace_comp2[j][k+1];

                              Fspace_comp[dim][LL+(nu+1)*(ne-npart)]=1;     

                              dim=dim+1;
                              }

                     }
                  
              

                 }




         }


     }











void bin_complement(vector< vector <int> > bin, vector< vector <int> > & bin_comp, int ne, int dim, int LL)
{
int i,j,k;


bin_comp.resize(dim);
for(i=0;i<dim;i++)
bin_comp[i].resize(ne+LL);


for (i=0;i<dim;++i)
    {
     if (LL==1)
        bin_comp[i][0]=ne-bin[i][0];
     if (LL==2)
        {
         bin_comp[i][0]=abs(ne/2-1 - bin[i][0]);
         bin_comp[i][1]=abs(ne/2+1 - bin[i][1]);
         }
     if (LL==3)
        {
         bin_comp[i][0]=abs(ne/3-2 - bin[i][0]);
         bin_comp[i][1]=abs(ne/3 - bin[i][1]);
         bin_comp[i][2]=abs(ne/3+2 - bin[i][2]);
         }


    for (j=0+LL;j<ne+LL;++j)
        {
         if (bin[i][j]==1)
            { 
             bin_comp[i][j]=0;
             }
          else if (bin[i][j]==0)
                  {
                   bin_comp[i][j]=1;
                   }
         }
      }

}



void sign_binary(vector< vector <int> > bin, int *sgn_bin, int ne, int dim, int LL)
{
int i,j,k,signo,aux;
double npart;



for(i=0;i<dim;i++)
   {

    if (LL==1)
       {npart=bin[i][0];}
    else if (LL==2)
            {npart=bin[i][0]+bin[i][1];}
    else if (LL==3)
            {npart=bin[i][0]+bin[i][1]+bin[i][2];}


    signo=1;
    aux=-1;

    for(j=0+LL;j<ne+LL;j++)
       {        
        if (bin[i][j]==1 && aux<npart-1)           
            {
             aux=aux+1;
             signo=signo*pow(-1.0,j-LL-aux);
             }
       
         }
     
    sgn_bin[i]=signo;
    }

}


void inv_binary(vector< vector <double> > &Fspace, vector< vector <int> > bin, int ne, int dim, int LL)
{     
int i,j,k,count,angular;
double nee=ne;


Fspace.resize(dim);
for(i=0;i<dim;i++)
Fspace[i].resize(ne+LL+1);


if (LL==1)
   {

    for (i=0;i<dim;i++)
        {
         count=0;
         angular=0;
         Fspace[i][0]=bin[i][0];

         for(j=0+LL;j<ne+LL;++j)
            {
             if (bin[i][j]==1)         
                {
                 count=count+1;             
                 Fspace[i][count+LL-1]=angular-(nee-1)/2.0;
                 }
             angular=angular+1;
             }
          }        

    }
else if (LL==2)
        {

         for (i=0;i<dim;i++)
             {
              count=0;
              angular=0;
              Fspace[i][0]=bin[i][0];
              Fspace[i][1]=bin[i][1];

              for(j=0+LL;j<ne+LL;++j)
                 {
                  if (bin[i][j]==1)         
                     {
                      count=count+1;             
                      if (count<=Fspace[i][0])
                         {
                          Fspace[i][count+LL-1]=angular-(nee/4.0-1.0);
                          }
                      else
                         {
                          Fspace[i][count+LL-1]=angular-(nee/2.0-1.0)-nee/4.0;   
                          }
                      }
                  angular=angular+1;
                 }
              }
         }


////////////////// SEGUIR DESDE ACA ///////////////////////////




else if (LL==3)
        {

         double LLd=LL;
         double QJain_0th=(nee/LLd-LLd)/2.0;
         double QJain_1st=(nee/LLd-LLd+2)/2.0;
         double QJain_2nd=(nee/LLd-LLd+4)/2.0;



         for (i=0;i<dim;++i)
             {

              count=0;
              angular=0;
              Fspace[i][0]=bin[i][0];
              Fspace[i][1]=bin[i][1];
              Fspace[i][2]=bin[i][2];

              for(j=0+LL;j<ne+LL;j++)
                 {
                  if (bin[i][j]==1)         
                     {
                      count=count+1;             
                      if (count<=Fspace[i][0])
                         {
                          Fspace[i][count+LL-1]=angular-QJain_0th;
                    
                          }
                      else if (count>Fspace[i][0] && count <=(Fspace[i][0]+Fspace[i][1]))
                         {
                          Fspace[i][count+LL-1]=angular-2.0*QJain_0th-1-QJain_1st;   
                    
                          }
                      else if (count>Fspace[i][0]+Fspace[i][1])
                         {
                         Fspace[i][count+LL-1]=angular-2.0*QJain_0th-1-2.0*QJain_1st-1-QJain_2nd;   
                    
                          }

                      }
                  angular=angular+1;



                 }
              }


         }






}



void binary(vector< vector <double> > Fspace, vector< vector <int> > &bin, int ne, int dim, int LL)
{

int i,j,k;
double nee=ne;


   bin.resize(dim+1);
   for (i=0;i<dim+1;++i)
        bin[i].resize(4*ne);

   for (i=0;i<dim+1;++i)
       for (j=0;j<4*ne;++j)
           bin[i][j]=0;


if (LL==1)
   { 
        
    for(i=0;i<dim;++i)
       {
       bin[i][0]=Fspace[i][0];

       for(j=0;j<Fspace[i][0];++j)
          {
           bin[i][(int)(Fspace[i][j+LL]+LL+(nee-1.0)/2.0)]=1;
           }
       }
    }

else if (LL==2)
        {
    
         for (i=0; i<dim; ++i)
             { 
              bin[i][0]=Fspace[i][0];
              bin[i][1]=Fspace[i][1];


              for (j=0; j<Fspace[i][0]; ++j)
                   bin[(int) i][(int)( Fspace[i][j+LL] + (nee/4.0 -1)+LL)]=1;
         
              for (j=0; j<Fspace[i][1]; ++j)
                   bin[(int) i][(int)( Fspace[i][j+LL+Fspace[i][0]] + LL + (nee/4.0+1)+(nee/2.0-2))]=1;          
              }

          }




else if (LL==3)
        {
    
         double LLd=LL;
         double QJain_0th=(nee/LLd-LLd)/2.0;
         double QJain_1st=(nee/LLd-LLd+2)/2.0;
         double QJain_2nd=(nee/LLd-LLd+4)/2.0;


         for (i=0; i<dim; ++i)
             { 
              bin[i][0]=Fspace[i][0];
              bin[i][1]=Fspace[i][1];
              bin[i][2]=Fspace[i][2];


              for (j=0; j<Fspace[i][0]; ++j)
                   bin[(int) i][(int)( Fspace[i][j+LL] + QJain_0th +LL)]=1;
         
              for (j=0; j<Fspace[i][1]; ++j)
                   bin[(int) i][(int)( Fspace[i][j+LL+Fspace[i][0]] + LL + 2.0*QJain_0th + 1 + QJain_1st)]=1;          

              for (j=0; j<Fspace[i][2]; ++j)
                   bin[(int) i][(int)( Fspace[i][j+LL+Fspace[i][0]+Fspace[i][1]] + LL + 2.0*QJain_0th + 1 + 2.0* QJain_1st + 1 +
                       QJain_2nd)]=1;          


              }

          }



}     






void partition(int part, int npart, int norbit, vector< vector <double> > &Fspace, int &dim_space)
{
int *parts;
int *ptr;
int i;
int idx = 0;
int tot = 0;
int cur = 1;
int max = 1;


dim_space=-1;

    

    while((max*(max + 1)) / 2 <= part) max++;

    Fspace.resize(1);
    Fspace[0].resize(max);

  
    
    ptr = parts = (int*) malloc(max*sizeof(int));
   

    for(;;) {
        if((tot += *(ptr++) = cur++) < part) continue;


        if(tot == part && ptr-parts==npart && parts[npart-1]<=norbit) 
           {
    
            if (dim_space==-1)
                dim_space=0;    

            for(i = 0 ; i < ptr-parts ; ++i) {           
                                              Fspace[dim_space][i]=parts[i];                                              
                                              }                                                                
            dim_space+=1;    

            Fspace.resize(dim_space+1);          
            Fspace[dim_space].resize(max);
    
            }


        do {
            if(ptr == parts) {free(parts); return;}
         tot -= cur = *--ptr;

        } while(++cur + tot > part);

    }

}






