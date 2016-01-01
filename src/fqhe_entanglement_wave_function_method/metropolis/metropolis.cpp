
#include "metropolis.hpp"

using namespace std;
typedef std::complex<double> dcmplx;

const long double PI=3.1415926535897932384626433832795028841972;




metropolis::metropolis(dcmplx *zz1, int np1, int iter1, double width1, int np_A1, int nu1, int LL1, int nu_eff1, double cut_sph1)
{


int i, dim;
vector< vector <double> > Fspace, Fspace_comp;


zz=zz1;
np=np1;
iter=iter1;
width=width1;
np_A=np_A1;
nu=nu1;
LL=LL1;
nu_eff=nu_eff1;
cut_sph=cut_sph1;



zu=new dcmplx [np];
zv=new dcmplx [np];
zu_move=new dcmplx [np];
zv_move=new dcmplx [np];
zz_move=new dcmplx [np];


Matrix1=new dcmplx* [np];
for (i=0;i<np;i++)
    Matrix1[i]=new dcmplx[np]();


wf=new Wave_Functions(np, np, nu, LL, 'P');

EWF_Space(0, np, dim, np, LL, Fspace, Fspace_comp, nu);


Lz_EWF=new double [nu*np+ np + LL+1]();
for (i=0;i<nu*np+np+LL+1;++i)
     Lz_EWF[i]=Fspace[0][i];
     
    

}





metropolis::metropolis(dcmplx *zz1, int np1, int iter1, int width1):np_A(np1/2),nu(1),LL(1),nu_eff(0),cut_sph(PI/2.0)
{

int i, dim;
vector< vector <double> > Fspace, Fspace_comp;

zz=zz1;
np=np1;
iter=iter1;
width=width1;



zu=new dcmplx [np];
zv=new dcmplx [np];
zu_move=new dcmplx [np];
zv_move=new dcmplx [np];
zz_move=new dcmplx [np];


Matrix1=new dcmplx* [np];
for (i=0;i<np;i++)
    Matrix1[i]=new dcmplx[np]();


wf=new Wave_Functions(np, np, nu, LL, 'P');

EWF_Space(0, np, dim, np, LL, Fspace, Fspace_comp, nu);

Lz_EWF=new double [nu*np+ np +LL+1]();
for (i=0;i<nu*np+ np + LL+1;++i)
     Lz_EWF[i]=Fspace[0][i];

}





metropolis::~metropolis()
{

int i;

delete[] zu; delete[]zv; delete[] zu_move; delete[] zv_move, delete[] zz_move;
delete[] Lz_EWF;


for (i=0; i<np; ++i)
    delete [] Matrix1[i];
delete [] Matrix1;

wf-> ~Wave_Functions();


}













void metropolis::run_A()
{ 
	   
	int i,j,i1,t1,cut1, npart_red;
	dcmplx norm1_a,norm1_b;
	double darea_a, darea_b, ratio;    



         accepted=0; rejected=0;

     	
	//	Convert zz to spinor co-ordiantes
	for (i=0;i<np_A;++i)
	{
		zu[i]= cos(real(zz[i]/2.0))*exp(dcmplx(0,imag(zz[i])/2.0));
		zv[i]= sin(real(zz[i]/2.0))*exp(dcmplx(0,-imag(zz[i])/2.0));
	}


//           norm1_a=Jain_negative_2_LL_A (zu,zv);
          norm1_a=IQHE_Prob_A(zu, zv);

 

	 darea_a=1.0;
	 for (i=0; i<np_A; i++)
	 darea_a=darea_a*sin(real(zz[i]));



	for (i1=1;i1<=iter;i1++)
		{
		  for (i=0;i<np_A;i++)
		       zz_move[i]=zz[i];		  

                  t1=mt.genrand_int31()%np_A;       



//		  zz_move[t1] = dcmplx(cut_sph + nfmod(real(zz_move[t1]) + (width)*(mt.random()-0.5),PI-cut_sph) ,
//                                     nfmod(imag(zz_move[t1]) + (width)*(mt.random()-0.5), 2.0*PI));


		  zz_move[t1] = dcmplx(nfmod(real(zz_move[t1]) + (width)*(mt.random()-0.5),PI) ,
                                     nfmod(imag(zz_move[t1]) + (width)*(mt.random()-0.5), 2.0*PI));



	///////////////// INTEGRATION REGION //////////////////////
		
            cut1=0;
            for (i=0;i<np_A;i++)                
		{
                 if (real(zz_move[i])>=cut_sph && real(zz_move[i])<=PI)
	             cut1++;                 
		 }

        ///////////////////////////////////////////////////////////




             if (cut1==np_A)
		 {

		   darea_b=1.0;
		   for (i=0; i<np_A; i++)
		        darea_b=darea_b*sin(real(zz_move[i]));


      			 for (i=0;i<np_A;i++)
			 {
				 zu_move[i]= cos(real(zz_move[i]/2.0))*exp(dcmplx(0,imag(zz_move[i])/2.0));
				 zv_move[i]= sin(real(zz_move[i]/2.0))*exp(dcmplx(0,-imag(zz_move[i])/2.0));
			 }
		



//               norm1_b=Jain_negative_2_LL_A (zu_move, zv_move);

             norm1_b=IQHE_Prob_A(zu_move, zv_move);


             ratio=pow(abs(exp(norm1_b-norm1_a)),2.0)*(darea_b/darea_a);


				   
		   if ((ratio>=1.0 ) || ratio>mt.random())
			  {
			   accepted=accepted+1;
			   zz[t1]=zz_move[t1];
			   norm1_a=norm1_b;              
			   darea_a=darea_b;
			   }
		   else
		            rejected++;
	             }

	          else    
			    rejected++;
	}


	return;   
}








void metropolis::run_B()
    { 
	
        int i,j,i1,t1,cut1, npart_red;
	dcmplx norm1_a,norm1_b;
	double darea_a, darea_b, ratio;    


         accepted=0; rejected=0;  

        
	
	//	Convert zz to spinor co-ordiantes
	for (i=0;i<np-np_A;++i)
	{
		zu[i]= cos(real(zz[i]/2.0))*exp(dcmplx(0,imag(zz[i])/2.0));
		zv[i]= sin(real(zz[i]/2.0))*exp(dcmplx(0,-imag(zz[i])/2.0));
	}


//         norm1_a=Jain_negative_2_LL_B (zu, zv);
            norm1_a=IQHE_Prob_B(zu,zv);



	 darea_a=1.0;
	 for (i=0; i<np-np_A; i++)
	 darea_a=darea_a*sin(real(zz[i]));




	for (i1=1;i1<=iter;i1++)
		{
		  for (i=0;i<np-np_A;i++)
			  zz_move[i]=zz[i];		  

                  t1=mt.genrand_int31()%(np-np_A);       



		  zz_move[t1] = dcmplx(nfmod(real(zz_move[t1]) + (width)*(mt.random()-0.5),cut_sph) ,
                                     nfmod(imag(zz_move[t1]) + (width)*(mt.random()-0.5), 2.0*PI));


//		  zz_move[t1] = dcmplx(cut_sph + nfmod(real(zz_move[t1]) + (width)*(mt.random()-0.5),PI-cut_sph) ,
  //                                   nfmod(imag(zz_move[t1]) + (width)*(mt.random()-0.5), 2.0*PI));


//		  zz_move[t1] = dcmplx(nfmod(real(zz_move[t1]) + (width)*(mt.random()-0.5),PI) ,
  //                                   nfmod(imag(zz_move[t1]) + (width)*(mt.random()-0.5), 2.0*PI));




     
                  cut1=0;
                  for (i=0;i<np-np_A;i++)                
		      {
                       if (real(zz_move[i])>=0.0 && real(zz_move[i])<=cut_sph)
	               cut1++;                 
                       }
     
        ///////////////////////////////////////////////////////////




             if (cut1==np-np_A)
		 {

		   darea_b=1.0;
		   for (i=0; i<np-np_A; i++)
		        darea_b=darea_b*sin(real(zz_move[i]));


      			 for (i=0;i<np-np_A;i++)
			 {
				 zu_move[i]= cos(real(zz_move[i]/2.0))*exp(dcmplx(0,imag(zz_move[i])/2.0));
				 zv_move[i]= sin(real(zz_move[i]/2.0))*exp(dcmplx(0,-imag(zz_move[i])/2.0));
			 }
		

  
//               norm1_b=Jain_negative_2_LL_B (zu_move, zv_move);
                    norm1_b=IQHE_Prob_B(zu_move, zv_move);            

            
 
                   ratio=pow(abs(exp(norm1_b-norm1_a)),2.0)*(darea_b/darea_a);



				   
		   if ((ratio>=1.0 ) || ratio>mt.random())
			  {
			   accepted=accepted+1;
			   zz[t1]=zz_move[t1];
			   norm1_a=norm1_b;              
			   darea_a=darea_b;
			   }
		   else
		            rejected++;
	             }

	          else    
			    rejected++;
	}

//	printf("accepted %d, rejected %d\n",accepted, rejected);
	return;
}












void metropolis::run_Sphere_Laughlin(int nu1)
{ 
	   
	int i,j,i1,cut1, npart_red, t1[np];
	dcmplx norm1_a,norm1_b;
	double darea_a, darea_b, ratio;    



         accepted=0; rejected=0;        
	
	//	Convert zz to spinor co-ordiantes
	for (i=0;i<np;++i)
	{
		zu[i]= cos(real(zz[i]/2.0))*exp(dcmplx(0,imag(zz[i])/2.0));
		zv[i]= sin(real(zz[i]/2.0))*exp(dcmplx(0,-imag(zz[i])/2.0));
	}


          norm1_a=Laugh_Sphere(zu, zv, nu1);

 

	 darea_a=1.0;
	 for (i=0; i<np; i++)
             darea_a=darea_a*sin(real(zz[i]));



	for (i1=1;i1<=iter;i1++)
		{
		  for (i=0;i<np;i++)
		      zz_move[i]=zz[i];		  

                  

                      t1[0]=mt.genrand_int31()%np;

//                  for(i=0;i<np;++i)
           	      zz_move[t1[0]] = dcmplx(nfmod(real(zz_move[t1[0]]) + (width)*(mt.random()-0.5),PI) ,
                                     nfmod(imag(zz_move[t1[0]]) + (width)*(mt.random()-0.5), 2.0*PI));




		   darea_b=1.0;
		   for (i=0; i<np; i++)
		        darea_b=darea_b*sin(real(zz_move[i]));


      			 for (i=0;i<np;i++)
			 {
				 zu_move[i]= cos(real(zz_move[i]/2.0))*exp(dcmplx(0,imag(zz_move[i])/2.0));
				 zv_move[i]= sin(real(zz_move[i]/2.0))*exp(dcmplx(0,-imag(zz_move[i])/2.0));
			 }
		


                   norm1_b=Laugh_Sphere(zu_move, zv_move, nu1);


 
                   ratio=pow(abs(exp(norm1_b-norm1_a)),2.0)*(darea_b/darea_a);



				   
		   if ((ratio>=1.0 ) || ratio>mt.random())
			  {
			   accepted=accepted+1;
//                           for (i=0; i<np; ++i)
			       zz[t1[0]]=zz_move[t1[0]];
//			   zz[t2]=zz_move[t2];
			   norm1_a=norm1_b;              
			   darea_a=darea_b;
			   }
		   else
		            rejected++
;
	             }

	


	return;
}



















dcmplx metropolis::IQHE_Prob_A(dcmplx *zu1, dcmplx *zv1)
     {
      int i,j, nu_1;
      double QLaugh, neff;
      dcmplx expo, norm;


      QLaugh=(np-1.0)/2.0;

   
      norm=0.0;
      for (i=0;i<np_A-1;++i)
          for (j=i+1;j<np_A;++j)
              norm+=log(zu1[i]/zv1[i]-zu1[j]/zv1[j]);


      expo=0.0;
      for (i=0;i<np_A;++i)
          expo+=2.0*QLaugh*log(zv1[i]);


      nu_1=1;            

        norm=double(nu_1)*(norm+expo);

      return norm;
  
     }
 







dcmplx metropolis::IQHE_Prob_B(dcmplx *zu1, dcmplx *zv1)
     {
      int i,j,nu_1;
      double QLaugh, leff, neff;
      dcmplx expo, shift=0.0, norm;


      QLaugh=(np-1.0)/2.0;

   
      norm=0.0;
      for (i=0;i<(np-np_A)-1;++i)
          for (j=i+1;j<np-np_A;++j)
              norm+=log(zu1[i]/zv1[i]-zu1[j]/zv1[j]);


  

      shift=0.0;      //CHANGE
    
      for (i=0;i<np-np_A;++i)
          shift+=log(zu1[i]/zv1[i]);



      shift = (double)(np_A)*shift;
    

      expo=0.0;
      for (i=0;i<np-np_A;++i)
          expo+=2.0*QLaugh*log(zv1[i]);



      nu_1=1;            

        norm=double(nu_1)*(norm+expo+shift);

     return norm;
  
     }












dcmplx metropolis::Laugh_Sphere(dcmplx *zu1, dcmplx *zv1, int nu1)
{

      int i,j;
      double QLaugh;
      dcmplx expo, norm;

      QLaugh=(np-1.0)/2.0;


   
      norm=0.0;
      for (i=0;i<np-1;++i)
          for (j=i+1;j<np;++j)
              norm+=log(zu1[i]/zv1[i]-zu1[j]/zv1[j]);


      expo=0.0;
      for (i=0;i<np;++i)
          expo+=2.0*QLaugh*log(zv1[i]);

            
      norm=double(nu1)*(norm+expo);

     return norm;

}






dcmplx metropolis::EWF_Expansion_State(dcmplx *zu1, dcmplx *zv1, int npart_A)
{

         int lz,i,j, dim, total_num, Lz;
         double *Lz_EWF1, *Lz_EWF2, QJain, N1_Jain, N2_Jain, LLd, QJain_0th, QJain_1st, QJain_2nd, N_0th, N_1st, N_2nd, Lzz, Lz_0;
         dcmplx state_A, state_B, state_EWF=0.0, *zu_AA, *zv_AA, *zu_BB, *zv_BB;

         
         vector< vector <double> > Fspace, Fspace_comp;
    
         Wave_Functions *wf_EWF = new Wave_Functions(np, npart_A, nu, LL, 'P');




   


//////////////////////////////////// COMPUTE Lz0 (Lz0=Lz_min) ///////////////////////////////////////////////////////

          QJain=(np/2.0-2.0)/2.0; 
          N1_Jain=int((2.0*QJain+1.0)/2.0); 
          N2_Jain=npart_A-N1_Jain;
	  LLd=LL;
	  QJain_0th=(np/LLd-LLd)/2.0;
	  QJain_1st=(np/LLd-LLd+2)/2.0;
	  QJain_2nd=(np/LLd-LLd+4)/2.0;
	  N_0th=(2.0*QJain_0th+1.0)/2.0;
	  N_1st=(2.0*QJain_1st+1.0)/2.0;
	  N_2nd=npart_A-N_0th-N_1st;


	if (LL==1)
	   Lz_0=(npart_A*(npart_A-1)/2.0 - (np-1)*npart_A/2.0)*nu;   //FOR LAUGHLIN filling 1/nu
	else if (LL==2)
	        {
	         Lz_0 = N1_Jain*(N1_Jain-1.0)/2.0 - QJain*N1_Jain + N2_Jain*(N2_Jain-1.0)/2.0 - (QJain+1.0)*N2_Jain + 
	         (npart_A*(npart_A-1)/2.0 - (np-1)*npart_A/2.0)*nu;   // FOR JAIN nu=2/5 or nu=2/3
	          }
	else if (LL==3)
	        {
	         Lz_0 = N_0th*(N_0th-1.0)/2.0 - QJain_0th*N_0th + N_1st*(N_1st-1.0)/2.0 - QJain_1st*N_1st + N_2nd*(N_2nd-1.0)/2.0 
		- QJain_2nd*N_2nd + (npart_A*(npart_A-1)/2.0 - (np-1)*npart_A/2.0)*nu;   // FOR JAIN nu=3/7 or nu=3/5
	          }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////




         
        zu_AA=new dcmplx[np];
        zv_AA=new dcmplx[np];
        zu_BB=new dcmplx[np];
        zv_BB=new dcmplx[np];


        for (i=0; i<npart_A; ++i)
            {
             zu_AA[i]=zu1[i];
             zv_AA[i]=zv1[i];
            }

        for (i=npart_A; i<np; ++i)
            {
             zu_BB[i-npart_A]=zu1[i];
             zv_BB[i-npart_A]=zv1[i];
             }





      if (LL==1)
         total_num=nu*npart_A+LL+1;
      else 
          total_num=nu*npart_A+npart_A+LL+1;
          



         Lz_EWF1=new double [nu*np+np+LL+4]();
         Lz_EWF2=new double [nu*np+np+LL+4]();


 

	  for (Lz=0; Lz<=-2*Lz_0; ++Lz)
	      {
               Lzz=Lz+Lz_0;



	       EWF_Space(Lzz, np, dim, npart_A, LL, Fspace, Fspace_comp, nu);


 
        for (i=0; i<dim; ++i)
            {

            for (j=0;j<total_num; ++j)
                {
                 Lz_EWF1[j]=Fspace[i][j];
                 Lz_EWF2[j]=Fspace_comp[i][j];
                 }
            Lz_EWF2[total_num-1]=1;

            wf_EWF-> w_func (state_A, Lz_EWF1, np, zu_AA, zv_AA, binom, nu, LL);
            wf_EWF-> w_func (state_B, Lz_EWF2, np, zu_BB, zv_BB, binom, nu, LL);

            
            
            state_EWF+=exp(state_A)*exp(state_B);

            


             }
      }





delete [] Lz_EWF1; delete [] Lz_EWF2; delete [] zu_AA; delete [] zv_AA; delete [] zu_BB; delete [] zv_BB;

return state_EWF;

wf_EWF-> ~Wave_Functions();

}                







dcmplx metropolis::Jain_negative_2_LL_A (dcmplx *zu1, dcmplx *zv1)
     {

     int i,j,k, dim_M1, n1, n2;
     dcmplx FJ[20][20][20],PJ[20][20][20],laugh,expo,wfuncCF,wfuncCF1,wfuncCF2, norm1, norm;
     double p,QJain,QLaugh, Lz_EWF1[20];

     p=1.0;



 Lz_EWF1[0]=2;
 Lz_EWF1[1]=1;
 Lz_EWF1[2]=-0.5; Lz_EWF1[3]=0.5; Lz_EWF1[4]=1.5;

  

 n1=Lz_EWF1[0];
 n2=Lz_EWF1[1];
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
              FJ[1][0][i]= FJ[1][0][i] + pow( (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[0][1][i]= FJ[0][1][i] + pow( (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);

              FJ[2][0][i]= FJ[2][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[1][1][i]= FJ[1][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[0][2][i]= FJ[0][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);


              FJ[3][0][i]= FJ[3][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[2][1][i]= FJ[2][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[1][2][i]= FJ[1][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[0][3][i]= FJ[0][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);


              FJ[4][0][i]= FJ[4][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[3][1][i]= FJ[3][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[2][2][i]= FJ[2][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[1][3][i]= FJ[1][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);
              FJ[0][4][i]= FJ[0][4][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0);


              FJ[5][0][i]= FJ[5][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),5.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[4][1][i]= FJ[4][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[3][2][i]= FJ[3][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[2][3][i]= FJ[2][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);
              FJ[1][4][i]= FJ[1][4][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0);
              FJ[0][5][i]= FJ[0][5][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),5.0); 
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





       QJain=(np/2.0 - 2.0)/2.0;




        for (i=0;i<n1;i++)
            for (j=0;j<dim_M1;j++)
                 Matrix1[i][j]= PJ[(int) (QJain+Lz_EWF1[i+2])][(int) (QJain-Lz_EWF1[i+2])][j];  
                

                    
        for (i=n1;i<dim_M1;i++)
            for (j=0;j<dim_M1;j++)
                 {
                  Matrix1[i][j]= binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF1[i+2])]*zv1[j]*
                                PJ[(int) (QJain+Lz_EWF1[i+2])][(int) (QJain-Lz_EWF1[i+2]+1)][j] - 
                                    binom[(int) (2.0*QJain+1.0)][(int) (QJain-Lz_EWF1[i+2])]*zu1[j]*
                                PJ[(int) (QJain+Lz_EWF1[i+2]+1)][(int) (QJain-Lz_EWF1[i+2])][j];    

    //             cout<<(int) (QJain+Lz_EWF1[i+2]) <<" "<<(int) (QJain-Lz_EWF1[i+2]+1)<<endl;
      //           cout<<binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF1[i+2])]<<endl;
             //    cin.get();
     
                 } 






                        norm=LogDeterminant(Matrix1,dim_M1);
          //   norm=LogDeterminant<dcmplx>(Matrix1,dim_M1);





      QLaugh=(np-1.0)/2.0;

   
      norm1=0.0;
      for (i=0;i<dim_M1-1;++i)
          for (j=i+1;j<dim_M1;++j)
              norm1+=log(zu1[i]/zv1[i]-zu1[j]/zv1[j]);


      expo=0.0;
      for (i=0;i<dim_M1;++i)
          expo+=2.0*QLaugh*log(zv1[i]);



        norm1=(norm1+expo);


    norm+=norm1;


     return norm;

      }













dcmplx metropolis::Jain_negative_2_LL_B (dcmplx *zu1, dcmplx *zv1)
     {

     int i,j,k, dim_M1, n1, n2;
     dcmplx FJ[20][20][20],PJ[20][20][20],laugh,expo,wfuncCF,wfuncCF1,wfuncCF2, norm1, norm;
     double p,QJain,QLaugh, Lz_EWF1[20];


 p=1.0;


 Lz_EWF1[0]=2;
 Lz_EWF1[1]=1;
 Lz_EWF1[2]=-0.5; Lz_EWF1[3]=0.5; Lz_EWF1[4]=1.5;

  

 n1=Lz_EWF1[0];
 n2=Lz_EWF1[1];
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
              FJ[1][0][i]= FJ[1][0][i] + pow( (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[0][1][i]= FJ[0][1][i] + pow( (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);


              FJ[2][0][i]= FJ[2][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[1][1][i]= FJ[1][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[0][2][i]= FJ[0][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);


              FJ[3][0][i]= FJ[3][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[2][1][i]= FJ[2][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[1][2][i]= FJ[1][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[0][3][i]= FJ[0][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);


              FJ[4][0][i]= FJ[4][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[3][1][i]= FJ[3][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[2][2][i]= FJ[2][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[1][3][i]= FJ[1][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);
              FJ[0][4][i]= FJ[0][4][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0);


              FJ[5][0][i]= FJ[5][0][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),5.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0);
              FJ[4][1][i]= FJ[4][1][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0);
              FJ[3][2][i]= FJ[3][2][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0);
              FJ[2][3][i]= FJ[2][3][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),2.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),3.0);
              FJ[1][4][i]= FJ[1][4][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),1.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),4.0);
              FJ[0][5][i]= FJ[0][5][i] + pow(  (zv1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),0.0) *  pow( (-zu1[k]/(zu1[i]*zv1[k]-zv1[i]*zu1[k])),5.0); 
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





       QJain=(np/2.0 - 2.0)/2.0;




        for (i=0;i<n1;i++)
            for (j=0;j<dim_M1;j++)
                 Matrix1[i][j]= PJ[(int) (QJain+Lz_EWF1[i+2])][(int) (QJain-Lz_EWF1[i+2])][j];  
                

                    
        for (i=n1;i<dim_M1;i++)
            for (j=0;j<dim_M1;j++)
                 {
                  Matrix1[i][j]= binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF1[i+2])]*zv1[j]*
                                PJ[(int) (QJain+Lz_EWF1[i+2])][(int) (QJain-Lz_EWF1[i+2]+1)][j] - 
                                    binom[(int) (2.0*QJain+1.0)][(int) (QJain-Lz_EWF1[i+2])]*zu1[j]*
                                PJ[(int) (QJain+Lz_EWF1[i+2]+1)][(int) (QJain-Lz_EWF1[i+2])][j];    

    //             cout<<(int) (QJain+Lz_EWF1[i+2]) <<" "<<(int) (QJain-Lz_EWF1[i+2]+1)<<endl;
      //           cout<<binom[(int) (2.0*QJain+1.0)][(int) (QJain+1.0-Lz_EWF1[i+2])]<<endl;
             //    cin.get();
     
                 } 






                      norm=LogDeterminant(Matrix1,dim_M1);
         //    norm=LogDeterminant<dcmplx>(Matrix1,dim_M1);







      QLaugh=(np-1.0)/2.0;



     dcmplx shift=0.0;      //CHANGE
    
      for (i=0;i<dim_M1;++i)
          shift+=log(zu1[i]/zv1[i]);


      shift = (double)((np-dim_M1))*shift;

   
      norm1=0.0;
      for (i=0;i<dim_M1-1;++i)
          for (j=i+1;j<dim_M1;++j)
              norm1+=log(zu1[i]/zv1[i]-zu1[j]/zv1[j]);


      expo=0.0;
      for (i=0;i<dim_M1;++i)
          expo+=2.0*QLaugh*log(zv1[i]);



        norm1=(norm1+expo+shift);


    norm+=norm1;


    return norm;

      }







 










vector<double>  metropolis::correl(vector<double> Obs, int num_iter, double mean, int range)
      {

      int i,l;
      double sum_corr, norm_corr;
      vector<double> Corr;

      for (i=0; i<num_iter; ++i)
          Obs[i]=Obs[i]-mean;
      

      for(l=0;l<range;++l)
         {
          sum_corr=0.0;
          for (i=0;i<num_iter-l;++i)
              sum_corr+=(Obs[i]*Obs[i+l])*1.0/(1.0*num_iter-1.0*l);
          Corr.push_back(sum_corr);
          }


      norm_corr=Corr[0];
      for(l=0;l<range;++l)
          Corr[l]=Corr[l]/norm_corr;

      return Corr;


      }





void metropolis::plot(vector<double> Corr, int range)
{

    vector<double> gnu_x, gnu_y;
    int i;

    cout << "*** start of GNUPLOT ***" << endl;

    try {
         Gnuplot g1 = Gnuplot("lines");

    for (i = 0; i < range; i++)
      {
        gnu_x.push_back((double)i);
        gnu_y.push_back(Corr[i]);
  //    gnu_y.push_back(Obs[i]);
      }


    g1.set_style("linespoints");
    g1.plot_xy(gnu_x,gnu_y,"Metropolis Correlations.");
    g1.set_style("lines");
    g1.plot_equation("exp(-1)","1/e");



    cin.get();
     g1.reset_plot();


     } catch (GnuplotException ge) {
         cout << ge.what() << endl;
     }

     cout << endl << "*** end of GNUPLOT ***** " << endl;

}









void metropolis::run_Sphere()
{ 
	   
	int i,j,i1,cut1, npart_red, t1[np];
	dcmplx norm1_a,norm1_b;
	double darea_a, darea_b, ratio;    



         accepted=0; rejected=0;        
	
	//	Convert zz to spinor co-ordiantes
	for (i=0;i<np;++i)
	{
		zu[i]= cos(real(zz[i]/2.0))*exp(dcmplx(0,imag(zz[i])/2.0));
		zv[i]= sin(real(zz[i]/2.0))*exp(dcmplx(0,-imag(zz[i])/2.0));
	}


  
          wf->w_func (norm1_a, Lz_EWF, np, zu, zv, binom, nu, LL);

//           norm1_a = log(EWF_Expansion_State(zu, zv, np/2));




 

 

	 darea_a=1.0;
	 for (i=0; i<np; i++)
             darea_a=darea_a*sin(real(zz[i]));



	for (i1=1;i1<=iter;i1++)
		{
		  for (i=0;i<np;i++)
		      zz_move[i]=zz[i];		  

                  

                      t1[0]=mt.genrand_int31()%np;

//                  for(i=0;i<np;++i)
           	      zz_move[t1[0]] = dcmplx(nfmod(real(zz_move[t1[0]]) + (width)*(mt.random()-0.5),PI) ,
                                     nfmod(imag(zz_move[t1[0]]) + (width)*(mt.random()-0.5), 2.0*PI));




		   darea_b=1.0;
		   for (i=0; i<np; i++)
		        darea_b=darea_b*sin(real(zz_move[i]));


      			 for (i=0;i<np;i++)
			 {
				 zu_move[i]= cos(real(zz_move[i]/2.0))*exp(dcmplx(0,imag(zz_move[i])/2.0));
				 zv_move[i]= sin(real(zz_move[i]/2.0))*exp(dcmplx(0,-imag(zz_move[i])/2.0));
			 }
		


                 wf->w_func (norm1_b, Lz_EWF, np, zu_move, zv_move, binom, nu, LL);
//                   norm1_b = log(EWF_Expansion_State(zu_move, zv_move, np/2));



                   ratio=pow(abs(exp(norm1_b-norm1_a)),2.0)*(darea_b/darea_a);



				   
		   if ((ratio>=1.0 ) || ratio>mt.random())
			  {
			   accepted=accepted+1;
//                           for (i=0; i<np; ++i)
			       zz[t1[0]]=zz_move[t1[0]];
//			   zz[t2]=zz_move[t2];
			   norm1_a=norm1_b;              
			   darea_a=darea_b;
			   }
		   else
		            rejected++
;
	             }

	


	return;
}








//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// CALCULATION OF NORM USING NU=1 AS PROB. DISTRIBUTION IN METROPOLIS //////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////




dcmplx metropolis::Wfunction_Norm (double width_norm, int iter_norm, int init_norm, int metrop_norm, dcmplx &alpha)
{

double norm_nu_1, norm_state,  norma_state=0.0, norma_state_sq=0.0, norma_state_aux=0.0, norma_state_sq_aux=0.0, error, 
error_block, norm_Gstate, norm_Error;
int i,j,k, l, ne, num_M, num_L, rank, numtasks, aux=0;
dcmplx state_0, state_A, sign, *zu_A, *zv_A;
char flux='P';


MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);


ne=np;
width=width_norm;

//////////////////////// log of norm of nu=1 IQH state used in metropolis_sphere //////////////////

norm_nu_1=0.0;
for (i=0;i<ne;++i)
     norm_nu_1+=log(4.0*PI*Gamma((double)(1.0+i))*Gamma((double)(ne-i))/Gamma((double)(ne+1.0)));

norm_nu_1+=log(factorial(ne));


///////////////////////////////////////////////////////////////////////////////////////////////////



zu_A=new dcmplx [ne]();
zv_A=new dcmplx [ne]();


///////////////////// COORDINATE INITIALIZATION FOR METROPOLIS //////////////////////////////////

for (i=0;i<ne;i++)
    zz[i]=dcmplx((PI*(i+1)/ne), nfmod(mt.random(),2.0*PI));

//////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////// INITIALIZE METROPOLIS ////////////////////////////////////////////


for (k=1; k<=init_norm ; k++)
    {
     run_Sphere_Laughlin(1);    

     if (rank==0 && k%400==0)
         cout<<"Norm initial iteration "<<k<<" accepted "<<accepted<<" rejected "<<rejected<<endl;



     equilibration (accepted, rejected);

    }

//////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////// COMPUTE SCALING ALPHA TO MAKE NORM NUMERICALLY CONVERGENT //////////////////





if (rank==0)
{

for (k=1; k<=4; k++)
    {

     run_Sphere_Laughlin(1);


           for (j=0;j<ne;++j)
	       {
	        zu_A[j]= cos(real(zz[j]/2.0))*exp(dcmplx(0,imag(zz[j])/2.0));
	        zv_A[j]= sin(real(zz[j]/2.0))*exp(dcmplx(0,-imag(zz[j])/2.0));
      	        }


          state_0 = Laugh_Sphere(zu_A, zv_A, 1);


          wf-> w_func (state_A, Lz_EWF, ne, zu_A, zv_A, binom, nu, LL);
//             state_A=log(EWF_Expansion_State(zu_A, zv_A, np/2)); 

          alpha = norm_nu_1 + state_A + conj(state_A) - state_0 - conj(state_0) ;

          }

}




MPI_Bcast(&alpha, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

//alpha=0.0;

////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////// NORM COMPUTATION BY MONTECARLO  //////////////////////////////////



num_M=iter_norm;
num_L=metrop_norm/num_M;

double *error_pointer = new double[metrop_norm + 1]();
double *means_pointer = new double[num_L + 1]();


for (k=1; k<=metrop_norm; k++)
    {

    run_Sphere_Laughlin(1);

     if (rank==0 && k%500==0)
         cout<<"Norm metrop. iteration: "<<k<<" of  "<<metrop_norm<<endl;
  

           for (j=0;j<ne;++j)
	       {
	        zu_A[j]= cos(real(zz[j]/2.0))*exp(dcmplx(0,imag(zz[j])/2.0));
	        zv_A[j]= sin(real(zz[j]/2.0))*exp(dcmplx(0,-imag(zz[j])/2.0));              
      	        }


           state_0 = Laugh_Sphere(zu_A, zv_A,1);



               wf-> w_func (state_A, Lz_EWF, ne, zu_A, zv_A, binom, nu, LL);                  
    //         cout<<state_A<<endl;
    //         state_A=log(EWF_Expansion_State(zu_A, zv_A, np/2)); 
    //         cout<<state_A<<endl;
    //         state_A=log(EWF_Expansion_State(zu_A, zv_A, np/2-1)); 
    //         cout<<state_A<<endl;
     //        cin.get();




           error_pointer[k-1] = exp(real(norm_nu_1 + state_A + conj(state_A) - state_0 - conj(state_0)  - alpha  )) ;
                       
    }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// error computed by brute force  (valid with independent measures)//////
//norma=norma_state*exp(real(alpha));
//error=sqrt((exp(2.0*real(alpha))*norma_state_sq-pow(exp(real(alpha))*norma_state,2.0))/(1.0*metrop_norm-1.0));

//cout<<"first method:  "<<norma<<" , "<<error<<endl;
//cout<<"relative error:  "<<1.0*error/(1.0*norma)<<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////// error computed with blocking method //////////////////////////////////


for (i=0;i<metrop_norm;++i)
     {
      if ( ((i+1) % num_M) == 0)
          aux=aux+1;
     means_pointer[aux]+=error_pointer[i]/(1.0*num_M);
     }


//for (i=0; i<num_L; ++i)
//     means_pointer[i]*=exp(real(alpha));




for (i=0;i<num_L;++i)
    {
     norma_state+=means_pointer[i]/(1.0*num_L);
     norma_state_sq+=pow(means_pointer[i],2.0)/(1.0*num_L);
    }

      
  MPI_Reduce (&norma_state, &norma_state_aux,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce (&norma_state_sq, &norma_state_sq_aux,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     
    
      if (rank==0)
         {

          norma_state = norma_state_aux/(1.0*numtasks);
          norma_state_sq = norma_state_sq_aux/(1.0*numtasks);
         
          norm_Gstate=norma_state;
          norm_Error=sqrt( (norma_state_sq - pow(norma_state,2.0))/(1.0*num_L*numtasks -1.0) );

         }


delete [] zu_A; delete [] zv_A; delete [] error_pointer; delete [] means_pointer;



MPI_Bcast(&norm_Gstate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
MPI_Bcast(&norm_Error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


MPI_Barrier(MPI_COMM_WORLD);


return dcmplx(norm_Gstate, norm_Error); 


}

                                 






void metropolis::Wfunction_Corr_IQHE_N(int num_metrop_corr, int range)

    {


     vector <double> Obs, Corr;
     double sum_corr, mean_zz;
     int i;

     int iter_aux=iter;   /// To use nu=1 for prob. distribution in run_sphere_Laughlin()  /////
     

     iter=1;    ///// only for correlations. I restore variable at the end of function //////

//////////////////////////////////////// COORDINATES INITIALIZATION //////////////////////////////

     for (i=0;i<np;i++)
          zz[i]=dcmplx((PI*(i+1)/np), nfmod(mt.random(),2.0*PI));


//////////////////////////// METROPOLIS CORRELATION CALCULATIONS ////////////////////////////////



     for (i=0;i<num_metrop_corr;++i)
        {
         run_Sphere_Laughlin(1);
         cout<<"accepted, rejected:  "<<accepted<<"  "<<rejected<<endl;    
 
         Obs.push_back(real(zz[0]));

         mean_zz+=pow(real(zz[0]),1.0)/(1.0*num_metrop_corr);

         }


Corr=correl(Obs,num_metrop_corr,mean_zz,range);

plot(Corr, range);


iter=iter_aux;

return;


}








void metropolis::Wfunction_Corr_IQHE_NA(int num_metrop_corr, int range)

    {


     vector <double> Obs, Corr;
     double sum_corr, mean_zz;
     int i;
    

//////////////////////////////////////// COORDINATES INITIALIZATION //////////////////////////////

    for (i=0;i<np_A;i++)
          zz[i]=dcmplx(((PI-0.02)*i/(1.0*np_A)) + cut_sph + 0.01, nfmod(mt.random(),2.0*PI));

//////////////////////////// METROPOLIS CORRELATION CALCULATIONS ////////////////////////////////



     for (i=0;i<num_metrop_corr;++i)
        {
         run_A();
         cout<<"accepted, rejected:  "<<accepted<<"  "<<rejected<<endl;    
 
         Obs.push_back(real(zz[0]));

         mean_zz+=pow(real(zz[0]),1.0)/(1.0*num_metrop_corr);

         }


Corr=correl(Obs,num_metrop_corr,mean_zz,range);

plot(Corr, range);


return;


}








void metropolis::Wfunction_Corr_IQHE_NB(int num_metrop_corr, int range)

    {


     vector <double> Obs, Corr;
     double sum_corr, mean_zz;
     int i;
    

//////////////////////////////////////// COORDINATES INITIALIZATION //////////////////////////////

     for (i=0;i<np-np_A;i++)
          zz[i]=dcmplx((cut_sph-0.02)*i/(1.0*np-1.0*np_A) + 0.01, nfmod(mt.random(),2.0*PI));

//////////////////////////// METROPOLIS CORRELATION CALCULATIONS ////////////////////////////////



     for (i=0;i<num_metrop_corr;++i)
        {
         run_B();
         cout<<"accepted, rejected:  "<<accepted<<"  "<<rejected<<endl;    
 
         Obs.push_back(real(zz[0]));

         mean_zz+=pow(real(zz[0]),1.0)/(1.0*num_metrop_corr);

         }


Corr=correl(Obs,num_metrop_corr,mean_zz,range);

plot(Corr, range);


return;


}




void metropolis::Wfunction_Corr_N(int num_metrop_corr, int range)

    {


     vector <double> Obs, Corr;
     double sum_corr, mean_zz;
     int i;
    

//////////////////////////////////////// COORDINATES INITIALIZATION //////////////////////////////

     for (i=0;i<np;i++)
          zz[i]=dcmplx(((cut_sph-0.02)*i/(1.0*np)) + 0.01, nfmod(mt.random(),2.0*PI));

//////////////////////////// METROPOLIS CORRELATION CALCULATIONS ////////////////////////////////



     for (i=0;i<num_metrop_corr;++i)
        {
         run_Sphere();
         cout<<"accepted, rejected:  "<<accepted<<"  "<<rejected<<endl;    
 
         Obs.push_back(real(zz[0]));

         mean_zz+=pow(real(zz[0]),1.0)/(1.0*num_metrop_corr);

         }


Corr=correl(Obs,num_metrop_corr,mean_zz,range);

plot(Corr, range);


return;


}










dcmplx metropolis::Wfunction_Coulomb_energy (int initial, int num_metrop_iter)
     {

      int i, j, k;
      double radio_sp[np][np], Coulomb=0.0, QLaugh, Coulomb_sq=0.0;

      QLaugh=nu*(np -1.0)/2.0;



      for (i=0;i<initial;++i)
         { 
          run_Sphere();
          
          cout<<"Coulomb initial:  "<<accepted<<" "<<rejected<<endl;

//          equilibration(accepted, rejected);

           }





      for (k=0;k<num_metrop_iter;++k)
         { 
          run_Sphere();
          
          cout<<"Coulomb calculation:  "<<accepted<<" "<<rejected<<endl;


               for(i=0;i<np-1;++i)             // For Energy
                for(j=i+1;j<np;++j)
                      radio_sp[i][j]=2.0*sqrt(QLaugh)*abs( cos(real(zz[i]/2.0))*sin(real(zz[j]/2.0)) - 
                      cos(real(zz[j]/2.0))*sin(real(zz[i]/2.0))*exp(dcmplx(0.0,imag(zz[i]-zz[j]))) ); 

              double interaction=0.0;                   // For Energy
              for(i=0;i<np-1;i++)
                 for(j=i+1;j<np;j++)
                     interaction+=(1.0/(1.0*np))*1.0/radio_sp[i][j];    

             Coulomb+=interaction/num_metrop_iter;  
 
             Coulomb_sq+=pow(interaction,2.0)/num_metrop_iter;  


           }


      cout<< "Coulomb Energy Unnomrmalized:  "<<Coulomb<<endl;

      return dcmplx (Coulomb - 0.5*np/sqrt(QLaugh),sqrt( (Coulomb_sq - pow(Coulomb, 2.0))/num_metrop_iter ) );


      }


























void metropolis::Coulomb_Thermal (double width_norm, int metrop_norm, int range)
{

double radio_sp[np][np], Coulomb=0.0, QLaugh;
int i,j,k, l, ne;
char flux='P';
vector <double> Obs;




ne=np;
width=width_norm;


///////////////////// COORDINATE INITIALIZATION FOR METROPOLIS //////////////////////////////////

for (i=0;i<ne;i++)
    zz[i]=dcmplx((PI*(i+1)/ne), nfmod(mt.random(),2.0*PI));

/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////  MONTECARLO THERMALIZATION  /////////////////////////////////////

int iter_aux=iter;
iter=1;


     QLaugh=nu*(np -1.0)/2.0;



      for (k=0;k<metrop_norm;++k)
         { 
          run_Sphere();
          
          cout<<"Energy Thermalization:  "<<accepted<<" "<<rejected<<endl;


               for(i=0;i<np-1;++i)           
                for(j=i+1;j<np;++j)
                      radio_sp[i][j]=2.0*sqrt(QLaugh)*abs( cos(real(zz[i]/2.0))*sin(real(zz[j]/2.0)) - 
                      cos(real(zz[j]/2.0))*sin(real(zz[i]/2.0))*exp(dcmplx(0.0,imag(zz[i]-zz[j]))) ); 

              double interaction=0.0;        
              for(i=0;i<np-1;i++)
                 for(j=i+1;j<np;j++)
                     interaction+=(1.0/(1.0*np))*1.0/radio_sp[i][j];    

             Coulomb+=interaction/metrop_norm;  

             Obs.push_back(interaction - 0.5*np/sqrt(QLaugh));


           }



plot(Obs, range);


iter=iter_aux;




return; 

}























void metropolis::Wfunction_Norm_Thermal (double width_norm, int metrop_norm, int range)
{

double norm_nu_1, norm_state, alpha;
int i,j,k, l, ne;
dcmplx state_0, state_A, sign, *zu_A, *zv_A;
char flux='P';
vector <double> Obs;




ne=np;
width=width_norm;

//////////////////////// log of norm of nu=1 IQH state used in metropolis_sphere //////////////////

norm_nu_1=0.0;
for (i=0;i<ne;++i)
     norm_nu_1+=log(4.0*PI*Gamma((double)(1.0+i))*Gamma((double)(ne-i))/Gamma((double)(ne+1.0)));

norm_nu_1+=log(factorial(ne));

///////////////////////////////////////////////////////////////////////////////////////////////////



zu_A=new dcmplx [ne]();
zv_A=new dcmplx [ne]();


///////////////////// COORDINATE INITIALIZATION FOR METROPOLIS //////////////////////////////////

for (i=0;i<ne;i++)
    zz[i]=dcmplx((PI*(i+1)/ne), nfmod(mt.random(),2.0*PI));

/////////////////////////////////////////////////////////////////////////////////////////////////


int iter_aux=iter;   

iter=1;     //// only for thermalization. I restore this variable to its original value at the end of the function ///////



///////////////////////// COMPUTE SCALING ALPHA TO MAKE NORM NUMERICALLY CONVERGENT //////////////////



for (k=1; k<=4; k++)
    {

     run_Sphere_Laughlin(1);


           for (j=0;j<ne;++j)
	       {
	        zu_A[j]= cos(real(zz[j]/2.0))*exp(dcmplx(0,imag(zz[j])/2.0));
	        zv_A[j]= sin(real(zz[j]/2.0))*exp(dcmplx(0,-imag(zz[j])/2.0));
      	        }


          state_0 = Laugh_Sphere(zu_A, zv_A, 1);

          wf-> w_func (state_A, Lz_EWF, ne, zu_A, zv_A, binom, nu, LL);

          alpha = real(norm_nu_1 + state_A + conj(state_A) - state_0 - conj(state_0)) ;


          }


//alpha=0.0;

////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////  MONTECARLO THERMALIZATION  /////////////////////////////////////



for (k=1; k<=metrop_norm; k++)
//  for (k=1; k<=10000; k++)
    {

    run_Sphere_Laughlin(1);

    cout<<"metrop: "<<k<<" accepted,  rejected  :"<<accepted<<" "<<rejected<<endl;


           for (j=0;j<ne;++j)
	       {
	        zu_A[j]= cos(real(zz[j]/2.0))*exp(dcmplx(0,imag(zz[j])/2.0));
	        zv_A[j]= sin(real(zz[j]/2.0))*exp(dcmplx(0,-imag(zz[j])/2.0));              
      	        }


           state_0 = Laugh_Sphere(zu_A, zv_A, 1);


           wf-> w_func (state_A, Lz_EWF, ne, zu_A, zv_A, binom, nu, LL);


           Obs.push_back(real(exp(norm_nu_1 + state_A + conj(state_A) - state_0 - conj(state_0) -alpha )));

//          cout<<"Obs:  "<<Obs[k-1]<<endl;
  //        cout<<"Obs:  "<<pow(10.0,4.0)*exp(real(norm_nu_1 + state_A + conj(state_A) - state_0 - conj(state_0) -alpha ))<<endl;
    //      cin.get();

          }

iter=iter_aux;


plot(Obs, range);



delete [] zu_A; delete [] zv_A;




return; 

}







void metropolis::equilibration(int accep, int rejec)
     {

      if ((accep-rejec)<6)
         width+=0.1;
      else if ((accep-rejec)>=6)
               width-=0.1;
        
      }




















/*
          for(i=0;i<ne-1;i++)             // For Energy
                for(j=i+1;j<ne;j++)
                    {radio_sp[i][j]=2.0*sqrt(QLaugh)*abs( cos(real(zz[i]/2.0))*sin(real(zz[j]/2.0)) - 
                                    cos(real(zz[j]/2.0))*sin(real(zz[i]/2.0))*exp(dcmplx(0.0,imag(zz[i]-zz[j]))) );} 

              double interaction=0.0;                   // For Energy
              for(i=0;i<ne-1;i++)
                 for(j=i+1;j<ne;j++)
                     interaction+=(1.0/ne)*1.0/radio_sp[i][j];    




         Coul+=norm_nu_1*interaction*(abs(exp(state_A)*conj(exp(state_A))))/(abs(exp(state_0)*conj(exp(state_0)))*
               (double)metrop_norm);         */


/*     metropolis_Sphere(ne, zz, width_norm, iter_norm, nu);


               for(i=0;i<ne-1;i++)             // For Energy
                for(j=i+1;j<ne;j++)
                    {radio_sp[i][j]=2.0*sqrt(QLaugh)*abs( cos(real(zz[i]/2.0))*sin(real(zz[j]/2.0)) - 
                                    cos(real(zz[j]/2.0))*sin(real(zz[i]/2.0))*exp(dcmplx(0.0,imag(zz[i]-zz[j]))) );} 

              double interaction=0.0;                   // For Energy
              for(i=0;i<ne-1;i++)
                 for(j=i+1;j<ne;j++)
                     interaction+=(1.0/ne)*1.0/radio_sp[i][j];    

             Coulomb+=interaction/metrop_norm;  */


