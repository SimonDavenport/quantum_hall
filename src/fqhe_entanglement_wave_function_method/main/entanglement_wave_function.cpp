#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <memory.h>
#include <complex>
#include <fstream>      //  for the std::fstream variable type
#include <sstream>      //  for the std::stringstream variable type
#include <unistd.h>
#include <ctime>
#include "../gnuplot/gnuplot_i.hpp"   /// gnuplot
#include <string>
#include <algorithm>

#include "../../utilities/mathematics/mt.hpp"
#include "../wave_functions/wave_functions.hpp"
#include "../hilbert_ewf/hilbert_ewf.hpp"	
#include "../algebra/algebra.hpp"
//#include "Algebra/determinant.hpp"
#include "../algebra/deter.hpp"  
#include "../metropolis/metropolis.hpp"
//#include <boost/math/special_functions/gamma.hpp>










//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////PARAMETERS OPTIONS/////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////


 
 
#include <boost/program_options.hpp>
#include <sstream>
 
struct Params
{
    std::string path;
    int ne;    
    int npart_i;
    int npart_f;
    int npart;
    int nu;   
    int LL;
    int num_Lz_sectors;
    int num_iter_metrop;
    int metrop_initial;
    int iter;
    int iter_normal;
    int num_L;
    int initial_normal;
    char plot;



     
    //  Constructor from command line parsing
     
    Params(int argc, char* argv[])
    {
        namespace po = boost::program_options;
 
        po::variables_map vm;
     
        int rank;
        bool exitFlag = false;
        path = "";
   
     
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         
        if(rank == 0)   // FOR THE MASTER NODE
        {      
            //  General options
            po::options_description general_opt("General Options");
            general_opt.add_options()
            ("help,h",
             "Display this message\n")
            ("path,p",po::value<std::string>()->default_value("./Eigen/"),
             "Set path to output files\n")
            ("nbr-n,n",po::value<int>()->default_value(12),
             "Set number of particles\n")
            ("nbr-i,i",po::value<int>()->default_value(6),
             "Set number of particles in the A subsystem (default 6). Set to -1 for all NA sectors.\n")
            ("nbr-f,f",po::value<int>()->default_value(6),
             "Set number of particles in the A subsystem (default 6). Set to -1 for all NA sectors\n")

            ("nu",po::value<int>()->default_value(2),
             " Laughlin state = prod_{i<j} (z_i - z_j)^nu  -- Jain state = IQH state(with filling = LL) * prod_{i<j} (z_i - z_j)^nu \n")
            ("LL",po::value<int>()->default_value(1),
             "Number of effective Landau levels\n")
            ("sector,s",po::value<int>()->default_value(8),
             "number of Lz_A sectors. Set to -1 for all LzA sectors\n")

            ("mp-inter-rho",po::value<int>()->default_value(40),
             "Metropolis internal iterations for P and Q matrices\n")
            ("mp-iter-rho",po::value<int>()->default_value(40),
             "Metropolis iterations for P and Q matrices = mp-inter-rho x mp-iter-rho\n")
            ("mp-therm-rho",po::value<int>()->default_value(200),
             "Metropolis thermalization for Q and P matrices\n")
            ("mp-inter-gs",po::value<int>()->default_value(40),
             "Metropolis internal iterations for g.state norm\n")
            ("mp-iter-gs",po::value<int>()->default_value(40),
             "Metropolis iterations for g.state norm = mp-inter-gs x mp-iter-gs\n")
            ("mp-therm-gs",po::value<int>()->default_value(200),
             "Metropolis thermalization for g.state norm\n")

            ("plot,p",po::value<char>()->default_value('N'),
             "y: to save plot (default = N) \n");

 
            //  Map the command line argument onto the option set
            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv,general_opt), vm);
            po::notify(vm);
 
            //  Respond to help request
            if(vm.count("help"))
            {
                std::cout << general_opt << "\n";
                exitFlag=true;
            }


 
            path             = vm["path"].as<std::string>();
            ne               = vm["nbr-n"].as<int>();

            npart_i          = vm["nbr-i"].as<int>();
            npart_f          = vm["nbr-f"].as<int>();

            nu               = vm["nu"].as<int>();
            LL               = vm["LL"].as<int>();
            num_Lz_sectors   = vm["sector"].as<int>();

            iter             = vm["mp-inter-rho"].as<int>();   
            num_iter_metrop  = vm["mp-iter-rho"].as<int>();   
            metrop_initial   = vm["mp-therm-rho"].as<int>();

            iter_normal  = vm["mp-inter-gs"].as<int>();   
            num_L  = vm["mp-iter-gs"].as<int>();   
            initial_normal   = vm["mp-therm-gs"].as<int>();

            plot = vm["plot"].as<char>();




 
 
        }
         
        MPI_Bcast(&exitFlag,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
         
        if(exitFlag)
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
         
        //  Print parameters
        if(0 == rank)   // FOR THE MASTER NODE
        {
            std::cout<<"path to save Eigenvalues\t\t\t\t= "<<path<<std::endl;
            std::cout<<"ne \t\t\t\t\t\t\t= "<<ne<<std::endl;
            std::cout<<"NA initial \t\t\t\t\t\t= "<<npart_i<<std::endl;
            std::cout<<"NA final \t\t\t\t\t\t= "<<npart_f<<std::endl;
            std::cout<<"nu \t\t\t\t\t\t\t= "<<nu<<std::endl;
            std::cout<<"LL \t\t\t\t\t\t\t= "<<LL<<std::endl;
            std::cout<<"num_Lz_sectors\t\t\t\t\t\t= "<<num_Lz_sectors<<std::endl;

            std::cout<<"metropolis thermalization for P and Q matrices \t\t= "<<metrop_initial<<std::endl;
            std::cout<<"number of metroplis iteartions per PC for P and Q\t= "<<num_iter_metrop*iter<<std::endl;


            std::cout<<"metropolis thermalization for G.state \t\t\t= "<<initial_normal<<std::endl;
            std::cout<<"number of metroplis iteartions per PC for G.state\t= "<<iter_normal*num_L<<std::endl;

            std::cout<<"save plot to ./Plots/ \t\t\t\t\t= "<<plot<<std::endl;
     
        }
         
        //  MPI synchronize all parameters except the file path

         
        MPI_Bcast(&ne,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&npart_i,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&npart_f,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&nu,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&LL,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&num_Lz_sectors,1,MPI_INT,0,MPI_COMM_WORLD);

        MPI_Bcast(&iter,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&num_iter_metrop,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&metrop_initial,1,MPI_INT,0,MPI_COMM_WORLD);

        MPI_Bcast(&iter_normal,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&num_L,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&initial_normal,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&plot,1,MPI_CHAR,0,MPI_COMM_WORLD);

    }
};



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////








//#define DEBUG_NEW new(__FILE__, __LINE__)
//#define new DEBUG_NEW




using namespace std;
typedef complex<double> dcmplx;
const long double PI=3.1415926535897932384626433832795028841972;


MersenneTwister mt;
double ** binom;

extern "C"{
	//	Define external lapack routine to get the singular value decomposition
//	int zgesvd_(char*, char*, int*, int*, dcmplx*,int*, double*,dcmplx*, int*, dcmplx*, int*,dcmplx*, int*, double*, int*);
//        int zgetrf_(int*, int*, double*, int*, int*, int* );
        int zgeev_(char*, char*, int*, dcmplx* , int* , dcmplx* , dcmplx*, int*, dcmplx*, int*, dcmplx*, int*, double*, int* );
        int zgetrf_ ( int* M, int* N, dcmplx *A, int* LDA, int *IPIV, int *INFO );
          }


int main(int argc, char* argv[])
{
  int npart, i, j, k, l, m, dim, LL, ne, nu, lwork, ldu, info, iter, metrop_initial, Lz_sector,
  num_Lz_sectors, nu_eff, *dim_M_Lz_EWF, numtasks, rank, number_Lz_sectors, iter_normal, num_L,
  initial_normal, npart_i, npart_f;

  double Lz, nee, *Lz_EWF, *Lz_EWF_0, topbinom, Lz_0, rwork[5000], norma_rhoA=0.0, num_iter_metrop, QJain, N1_Jain, N2_Jain, 
  *** M_Lz_EWF, cut_sph, norm_EWF_A, norm_EWF_B, t1, t2, norm_Gstate, norm_Error, **Eigenvalues_Lz; 

  dcmplx *zu_A, *zv_A, *zu_B, *zv_B, state_A, state_B, state_A_comp, state_B_comp, *zz_normal, *zz, ***P, ***Q, *Mat, 
  *Eigen, *work, dummy[1][1], state_0, state_0_comp, *M_state, *P_send, *P_recv, *Q_send, *Q_recv, alpha, norma_state, norm_gstate,
   means_distrib;

  string path0, line, name; 
  char plot;

  vector< vector <double> > Fspace, Fspace_comp;
  vector<double> gnu_x;   /// for gnuplot
  vector<double> gnu_y;   /// for gnuplot


//  clock_t t1=clock(),t2;

  stringstream fileName, fileName_1;   
  ifstream f_in, f_in_1;
  ofstream f_out, f_out_1;

  Wave_Functions *wf;





///////////////////////////////////////////////////////////////////////////////////////////
//////////////////// INITIALIZATION OF MPI AND RANDOM NUMBERS /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

t1=MPI_Wtime();

rand_init(argv);

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////













////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// DEFINE PHYSICAL VARIABLES ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////



/*
ne=10;     //// number of electrons
npart_i=5; //// npart_i=NA_initial, npart_f=NA_final. SET one of them or both TO -1 IF YOU WANT ALL THE NA SECTORS
npart_f=5; 
LL=1;     //// Landau levels  LL=1,2,3,4,... 
nu=1;     //// see comment below

// Comment: I am using this notation :
           //// Laughlin state = prod_{i<j} (z_i - z_j)^nu    (LL=1, nu , i.e. for Laughlin nu is the inverse of the filling)
           //// Jain state = IQH state(with filling = LL) * prod_{i<j} (z_i - z_j)^nu  ( for Jain nu are the fluxes ) 

number_Lz_sectors=6; //// number of Lz_A sectors. SET TO -1 IF YOU WANT ALL THE LZA SECTORS

/////// Metropolis variables for P and Q matrices/////////

iter=40;                /// number of internal iterations in the Metropolis. Iter >~ 2*Corr_time
num_iter_metrop=40;  /// number of metropolis interations. Using MPI with let's say 'np' processors you will compute  np*num_iter_metrop 
                        /// metropolis iterations in the Montecarlo. 
metrop_initial=1500;    /// number of iterations for Metropolis thermalization. 

///////////////////////////////////////////////////////////////////////////


/////// Metropolis variables for G.state norm /////////

/// Comment: number of metropolis iteration = iter_normal*num_L*np , with np the number of MPI processors 


iter_normal=40; 
num_L=30; 
initial_normal=500;

///////////////////////////////////////////////////////
*/




////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// SET VARIABLES ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//system("clear");

Params params(argc,argv);



cin.get();

ne = params.ne;
npart_i = params.npart_i;
npart_f = params.npart_f;
nu = params.nu;
LL = params.LL;
number_Lz_sectors = params.num_Lz_sectors;
iter = params.iter;
num_iter_metrop = params.num_iter_metrop;
metrop_initial = params.metrop_initial;
iter_normal = params.iter_normal;
num_L = params.num_L;
initial_normal = params.initial_normal;
path0 = params.path;
plot = params.plot;



cut_sph=PI/2.0+0.0; //// this is the cut on the sphere. This part of the code is not finished yet!

char flux='P'; //// P if you want Jain in positive field and N for negative field. The negative case is not finished yet!

string norm_gs="Y" ; /// Y: if you want to compute norm of G.State (normalized rho_A). 
                   /// N: for unnormalized rho_A.
                   /// ext_file: to get normalization from external file from Norm_Gstate/Norm_Ne_flux_LL_Num_iter.txt



////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// BINOMIAL //////////////////////////////////////////////////

 topbinom=60;

 binom=new double* [int(topbinom+2)];
   for (i=0;i<topbinom+1;i++)
 binom[i]=new double[int(topbinom+2)]();

 binomial (binom, topbinom);

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////









///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// COMPUTE NORM G.STATE ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


/// Comment : I'm not computing the exact norm. I've introduced an scaling alpha (I cancell that scaling 
///            at the end of the calculations) otherwise the norm explode for N~40.   


if (norm_gs=="ext_file")
{

 fileName.str("");   
 fileName<<path0<<"/Norm_Gstate/Norm_N_"<<ne<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<iter_normal*num_L*numtasks<<".dat";
 name = fileName.str();

f_in.open (name.c_str(), ios::in);

if (f_in.is_open())
{
getline (f_in, line);
f_in>>norm_Gstate;
f_in>>alpha;

f_in.close();
}
else
{
 cout<<"file not found"<<endl;
 MPI_Finalize();
 return 0;
 }


}
else if (norm_gs=="Y")
      {

      if (ne<=60)                       
        {
         zz_normal=new dcmplx [ne]();   

         for (i=0;i<ne;i++)
             zz_normal[i]=dcmplx( nfmod(mt.random(),PI), nfmod(mt.random(),2.0*PI));

         metropolis *metrop1 = new metropolis(zz_normal, ne, iter, 0.0, ne, nu, LL, 1, PI/2.0);



//    metrop1->width=1.5;
//    cout<<metrop1 -> Wfunction_Coulomb_energy (300, 100)<<endl;
//    cin.get();




/////////////////////////// Autocorrelation and Thermalization ///////////////////////////////////
//        if (rank==0)                                                                          //
//           {                                                                                  // 
//           metrop1 -> Wfunction_Norm_Thermal (1.5, 5000, 1600);                               //
//           metrop1 -> Wfunction_Corr_IQHE_N(5000, 1600);                                      // 
//                                                                                              //
//           metrop1->width=1.5;                                                                //
//           cout<<metrop1 -> Wfunction_Coulomb_energy (300, 100)<<endl;                        //
//           cin.get();  									//
//           return 0;                                                        			//
//           }                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////
 
      norm_gstate = metrop1 -> Wfunction_Norm(0.0, iter_normal, initial_normal, iter_normal*num_L, alpha); //dcmplx(norm, error)

      norm_Gstate=real(norm_gstate);    //// this is not the real value but including scaling alpha //// 
      norm_Error=imag(norm_gstate);     //// this is not the real value but including scaling alpha ////
 

         metrop1 -> ~metropolis();
        
         delete [] zz_normal;



     if (rank==0)
        {
         fileName<<path0<<"/Norm_Gstate/Norm_N_"<<ne<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<iter_normal*num_L*numtasks<<".dat";
         fileName.str("");   
         fileName<<path0<<"/Norm_Gstate/Norm_N_"<<ne<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<iter_normal*num_L*numtasks<<".dat";
         name = fileName.str();
        
         f_out.open (name.c_str());

         f_out<<"Rescale values: norm, alpha, relative_error: "<<endl;
         f_out <<norm_Gstate<<endl;
         f_out <<alpha<<endl;
         f_out <<norm_Error/norm_Gstate<<endl;

         f_out<<"Real values: Norm, relative_error: "<<endl;
         f_out <<exp(alpha)*norm_Gstate<<endl;
         f_out <<norm_Error/norm_Gstate<<endl;

         f_out.close();
         }


        }

    else
        {  
        norm_Gstate=1.0;
        alpha=0.0;
        }

   }

else if (norm_gs=="N")
     {        
      norm_Gstate=1.0;
      alpha=0.0; 
      }





///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////








if (npart_i==-1 || npart_f==-1) 
   {
    npart_i=0;
    npart_f=ne;
    }


for (npart=npart_i; npart<=npart_f; ++npart)
    {




/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////// FIND L_0 ; THE MINIMUM VALUE OF LZ_A /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


nee=ne;
QJain=(nee/2.0-2.0)/2.0; N1_Jain=int((2.0*QJain+1.0)/2.0); N2_Jain=npart-N1_Jain;


if (LL==1)
nu_eff=0;
else 
nu_eff=1;



double LLd=LL;
double QJain_0th=(nee/LLd-LLd)/2.0;
double QJain_1st=(nee/LLd-LLd+2)/2.0;
double QJain_2nd=(nee/LLd-LLd+4)/2.0;
int N_0th=(2.0*QJain_0th+1.0)/2.0;
int N_1st=(2.0*QJain_1st+1.0)/2.0;
int N_2nd=npart-N_0th-N_1st;




if (LL==1)
   Lz_0=(npart*(npart-1)/2.0 - (ne-1)*npart/2.0)*nu;   //FOR LAUGHLIN filling 1/nu
else if (LL==2)
        {
         Lz_0 = N1_Jain*(N1_Jain-1.0)/2.0 - QJain*N1_Jain + N2_Jain*(N2_Jain-1.0)/2.0 - (QJain+1.0)*N2_Jain + 
         (npart*(npart-1)/2.0 - (ne-1)*npart/2.0)*nu;   // FOR JAIN nu=2/5 or nu=2/3
          }
else if (LL==3)
        {
         Lz_0 = N_0th*(N_0th-1.0)/2.0 - QJain_0th*N_0th + N_1st*(N_1st-1.0)/2.0 - QJain_1st*N_1st + N_2nd*(N_2nd-1.0)/2.0 - QJain_2nd*N_2nd +
         (npart*(npart-1)/2.0 - (ne-1)*npart/2.0)*nu;   // FOR JAIN nu=3/7 or nu=3/5
          }




 if (number_Lz_sectors==-1)
     num_Lz_sectors=2.0*abs(Lz_0)+1;
 else if(number_Lz_sectors>2.0*abs(Lz_0)+1)
        {
         //cout<<"ERROR - MAX SECTOR IS "<<2.0*abs(Lz_0)+1<<endl;         
         //MPI_Finalize();         
         //return -2;
         
         cout<<"WARNING - MAX SECTOR IS "<<2.0*abs(Lz_0)+1<<endl;
         num_Lz_sectors=2.0*abs(Lz_0)+1;
         
         }
    else 
        num_Lz_sectors=number_Lz_sectors;


////////////////////////////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////

















//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// LOG OF NORMALIZATION OF GROUND STATE AND METROPOLIS PROBABILITY   /////////////////////
//////////////////////////////         TO COMPUTE THE NORMALIZED RHO_A                    ////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



  if (ne<=60)
   {
    norm_EWF_A=0.0;
    for (i=0;i<npart;++i)
         norm_EWF_A+=log(4.0*PI*Beta_I ( pow(cos(cut_sph/2.0),2.0) , 1.0+i, (double)(ne-i)));  
    norm_EWF_A+=log(factorial((int)npart));

    norm_EWF_B=0.0;
    for (i=npart;i<ne;++i)
        norm_EWF_B+=log(4.0*PI*Gamma(1.0+i) *
        ( (Gamma((double)(ne-i))/Gamma(ne+1.0)) - (1.0+i)*Beta_I(pow(cos(cut_sph/2.0),2.0) , 1.0+i,(double)ne-i)/Gamma(2.0 + i) ));
    norm_EWF_B+=log(factorial((int)(ne -npart))); 
   }
else
    {
     norm_EWF_A=0.0;
     norm_EWF_B=0.0;
     }



////////////////////////////////////////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////   ENTANGLEMENT SPECTRUM CALCULATION        /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// CREATE SOME POINTERS ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

     zu_A=new dcmplx [ne]();
     zv_A=new dcmplx [ne]();
     zu_B=new dcmplx [ne]();
     zv_B=new dcmplx [ne]();
     zz=new dcmplx [ne]();   
     Lz_EWF=new double [3*ne+LL]();
     Lz_EWF_0=new double [3*ne+LL]();

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////









/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// GENERATE CONFIG OF EWF FOR P AND Q MATRIX ////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////// 



 dim_M_Lz_EWF=new int [num_Lz_sectors+1];
 M_Lz_EWF=new double** [num_Lz_sectors+1];



 
 for (Lz_sector=0; Lz_sector<num_Lz_sectors; ++Lz_sector)  
     {
      Lz = Lz_0 + Lz_sector;


      if (npart==0)
         {
          dim=1;
          Fspace.resize(dim);
          for (i=0;i<dim;++i)
          Fspace[i].resize(dim);
          Fspace[0][0]=0.0;

          }
       else if ((ne-npart)==0)
         {
          EWF_Space(Lz, ne, dim, ne, LL, Fspace, Fspace_comp, nu);
          dim=1;
          }
        else
            EWF_Space(Lz, ne, dim, npart, LL, Fspace, Fspace_comp, nu);


      dim_M_Lz_EWF[Lz_sector]=dim;          

         
       M_Lz_EWF[Lz_sector]=new double *[dim+1];
       for (j=0;j<dim;j++)
             M_Lz_EWF[Lz_sector][j]=new double [(nu+1)*ne+LL]();





    for (i=0;i<dim;++i)
        {

         for (k=0;k<LL;++k)          
              M_Lz_EWF[Lz_sector][i][k]=Fspace[i][k];
              

         if (LL==1)
            {
             for (k=0;k<nu*npart+1;++k)
                  M_Lz_EWF[Lz_sector][i][k+LL]=Fspace[i][k+LL];

             }
         else if (LL==2 || LL==3)
                 {
                  for (k=0;k<nu*npart+npart+1;++k)
                       M_Lz_EWF[Lz_sector][i][k+LL]=Fspace[i][k+LL];
                  }
         }

     }




///////////////// FIND MAXIMUM DIMENSION OF HILBERT SPACE AND CREATE M_STATE TO USE IN P AND Q MATRICES ///////////////////////


    int *dim_M_Lz_EWF_aux;
    dim_M_Lz_EWF_aux=new int [num_Lz_sectors+2]();
    
    for (i=0;i<num_Lz_sectors;++i)
        dim_M_Lz_EWF_aux[i]=dim_M_Lz_EWF[i];


    int dim_aux;
    for (i=0;i<num_Lz_sectors-1;++i)
        for (j=i+1;j<num_Lz_sectors;++j)
            {
             if (dim_M_Lz_EWF_aux[i]<dim_M_Lz_EWF_aux[j])
                {
                 dim_aux=dim_M_Lz_EWF_aux[i];
                 dim_M_Lz_EWF_aux[i]=dim_M_Lz_EWF_aux[j];
                 dim_M_Lz_EWF_aux[j]=dim_aux;
                 }
             }


    M_state=new dcmplx[dim_M_Lz_EWF_aux[0]];
 
    delete [] dim_M_Lz_EWF_aux;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// P MATRIX CALCULATION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// CREATE P MATRIX ARRAY AND CONSTRUCTOR FOR W.FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////







        P=new dcmplx** [num_Lz_sectors];
        for (i=0;i<num_Lz_sectors;i++)
             P[i]=new dcmplx *[dim_M_Lz_EWF[i]];
        for (i=0;i<num_Lz_sectors;i++)
            for (j=0;j<dim_M_Lz_EWF[i];j++)
                P[i][j]=new dcmplx [dim_M_Lz_EWF[i]];





//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// INITIALIZE COORDINATES FOR METROPOLIS ALGORITHM /////////////////////

     for (i=0;i<npart;i++)
          zz[i]=dcmplx(((PI-cut_sph-0.02)*i/npart) + cut_sph + 0.01, nfmod(mt.random(),2.0*PI));
//          zz[i]=dcmplx((PI*(i+1)/ne)/2.0 + PI/2.0 - 0.01, nfmod(mt.random(),2.0*PI));

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CONSTRUCTOR FOR CLASS METROPOLIS AND W.FUNCTION//////////////////////

     wf = new Wave_Functions(ne, npart, nu, LL, flux);


     metropolis *metropA=new metropolis(zz, ne, iter, 0.0, npart, nu, LL, nu_eff, cut_sph);

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



      





/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// METROPOLIS FOR P MATRIX ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////// 




if (npart!=0)
   {


//////////////////////////////////////// INITIALIZATION /////////////////////////////////////////

metropA -> width=0.0;

     for (i=1; i<=metrop_initial; i++)
         {
          metropA -> run_A();



     if (rank==0 && i%400==0)
         cout<<"P initialization -> accepted,  rejected  :"<<metropA -> accepted<<" "<<metropA -> rejected<<endl;

 
          metropA -> equilibration (metropA->accepted, metropA->rejected);
         
          }

///////////////////////////////////////////////////////////////////////////////////////////////// 

	
////////////////////////////// MONTECARLO INTEGRATION FOR P MATRIX //////////////////////////////

      for (m=0;m<num_iter_metrop*iter;++m)
          {


           metropA -> run_A();
           if (rank==0 && m%500==0)
           cout<<"P Metrop. iteration: "<<m<<" of  "<<num_iter_metrop*iter<<" for LzA: "<<Lz_sector<<endl;
 
           for (i=0;i<npart;++i)
	       {
	        zu_A[i]= cos(real(metropA -> zz[i]/2.0))*exp(dcmplx(0,imag(metropA -> zz[i])/2.0));
	        zv_A[i]= sin(real(metropA -> zz[i]/2.0))*exp(dcmplx(0,-imag(metropA -> zz[i])/2.0));
      	        }         



          if (flux=='N')
             state_0=metropA -> Jain_negative_2_LL_A (zu_A, zv_A);
          else
             state_0=metropA -> IQHE_Prob_A(zu_A, zv_A);             
             
 
                

          for (l=0; l<num_Lz_sectors; ++l)
               {
      
               for (i=0;i<dim_M_Lz_EWF[l];++i)
                   {       
                    wf-> w_func (state_A, M_Lz_EWF[l][i], ne, zu_A, zv_A, binom, nu, LL);
                    M_state[i]=state_A; 

                    }

                            
           
               for (i=0;i<dim_M_Lz_EWF[l];++i)
                   for (j=i;j<dim_M_Lz_EWF[l];++j)           
                       P[l][i][j]+=exp(norm_EWF_A + M_state[i]-state_0 + conj(M_state[j])-conj(state_0)  -
                                    alpha/(2.0*1.0) )/(1.0*num_iter_metrop*iter*norm_Gstate);                         

                }       
             

           }

////////////////////////////////////////////////////////////////////////////////////////////////////////


    } 
    
      else    //// if (npart==0)
      P[0][0][0]=exp(-alpha/2.0)*1.0/(1.0*norm_Gstate);






///////////////////////// DESTRUCTOR FOR W.FUNCTION AND DELETE P MATRIX TO SAVE SPACE //////////////////


wf -> ~Wave_Functions();    
metropA -> ~metropolis();


     
         for (l=0;l<num_Lz_sectors; ++l)
             for (i=0;i<dim_M_Lz_EWF[l];++i)
                 for (j=0;j<i;++j)
                     P[l][i][j]=conj(P[l][j][i]);



          

/////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////





















/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// Q MATRIX CALCULATION //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// GENERATE CONFIG OF EWF FOR Q MATRIX //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////// 




 for (Lz_sector=0; Lz_sector<num_Lz_sectors; ++Lz_sector)
     {
      Lz=Lz_0 + Lz_sector;

      if (npart==0) 
         {
          EWF_Space(0, ne, dim, ne, LL, Fspace, Fspace_comp, nu);
          Fspace_comp=Fspace;
          dim=1;
          }
      else if ((ne-npart)==0)
          {
          dim=1;
          Fspace_comp.resize(dim);
          for (i=0;i<dim;++i)
          Fspace_comp[i].resize(dim);
          Fspace_comp[0][0]=0.0;
           }
      else
         EWF_Space(Lz, ne, dim, npart, LL, Fspace, Fspace_comp, nu);


     


     
     for (i=0;i<dim;++i)
         {
         for (k=0;k<LL;++k)          
              M_Lz_EWF[Lz_sector][i][k]=Fspace_comp[i][k];

         if (LL==1)
            {
             for (k=0;k<nu*(ne-npart)+1;++k)
                  M_Lz_EWF[Lz_sector][i][k+LL]=Fspace_comp[i][k+LL];
             }  
         else if (LL==2 || LL==3)
                 {
                  for (k=0;k<nu*(ne-npart)+ne-npart+1;++k)
                       M_Lz_EWF[Lz_sector][i][k+LL]=Fspace_comp[i][k+LL];
                  }  
         M_Lz_EWF[Lz_sector][i][(int)(nu*(ne-npart)+nu_eff*(ne-npart)+LL)]=1.0;

         }  
     }



///////////////////////////////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////// 




//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// CREATE Q MATRIX ARRAY AND CONSTRUCTOR FOR W.FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////






        Q=new dcmplx** [num_Lz_sectors];
        for (i=0;i<num_Lz_sectors;i++)
             Q[i]=new dcmplx *[dim_M_Lz_EWF[i]];
        for (i=0;i<num_Lz_sectors;i++)
            for (j=0;j<dim_M_Lz_EWF[i];j++)
                Q[i][j]=new dcmplx [dim_M_Lz_EWF[i]];




     for (i=0;i<ne-npart;i++)
          zz[i]=dcmplx(((cut_sph-0.02)*i/(ne-npart)) + 0.01, nfmod(mt.random(),2.0*PI));
//          zz[i]=dcmplx((PI*(i+1)/ne) -0.01, nfmod(mt.random(),2.0*PI));

     

     wf = new Wave_Functions(ne, ne-npart, nu, LL, flux);
     metropolis *metropB = new metropolis(zz, ne, iter, 0.0, npart, nu, LL, nu_eff, cut_sph);




//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////







/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// METROPOLIS FOR Q MATRIX ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////// 

  

//////////////////////////////////////// INITIALIZATION /////////////////////////////////////////


if ((ne-npart)!=0)
   {

metropB->width=0.0;

     for (i=1; i<=metrop_initial; i++)
         {
          metropB -> run_B();

     if (rank==0 && i%400==0)
         cout<<"Q initialization -> accepted,  rejected  :"<<metropB -> accepted<<" "<<metropB -> rejected<<endl;

          metropB -> equilibration (metropB->accepted, metropB->rejected);

          }





//////////////////////////////////////////////////////////////////////////////////////////////////////


   

////////////////////////////// MONTECARLO INTEGRATION FOR Q MATRIX //////////////////////////////



      for (m=0;m<num_iter_metrop*iter;++m)
          {

            metropB -> run_B();
           if (rank==0 && m%500==0)
           cout<<"Q Metrop. iteration: "<<m<<" of  "<<num_iter_metrop*iter<<" for LzA: "<<Lz_sector<<endl;
 
           for (i=0;i<ne-npart;++i)
	       {
	        zu_B[i]= cos(real(metropB -> zz[i]/2.0))*exp(dcmplx(0,imag(metropB -> zz[i])/2.0));
	        zv_B[i]= sin(real(metropB -> zz[i]/2.0))*exp(dcmplx(0,-imag(metropB -> zz[i])/2.0));
      	        }


           if (flux=='N')
           state_0_comp=metropB -> Jain_negative_2_LL_B (zu_B, zv_B);
           else
           state_0_comp=metropB -> IQHE_Prob_B(zu_B, zv_B);
                 

  
           for (l=0; l<num_Lz_sectors; ++l)
               {
                          
               for (i=0;i<dim_M_Lz_EWF[l];++i)
                   {
                    wf-> w_func (state_A_comp, M_Lz_EWF[l][i], ne, zu_B, zv_B, binom, nu, LL);
                    M_state[i]=state_A_comp; 
                    }

           
               for (i=0;i<dim_M_Lz_EWF[l];++i)
                   for (j=i;j<dim_M_Lz_EWF[l];++j)
                        Q[l][i][j]+=exp(norm_EWF_B + M_state[i]-state_0_comp + conj(M_state[j])-conj(state_0_comp)  - 
                                    alpha/(2.0*1.0))/((double)num_iter_metrop*iter);
                

                }
       
           }

   }
    
    else  /// if ((ne-npart)==0)    
        Q[0][0][0]=exp(-alpha/2.0);

//////////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////// DESTRUCTOR FOR W.FUNCTION AND DELETE Q MATRIX TO SAVE SPACE //////////////////


wf -> ~Wave_Functions();    
metropB -> ~metropolis();

     
         for (l=0;l<num_Lz_sectors; ++l)
             for (i=0;i<dim_M_Lz_EWF[l];++i)
                 for (j=0;j<i;++j)
                     Q[l][i][j]=conj(Q[l][j][i]);


/////////////////////////////////////////////////////////////////////////////////////////////////////////






/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////








////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// COMPUTE EIGENVALUES OF M=P*Q MATRIX /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 


//          
 

 
//            fileName<<path<<"/Eigenvalues_n_"<<ne<<"_na_"<<npart_i<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<num_MPI_iter<<".dat";
 
//            name = fileName.str();


     
//       char path1[100]; 
       int num_MPI_iter=numtasks*num_iter_metrop*iter;

   if (rank==0)
      {
       fileName.str("");   
       fileName<<path0<<"/Eigenvalues_N_"<<ne<<"_NA_"<<npart<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<num_MPI_iter<<".dat";           
       name = fileName.str();

       f_out.open (name.c_str());
       f_out.close();


        Eigenvalues_Lz = new double* [num_Lz_sectors];
        for (i=0;i<num_Lz_sectors; ++i)
            Eigenvalues_Lz[i]=new double[dim_M_Lz_EWF[i]];            

      }






   
   for (l=0 ; l<num_Lz_sectors; ++l)
       {

       dim=dim_M_Lz_EWF[l];


       P_send = new dcmplx [dim*dim];    
       Q_send = new dcmplx [dim*dim];

       if (rank==0)
          {
           Mat = new dcmplx [dim*dim];
           P_recv = new dcmplx [dim*dim];
           Q_recv = new dcmplx [dim*dim];  
           }
  



       for (i=0; i<dim; ++i)
           for (j=0; j<dim; ++j)
               P_send[j+dim*i]=P[l][j][i];

       for (i=0; i<dim; ++i)
           for (j=0; j<dim; ++j)
               Q_send[j+dim*i]=Q[l][j][i];



      MPI_Reduce (P_send,P_recv,dim*dim,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce (Q_send,Q_recv,dim*dim,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);


      delete [] P_send; delete [] Q_send;






 ////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////// DIAGONALIZATION IN MASTER NODE - RANK 0 ////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////





      if (rank==0)
        {


         for (i=0; i<dim; ++i)
             for (j=0; j<dim; ++j)
                 for (k=0;k<dim;++k)
                      {
                       if (ne<=60) 
                          Mat[j+dim*i]+=P_recv[j+dim*k]*Q_recv[k+dim*i]*binom[ne][npart]*pow(1.0*numtasks,-2);   
                       else
                           Mat[j+dim*i]+=P_recv[j+dim*k]*Q_recv[k+dim*i]*pow(1.0*numtasks,-2); 
                       }
            // the binom comes from normalization of rho_A and the pow(numtasks, -2) from the MPI_SUM in MPI_reduce //
     

       delete [] P_recv; delete [] Q_recv;






        Eigen = new dcmplx [dim];
     
    


        lwork=10000;
        work = (dcmplx*)malloc( lwork*sizeof(dcmplx) );
        ldu=1;


                       
        zgeev_((char*)"N",(char*)"N",&dim,Mat,&dim,Eigen,(dcmplx*)dummy,&ldu,(dcmplx*)dummy,&ldu,work,&lwork,rwork,&info);





        for (i=0; i<dim; ++i)
            {          



             if (abs(Eigen[i])<pow(10.0,-10.0))
                Eigen[i]=0.0;

            if (abs(Eigen[i])!=0.0)
            Eigenvalues_Lz[l][i]=-log((double)(real(Eigen[i])));
            else 
            Eigenvalues_Lz[l][i]=1000.0;

             }





       f_out.open (name.c_str(), ios::app);

       for (i=0;i<dim;++i)
            {
             f_out<< l <<" "<< Eigenvalues_Lz[l][i]<<endl;
             }
 
       f_out.close();



      if (npart==-1)
      for (i=0;i<dim;++i)
           norma_rhoA+=real(Eigen[i]);


      delete [] Mat; delete Eigen;
      
     
      }

    


  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// PLOT ENTANGLEMENT SPECTRUM TO EXTERNAL FILE //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if (rank==0 && plot=='y')
     {
 
     cout << "*** start of GNUPLOT ***" << endl;
 
     try {
          Gnuplot g1 = Gnuplot("dots");
    
     for (i = 0; i < num_Lz_sectors; ++i)
         for (j=0; j<dim_M_Lz_EWF[i]; ++j)       
       {
         gnu_x.push_back((double)i);
         gnu_y.push_back(Eigenvalues_Lz[i][j]);
       }


      
     fileName.str("");   
     fileName<<"set xrange [-0.5:"<<num_Lz_sectors+0.5<<"]\n";
     fileName<<"set yrange [-0.5:"<< *min_element(Eigenvalues_Lz[num_Lz_sectors-1],Eigenvalues_Lz[num_Lz_sectors-1]+ dim_M_Lz_EWF[num_Lz_sectors-1]) + 4  <<"]\n";



     name = fileName.str();
      
     g1.cmd(name.c_str());


      
     g1.cmd("set ylabel 'ES energies'\n");
     g1.cmd("set xlabel 'LzA sectors'\n");



     fileName.str("");   
     fileName<<"ES with Metropolis Method: N="<<ne<<", NA="<<npart<<", flux(nu)="<<nu<<", LL="<<LL<<", iter="<<num_MPI_iter;
     name = fileName.str();


     g1.cmd("set pointsize 1.2");
     g1.set_style("points");



     g1.plot_xy(gnu_x,gnu_y,"ES with Metropolis Method: ");



    fileName.str("");   
    fileName<<"set output "<<" 'Plots/Plot_N_"<<ne<<"_NA_"<<npart<<"_flux_"<<nu<<"_LL_"<<LL<<"_iter_"<<num_MPI_iter<<".png'\n";
    name = fileName.str();

    cout<<name<<endl;

    g1.cmd("set terminal pngcairo\n");
    g1.cmd(name.c_str());
    g1.cmd("replot\n");



    sleep(3);



//     cin.get();
     g1.reset_plot();


     } catch (GnuplotException ge) {
         cout << ge.what() << endl;
     }

     cout << endl << "*** end of GNUPLOT ***** " << endl;


    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// FREE POINTERS MEMORY ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



   for (i=0;i<num_Lz_sectors;i++)
       for (j=0;j<dim_M_Lz_EWF[i];j++)
            delete []  Q[i][j];
   for (i=0;i<num_Lz_sectors;i++)
       delete [] Q[i];
   delete [] Q;


   for (i=0;i<num_Lz_sectors;i++)
       for (j=0;j<dim_M_Lz_EWF[i];j++)
            delete []  P[i][j];
   for (i=0;i<num_Lz_sectors;i++)
       delete [] P[i];
   delete [] P;



   delete [] M_state;

 for (i=0;i<num_Lz_sectors;i++)
       for (j=0;j<dim_M_Lz_EWF[i];j++)
            delete []  M_Lz_EWF[i][j];
   for (i=0;i<num_Lz_sectors;i++)
       delete [] M_Lz_EWF[i];
   delete [] M_Lz_EWF;



 delete [] Lz_EWF_0; delete [] Lz_EWF; delete [] zu_A; delete [] zv_A; delete [] zu_B; delete [] zv_B; delete [] zz; delete [] dim_M_Lz_EWF;




   if (rank==0)
      { 
       for (i=0;i<num_Lz_sectors;i++)
           delete [] Eigenvalues_Lz[i];
       delete [] Eigenvalues_Lz;
       }




MPI_Barrier(MPI_COMM_WORLD);


}  /// end for npart ///



if (rank==0 || npart==-1)
   cout<<"norm_rho_A:  "<<norma_rhoA<<endl;



   for (i=0;i<topbinom+1;i++)
        delete [] binom[i];
   delete [] binom;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



MPI_Barrier(MPI_COMM_WORLD);


t2 = MPI_Wtime();

if (rank==0)
cout<<"time consume "<<t2-t1<<endl;

MPI_Finalize();


return 0;

///////////////////////////////////// END OF CODE ////////////////////////////////////////////////////


}
































