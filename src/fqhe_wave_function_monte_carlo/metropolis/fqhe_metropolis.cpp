////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 06/12/2013
//!
//!  \file 
//!		These functions perform metropolis moves, evaluate running average
//!		energy and write data files.
//!
//!                    Copyright (C) 2013 Simon C Davenport
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

#include "fqhe_metropolis.hpp"

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the Metropolis class
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::MonteCarloData::InitFromCommandLine(
    #if _ENABLE_MPI_
    boost::program_options::variables_map* options, //!<    Command line argument list
    const utilities::MpiWrapper& mpi)               //!<    Mpi wrapper class
    #else
    boost::program_options::variables_map* options)  //!<    Command line argument list
    #endif
{
    #if _ENABLE_MPI_
    if(0 == mpi.m_id)    //  For the master node
    #endif
    {
        nbrTherm     = (*options)["nbr-therm"].as<unsigned int>();
        nbrSamples   = (*options)["nbr-samples"].as<unsigned int>();
        runLabel     = (*options)["run"].as<unsigned int>();
        maxDist      = (*options)["maxd"].as<double>();
        sampleFreq   = (*options)["sample-freq"].as<unsigned int>();
        nbrToMove    = (*options)["nbr-move"].as<unsigned int>(); 
        wave         = (*options)["wave"].as<bool>();
        useInitFile  = (*options)["use-init-file"].as<bool>();  
        initFileName = (*options)["init-file"].as<std::string>();
        path         = (*options)["path"].as<std::string>();
    }
    
    #if _ENABLE_MPI_
    this->MpiSync(0,mpi);
    #endif
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

#if _ENABLE_MPI_

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Synchronize all data values with those on the master node
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::MonteCarloData::MpiSync(
    int syncId,                         //!<    Label of node to sync with
    const utilities::MpiWrapper& mpi)   //!<    Mpi wrapper class
{
    mpi.Sync(&nbrTherm,1,syncId);
    mpi.Sync(&nbrSamples,1,syncId);
    mpi.Sync(&runLabel,1,syncId);
    mpi.Sync(&maxDist,1,syncId);
    mpi.Sync(&sampleFreq,1,syncId);
    mpi.Sync(&wave,1,syncId);
    mpi.Sync(&useInitFile,1,syncId);
    mpi.Sync(initFileName,syncId);
    mpi.Sync(path,syncId);
}

#endif

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Constructor for the Metropolis class
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Metropolis::Metropolis(
	FQHE::WaveFunctionData *wavefunctionData,		//!<	Address of data structure containing basic wave function data
	FQHE::MonteCarloData *montecarloData)			//!<	Address of data structure containing basic Monte Carlo data
	: 
	theta(0),
	phi(0),
	r(0),
	thetaNew(0),
	phiNew(0),
	wfData(wavefunctionData),
	mcData(montecarloData),
	u(0),
	v(0),
	uNew(0),
	vNew(0),
	uTmp(0),
	vTmp(0),
	z(0),
	zNew(0),
	zTmp(0),
	toMove(0)
{
	mcData->nbrThermSteps=wfData->nbr*mcData->nbrTherm;
	mcData->nbrSteps=wfData->nbr*mcData->nbrSamples*mcData->sampleFreq;

	if(mcData->maxDist>1 || mcData->maxDist<=0)
	{
		std::cerr<<"\n\tERROR! Proportion of max move distance must 0<maxDist<=1"<<std::endl;
		exit(EXIT_FAILURE);;
	}
	
	//	Seed the standard random number generator
				
	int seed;
	
	if(getenv("PBS_JOBID")==NULL)	seed=time(0);
	else							seed=atoi(getenv("PBS_JOBID"));
	
	//  For testing
    //seed=100;
    
	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_
		seed=1;
		//	keep the same seed each time in order to repeat results
	#endif
	////////////////////////////////////////////////////////////////////////////////
	
	srand(seed);
	
	//	Now seed the Mersenne Twister random number generator using
	//	an array of these standard random numbers
	
	int initArraySize=100;
	long unsigned int initArray[initArraySize];
	
	for(int i=0;i<initArraySize;i++)
	{
		initArray[i]=rand();
	}
	
	mt.init_by_array(initArray,initArraySize);
	
	//	Initialise variables which depend on n (sphere and disc geometry)

	switch(wfData->geometry)
	{
		case _SPHERE_:
			utilities::cout.MainOutput()<<"\tInitializing Metropolis algorithm...\n\n\tUsing sphere geometry"<<std::endl<<std::endl;
			
			//	Initialise arrays for multi-particle MC moves
			
			uTmp = new (std::nothrow) dcmplx[mcData->nbrToMove];
			vTmp = new (std::nothrow)dcmplx[mcData->nbrToMove];
			uNew = new (std::nothrow) dcmplx[mcData->nbrToMove];
			vNew = new (std::nothrow) dcmplx[mcData->nbrToMove];
			thetaNew = new (std::nothrow) double[mcData->nbrToMove];
			phiNew = new (std::nothrow) double[mcData->nbrToMove];
			toMove =  new (std::nothrow) int[mcData->nbrToMove];
			
			phi = new (std::nothrow) double[wfData->nbr];
			theta = new (std::nothrow) double[wfData->nbr];
			u = new (std::nothrow) dcmplx[wfData->nbr];
			v = new (std::nothrow) dcmplx[wfData->nbr];
			
			//	Also, we now set the sphere radius based on the fillNumerator, fillDenominator 
			//	and shift defined in an individual wave function class
			
			//wfData->discRadius=sqrt((double)(wfData->flux/2.0));//*(wfData->nbr/(wfData->nbr-1.0)));
			
			//	Read initial configuration from file if selected
			
			if(mcData->useInitFile)
			{
				utilities::cout.MainOutput()<<"\tInitializing Metropolis algorithm...\n\n\tUsing initial configuration from file..."<<std::endl;
				thetaFile=mcData->initFileName;
				thetaFile.append("_theta.dat");
				f_init.open(thetaFile.c_str(), std::ios::binary);
				
				if(f_init.is_open()==0)
				{
					std::cerr<<"\tERROR! cannot open file: "<<thetaFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tReading initial theta values from file:\t"<<thetaFile<<std::endl;
				}
				
				f_init.read(reinterpret_cast<char*>(theta),wfData->nbr*SIZE_OF_DOUBLE);
				
				f_init.close();
				
				phiFile=mcData->initFileName;
				phiFile.append("_phi.dat");
				f_init.open(phiFile.c_str(), std::ios::binary);
				
				if(f_init.is_open()==0)
				{
					std::cerr<<"\tERROR! cannot open file: "<<phiFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					std::cerr<<"\tReading initial phi values from file:\t"<<phiFile<<std::endl;
				}
				
				f_init.read(reinterpret_cast<char*>(phi),wfData->nbr*SIZE_OF_DOUBLE);
				
				f_init.close();
				mcData->nbrThermSteps=0;
			}
			else
			{
				//	Otherwise, generate a random uniform distribution of points on the surface on a sphere.
				for (int i=0;i<wfData->nbr;i++)
				{
					phi[i]=2*PI*(mt.random())-PI;
					theta[i]=acos(mt.random()*pow(-1.0,(int)mt.genrand_int31()%2));
				}
			}
		
			for (int i=0;i<wfData->nbr;i++)
			{
				u[i]=dcmplx(sqrt((1.0+cos(theta[i]))/2)*cos(phi[i]/2),+sqrt((1.0+cos(theta[i]))/2)*sin(phi[i]/2));
				v[i]=dcmplx(sqrt((1.0-cos(theta[i]))/2)*cos(phi[i]/2),-sqrt((1.0-cos(theta[i]))/2)*sin(phi[i]/2));
			}
			
			if(mcData->useInitFile)	PrintConfigurationSphere();
			
			//	Open data files
			
			fileName.str("");
			if(mcData->path!="")
			{	
				fileName<<mcData->path<<"/";
			}
			fileName<<wfData->wfFileName<<"_run_"<<mcData->runLabel;
			
			phiFile=fileName.str();
			phiFile.append("_phi.dat");
			f_phi.open(phiFile.c_str(), std::ios::binary);
			
			if(f_phi.is_open()==0)
			{
				std::cerr<<"\tERROR! cannot open file: "<<phiFile<<std::endl<<std::endl;
				exit(EXIT_FAILURE); 
			}
			else
			{
				utilities::cout.MainOutput()<<"\tWriting phi values to file:\t\t"<<phiFile<<std::endl;
			}
			
			thetaFile=fileName.str();
			thetaFile.append("_theta.dat");
			f_theta.open(thetaFile.c_str(), std::ios::binary);
			
			if(f_theta.is_open()==0)
			{
				std::cerr<<"\tERROR! cannot open file: "<<thetaFile<<std::endl<<std::endl;
				exit(EXIT_FAILURE); 
			}
			else
			{
				utilities::cout.MainOutput()<<"\tWriting theta values to file:\t\t"<<thetaFile<<std::endl;
			}
			
			if(mcData->wave)
			{
				waveFile=fileName.str();
				waveFile.append("_wave.dat");
				f_wave.open(waveFile.c_str(), std::ios::binary);
				
				if(f_wave.is_open()==0)
				{
					std::cerr<<"\tERROR! cannot open file: "<<waveFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tWriting wave function values to file:\t"<<waveFile<<std::endl;
				}
			}
			
			//	Display a summary of the simulation information
			
			utilities::cout.MainOutput()<<"\n\tSphere Radius\t\t\t"<<wfData->discRadius<<"l_0";
			
			break;
		case _DISC_:
			utilities::cout.MainOutput()<<"\tUsing disc geometry"<<std::endl<<std::endl;
			
			//	Initialise arrays for multi-particle MC moves
			
			zTmp = new (std::nothrow)dcmplx[mcData->nbrToMove];
			zNew = new (std::nothrow) dcmplx[mcData->nbrToMove];
			toMove = new (std::nothrow) int[mcData->nbrToMove];
			
			r = new (std::nothrow) double[wfData->nbr];
			theta = new (std::nothrow) double[wfData->nbr];
			z = new (std::nothrow) dcmplx[wfData->nbr];
			
			//	Also, we now set the disc radius based on the fillNumerator, fillDenominator 
			//	and shift defined in an individual wavefunction class
				
			// disc radius is only a rough estimate
			//wfData->discRadius=sqrt(2.0*wfData->flux);		
			
			//	Read initial configuration from file if selected
			
			if(mcData->useInitFile)
			{
				utilities::cout.MainOutput()<<"\tUsing initial configuration from file..."<<std::endl;
				thetaFile=mcData->initFileName;
				thetaFile.append("_theta.dat");
				f_init.open(thetaFile.c_str(), std::ios::binary);
				
				if(f_init.is_open()==0)
				{
					std::cerr<<"\n\tERROR! cannot open file: "<<thetaFile<<std::endl<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\n\tReading initial theta values from file:\t"<<thetaFile<<std::endl<<std::endl;
				}
				
				for (int i=0;i<wfData->nbr;i++)
				{
					f_init.read(reinterpret_cast<char*>(theta),wfData->nbr*SIZE_OF_DOUBLE);
				}
				
				f_init.close();
				
				rFile=mcData->initFileName;
				rFile.append("_r.dat");
				f_init.open(rFile.c_str(), std::ios::binary);
				
				if(f_init.is_open()==0)
				{
					std::cerr<<"\n\tERROR! cannot open file: "<<rFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tReading initial r values from file:\t"<<rFile<<std::endl;
				}
				
				for (int i=0;i<wfData->nbr;i++)
				{
					f_init.read(reinterpret_cast<char*>(r),wfData->nbr*SIZE_OF_DOUBLE);
				}
				
				f_init.close();
				mcData->nbrThermSteps=0;
			}
			else
			{
				//	Otherwise, generate a random uniform distribution of points within a disc
				for (int i=0;i<wfData->nbr;i++)
				{
					*(r+i)=wfData->discRadius*sqrt(mt.random());
					*(theta+i)=2*PI*(mt.random())-PI;
				}
			}		
			
			for (int i=0;i<wfData->nbr;i++)
			{
				*(z+i)=dcmplx(*(r+i)*cos(*(theta+i)),*(r+i)*sin(*(theta+i)));
			}
			
			if(mcData->useInitFile)	PrintConfigurationSphere();
			
			//	Open data files
			
			fileName.str("");
			if(mcData->path!="")
			{	
				fileName<<mcData->path<<"/";
			}

			fileName<<wfData->wfFileName<<"_run_"<<mcData->runLabel;
			
			rFile=fileName.str();
			rFile.append("_r.dat");
			f_r.open(rFile.c_str(), std::ios::binary);
			
			if(f_r.is_open()==0)
			{
				std::cerr<<"\n\tERROR! cannot open file: "<<rFile<<std::endl;
				exit(EXIT_FAILURE); 
			}
			else
			{
				utilities::cout.MainOutput()<<"\tWriting r values to file:\t\t"<<rFile<<std::endl;
			}
			
			thetaFile=fileName.str();
			thetaFile.append("_theta.dat");
			f_theta.open(thetaFile.c_str(), std::ios::binary);
			
			if(f_theta.is_open()==0)
			{
				std::cerr<<"\n\tERROR! cannot open file: "<<thetaFile<<std::endl;
				exit(EXIT_FAILURE); 
			}
			else
			{
				utilities::cout.MainOutput()<<"\tWriting theta values to file:\t\t"<<thetaFile<<std::endl;
			}
			
			if(mcData->wave)
			{
				waveFile=fileName.str();
				waveFile.append("_wave.dat");
				f_wave.open(waveFile.c_str(), std::ios::binary);
				
				if(f_wave.is_open()==0)
				{
					std::cerr<<"\n\tERROR! cannot open file: "<<waveFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tWriting wave function values to file:\t"<<waveFile<<std::endl;
				}
			}
			
			//	Read in file containing a discretized version of the electron-background 
			//	function bgFunc(r/radius) given in PRB 67 075304.
			
			f_elbg.open("../background_energy_data/electron_disc_background_function.dat");
			if (f_elbg.is_open()==0)
			{
				utilities::cout.MainOutput()<<"\tERROR! cannot open file:"<<"background_energy_data/electron_disc_background_function.dat"<<std::endl;
				exit(EXIT_FAILURE);;
			}
			
			bgFunc = new double[BG_FUNC_DIM];
			
			for(int j=0;j<BG_FUNC_DIM;j++)
			{
				f_elbg>>*(bgFunc+j);
			}
			
			f_elbg.close();
			
			//	Display a summary of the simulation information
			
			utilities::cout.MainOutput()<<"\n\tDisc Radius\t\t\t"<<wfData->discRadius<<"l_0";
			break;	
		case _TORUS_:
			utilities::cout.MainOutput()<<"\tUsing torus geometry"<<std::endl<<std::endl;
						
			//	Initialise arrays for multi-particle MC moves
			
			mcData->maxTorusDist=mcData->maxDist*wfData->Lx;
			
			latticeConfigTmp = new (std::nothrow)std::complex<int>[mcData->nbrToMove];
			latticeConfigNew = new (std::nothrow)std::complex<int>[mcData->nbrToMove];
			toMove = new (std::nothrow) int[mcData->nbrToMove];
			
			latticeConfig = new (std::nothrow) std::complex<int>[wfData->nbr];
	
			//	Read initial configuration from file if selected
						
			if(mcData->useInitFile)
			{
				utilities::cout.MainOutput()<<"\tInitializing Metropolis algorithm...\n\n\tUsing initial configuration from file..."<<std::endl;
				latticeConfigFile=mcData->initFileName;
				latticeConfigFile.append(".dat");
				f_init.open(latticeConfigFile.c_str(), std::ios::binary);
				
				if(f_init.is_open()==0)
				{
					std::cerr<<"\tERROR! cannot open file: "<<latticeConfigFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tReading initial theta values from file:\t"<<latticeConfigFile<<std::endl;
				}
				
				double realVal,imagVal;
				
				for(int i=0;i<wfData->nbr;i++)
				{	
					f_init.read(reinterpret_cast<char*>(&realVal),SIZE_OF_DOUBLE);
					f_init.read(reinterpret_cast<char*>(&imagVal),SIZE_OF_DOUBLE);
					
					*(latticeConfig+i)=dcmplx(realVal,imagVal);
				}
				
				f_init.close();
				
				mcData->nbrThermSteps=0;
			}
			else
			{
				//	Otherwise, generate a random uniform distribution of points on a torus
				//	Also impose the maxOccupancy condition
				
				bool pass=false;
				
				while(!pass)
				{
					for (int i=0;i<wfData->nbr;i++)
					{
						*(latticeConfig+i)=std::complex<int>(mt.genrand_int31()%wfData->Lx,mt.genrand_int31()%wfData->Ly);
					}
					
					pass=CheckOccupancy(latticeConfig);
				}
			}
			
			if(mcData->useInitFile)	PrintConfigurationTorus();
						
			//	Open data files
			
			fileName.str("");

			if(mcData->path!="")
			{	
				fileName<<mcData->path<<"/";
			}

			fileName<<wfData->wfFileName<<"_run_"<<mcData->runLabel;
			
			latticeConfigFile=fileName.str();
			latticeConfigFile.append(".dat");
			f_latticeConfig.open(latticeConfigFile.c_str(), std::ios::binary);
			
			if(f_latticeConfig.is_open()==0)
			{
				std::cerr<<"\tERROR! cannot open file: "<<latticeConfigFile<<std::endl;
				exit(EXIT_FAILURE); 
			}
			else
			{
				utilities::cout.MainOutput()<<"\tWriting lattice configuration values to file:\t\t"<<latticeConfigFile<<std::endl;
			}
			
			if(mcData->wave)
			{
				waveFile=fileName.str();
				waveFile.append("_wave.dat");
				f_wave.open(waveFile.c_str(), std::ios::binary);
				
				if(f_wave.is_open()==0)
				{
					std::cerr<<"\tERROR! cannot open file: "<<waveFile<<std::endl;
					exit(EXIT_FAILURE); 
				}
				else
				{
					utilities::cout.MainOutput()<<"\tWriting wave function values to file:\t"<<waveFile<<std::endl;
				}
			}
			
			utilities::cout.MainOutput()<<"\n\tTorus dimensions \t\t"<<wfData->Lx<<"x"<<wfData->Ly;
			
			break;
	}
	
	utilities::cout.MainOutput()<<"\n\tNo. thermalising samples\t"<<mcData->nbrTherm;
	utilities::cout.MainOutput()<<"\n\tNo. config. samples\t\t"<<mcData->nbrSamples;
	utilities::cout.MainOutput()<<"\n\tSample frequency \t\t"<<mcData->sampleFreq*wfData->nbr<<" MC steps";
	utilities::cout.MainOutput()<<"\n\tParticle moves per MC step\t"<<mcData->nbrToMove;
	utilities::cout.MainOutput()<<"\n\tRandom seed\t\t\t"<<seed<<std::endl;	
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief	Destructor for the Metropolis class
//!
////////////////////////////////////////////////////////////////////////////////

FQHE::Metropolis::~Metropolis()
{	

	//	delete dynamic memory allocation, if it has been allocated
	//	The (std::nothrow) specification on new allows us to do the !=0 test.
	
	if(uTmp!=0)			delete[]	uTmp;
	if(vTmp!=0)			delete[]	vTmp;
	if(uNew!=0)			delete[]	uNew;
	if(vNew!=0)			delete[]	vNew;
	if(thetaNew!=0) 	delete[]	thetaNew;
	if(phiNew!=0) 		delete[]	phiNew;
	if(toMove!=0) 		delete[]	toMove;
	if(phi!=0) 			delete[]	phi;
	if(theta !=0) 		delete[]	theta;
	if(u!=0) 			delete[]	u;
	if(v!=0) 			delete[]	v;
	if(zTmp!=0) 		delete[]	zTmp;
	if(zNew!=0) 		delete[]	zNew;
	if(r!=0) 			delete[]	r;
	if(z!=0) 			delete[]	z;

	//	close open files
	if(f_theta.is_open()!=0)	f_theta.close();
	if(f_phi.is_open()!=0)		f_phi.close();
	if(f_wave.is_open()!=0)		f_wave.close();
	if(f_r.is_open()!=0)		f_r.close();

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief This function performs a multi-particle random walk move on a around the
//!	surface of a sphere.
//!
//! This function moves a random point a random distance around the surface of
//!	The sphere. To move the point in an unbiased way, we must make the move in
//!	the reference frame of the north pole:
//!
//!   - Pick some point to move, with co-ordinates theta (denote it by "t_old") 
//!	  	and phi (denote it by "p_old") on a unit sphere.
//!
//!	  - We want to pick a random distance and direction to move the point in. 
//!	  	This is done by making a random 3D Cartesian vector dx with components 
//!		between -a and +a, where "a" is the maximum move distance. Then I make 
//!	    dx+(0,0,1) and normalize it. This vector is now a unit vector pointing 
//!		to a point on the surface of the sphere almost at the north pole. 
//!		I think that there will be a uniform distribution of displacements of 
//!		that vector from (0,0,1).
//!
//!	  - Next we rotate this vector to the vicinity of the point that we want to 
//!     move (I will explain why this has to be done shortly). From the co-ordinates, 
//!     {t_old,p_old} we make a rotation matrix that will take the vector (0,0,1) 
//!     and map it to the direction {t_old,p_old}. The formula I use to get this 
//!     is given here: http://en.wikipedia.org/wiki/Rotation_matrix (scroll down 
//!     to where it says: "Rotation matrix from axis and angle").
//!     Let's call this R(t_old,p_old). I then use R(t_old,p_old) to rotate the 
//!     random unit vector (dx+(0,0,1) normalized) to the location of the point 
//!     we wanted to move. This gives some new Cartesian vector (x,y,z) and then
//!     from that I can get:
//!
//!     p_new = arctan(y,x)
//!		t_new = arccos(z)
//!
//!		which are very close to the original p_old and t_old. The point is that 
//!		the value of {x,y,z} that one would get for the same dx would depend on 
//!		R(t_old,p_old). So you could not just take the original point {t_old,p_old} ,
//!  	add dx to it and normalize the result to unit length. 
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::RandomMoveSphere()
{

	//	Determine which particles are to be moved
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		*(toMove+k)=mt.genrand_int31()%wfData->nbr; 
		for(unsigned int j=0;j<k;j++)
		{
			//	Discard the random number if it has been repeated. 
			if(*(toMove+k)==*(toMove+j))
				k--;
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////
	//	THE FOLLOWING METHOD DOES NOT SEEM TO GENERATE A UNIFORM DENSITY
	#if 0
	 
	for(k=0;k<mcData->nbrToMove;k++)
	{
		//	A randomly chosen increment in the theta and phi co-ordinates
		//	See e.g. R.K. Kamilla's phd thesis, appendix B.2.4
		
		*(thetaNew+k)=*(theta+*(toMove+k))+PI*mcData->maxDist*(mt.random()-0.5);
		*(phiNew+k)=*(phi+*(toMove+k))+2*PI*mcData->maxDist*(mt.random()-0.5);
	
		//	Return the angle to the range 0<theta<PI and 0<phi<2*PI
		//	According to the procedure given in appendix B.2.4.
		
		if(-PI/2<*(thetaNew+k) && *(thetaNew+k)<0)	*(thetaNew+k)=-*(thetaNew+k);
		else if(3*PI/2>*(thetaNew+k) && *(thetaNew+k)>PI) *(thetaNew+k)=2*PI-*(thetaNew+k);
		
		if(*(phiNew+k)>2*PI)	*(phiNew+k)=*(phiNew+k)-2*PI;
		else if(*(phiNew+k)<0)	*(phiNew+k)=*(phiNew+k)+2*PI;
		
	
		//	Convert to polar and spinor co-ordinates
		*(uNew+k)=dcmplx(sqrt((1.0+cos(*(thetaNew+k)))/2)*cos(*(phiNew+k)/2),+sqrt((1.0+ cos(*(thetaNew+k)))/2)*sin(*(phiNew+k)/2));
		*(vNew+k)=dcmplx(sqrt((1.0-cos(*(thetaNew+k)))/2)*cos(*(phiNew+k)/2),-sqrt((1.0-cos(*(thetaNew+k)))/2)*sin(*(phiNew+k)/2));
		
		//	Store old values and update global values
		
		*(vTmp+k)=*(v+*(toMove+k));
		*(uTmp+k)=*(u+*(toMove+k));
		
		*(u+*(toMove+k))=*(uNew+k);
		*(v+*(toMove+k))=*(vNew+k);
		
	}
	#endif
	////////////////////////////////////////////////////////////////////////////////
	
	double uvec[3];			//		Cartesian unit vector perpendicular to axis	
							//		of rotation: takes point to north pole
	double rotn[3][3];		//		Rotation matrix to take point to north pole
	double dx[3];			//		Randomly chosen Cartesian vector at north pole
	double norm,sum;		//		Variables to generate a Cartesian unit vector
	double newxyz[3];		//		Values of the random Cartesian vector rotated
							//		back to the point to move
	double cosTheta,sinTheta;//		To store values sin(theta) and cos(theta)
	double dist;			//		Distance to move the point
	
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		//  To change reference frame, we need to determine the rotation matrix taking the point to the north pole.
		
		cosTheta=cos(*(theta+*(toMove+k)));
		sinTheta=sin(*(theta+*(toMove+k)));
		
		//	Assign unit vector of rotation axis: (comes from the cross product of (x,y,z) with (0,0,r))
		uvec[0]=-sin(*(phi+*(toMove+k)));
		uvec[1]=cos(*(phi+*(toMove+k)));
		
		//	Calculate Cartesian rotation matrix:
		rotn[0][0]=cosTheta+pow(uvec[0],2.0)*(1-cosTheta);
		rotn[0][1]=uvec[0]*uvec[1]*(1-cosTheta);
		rotn[0][2]=uvec[1]*sinTheta;
		rotn[1][0]=uvec[1]*uvec[0]*(1-cosTheta);
		rotn[1][1]=cosTheta+pow(uvec[1],2.0)*(1-cosTheta);
		rotn[1][2]=-uvec[0]*sinTheta;
		rotn[2][0]=-uvec[1]*sinTheta;
		rotn[2][1]=uvec[0]*sinTheta;
		rotn[2][2]=cosTheta;
		
		//	Pick a random Cartesian unit vector:
		dx[0]=2*mt.random()-1;
		dx[1]=2*mt.random()-1;
		dx[2]=2*mt.random()-1;
		norm=sqrt(pow(dx[0],2.0)+pow(dx[1],2.0)+pow(dx[2],2.0));
		dist=mcData->maxDist*mt.random();
		dx[0]=dist*(dx[0] /norm);
		dx[1]=dist*(dx[1] /norm);
		dx[2]=dist*(dx[2] /norm);
		
		//	Shift to north pole and renormalise
		dx[2]=dx[2]+1.0;
		norm=sqrt(pow(dx[0],2.0)+pow(dx[1],2.0)+pow(dx[2],2.0));
		dx[0]=dx[0]/norm;
		dx[1]=dx[1]/norm;
		dx[2]=dx[2]/norm;
		
		//	Rotate back to vicinity of the point
		for(int i=0;i<3;i++)
		{
			sum=0;
			for(int j=0;j<3;j++)
			{
				sum=sum+rotn[i][j]*dx[j];
			}
			newxyz[i]=sum;
		}
		
		//	Convert to polar and spinor co-ordinates
		*(phiNew+k)=atan2(newxyz[1],newxyz[0]);
		*(thetaNew+k)=acos(newxyz[2]);
		*(uNew+k)=dcmplx(sqrt((1.0+newxyz[2])/2)*cos(*(phiNew+k)/2),+sqrt((1.0+ newxyz[2])/2)*sin(*(phiNew+k)/2));
		*(vNew+k)=dcmplx(sqrt((1.0-newxyz[2])/2)*cos(*(phiNew+k)/2),-sqrt((1.0-newxyz[2])/2)*sin(*(phiNew+k)/2));
		
		//	Store old values and update global values
		
		*(vTmp+k)=*(v+*(toMove+k));
		*(uTmp+k)=*(u+*(toMove+k));
		
		*(u+*(toMove+k))=*(uNew+k);
		*(v+*(toMove+k))=*(vNew+k);
		
	}

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief 	This function moves a random point a random distance around a disc.
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::RandomMoveDisc()
{

	double dir;				//	A random move direction (angle from 0 to 2Pi)
	double dist;			//	A random move distance

	//	Determine which particles are to be moved
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		*(toMove+k)=mt.genrand_int31()%wfData->nbr; 
		for(unsigned int j=0;j<k;j++)
		{
			//	Disgard the random number if it has been repeated. 
			if(*(toMove+k)==*(toMove+j))
				k--;
		}
	}		
	
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		//	pick random direction and distance to move
		
		dir=2*PI*(mt.random())-PI;
		dist=mt.random()*mcData->maxDist*wfData->discRadius;
		
		*(zNew+k)=*(z+*(toMove+k))+dcmplx(dist*cos(dir),dist*(sin(dir)));
		
		//	Store old values and update global values
		
		*(zTmp+k)=*(z+*(toMove+k));
		
		*(z+*(toMove+k))=*(zNew+k);
	}	
	
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

////////////////////////////////////////////////////////////////////////////////
//!	\brief 	This function moves a random point a random distance for a FQHE state
//!	defined on a lattice. This method assumes a torus geometry
//!
////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::RandomMoveTorus()
{

	int rand1,rand2;
	
	//	Determine which particles are to be moved
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		*(toMove+k)=mt.genrand_int31()%wfData->nbr; 
		for(unsigned int  j=0;j<k;j++)
		{
			//	Disgard the random number if it has been repeated. 
			if(*(toMove+k)==*(toMove+j))
				k--;
		}
	}
	
	//	pick random distance to move along the x and y directions
	//	keep proposing random moves until the maximum occupancy condition is satisfied
	
	for(unsigned int k=0;k<mcData->nbrToMove;k++)
	{
		*(latticeConfigTmp+k)=*(latticeConfig+*(toMove+k));
	}
		
	bool pass=false;
		
	while(!pass)
	{
		for(unsigned int k=0;k<mcData->nbrToMove;k++)
		{
			//	propose random moves
			
			//PrintConfigurationTorus();
			
			//std::cout<<(int)(mcData->maxTorusDist/2)<<std::endl;
			
			rand1=mt.genrand_int31()%mcData->maxTorusDist;
			rand2=mt.genrand_int31()%mcData->maxTorusDist;
			
			//	random sign
			
			if(mt.genrand_int31() & 1) rand1=-rand1;
			if(mt.genrand_int31() & 1) rand2=-rand2;
			
			//std::cout<<rand1<<" "<<rand2<<std::endl;

			*(latticeConfigNew+k)=std::complex<int>((real(*(latticeConfig+*(toMove+k)))+rand1)%wfData->Lx,(imag(*(latticeConfig+*(toMove+k)))+rand2)%wfData->Ly);
			
			//*(latticeConfigNew+k)=std::complex<int>(mt.genrand_int31()%wfData->Lx,mt.genrand_int31()%wfData->Ly);

			//	update configuration
				
			*(latticeConfig+*(toMove+k))=*(latticeConfigNew+k);
			//PrintConfigurationTorus();getchar();
		}

		pass=CheckOccupancy(latticeConfig);
		
		//std::cout<<pass<<std::endl;
		
		if(!pass)
		{
			//	reset the move
			for(unsigned int k=0;k<mcData->nbrToMove;k++)
			{
				*(latticeConfig+*(toMove+k))=*(latticeConfigTmp+k);
			}
		}
	}	

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief 	This function performs a Metropolis Monte Carlo test for a FQHE
//!	state in the sphere geometry.
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::MetropolisTestSphere()
{

	if (exp(2*real(wfValNew)-2*real(wfValOld))>mt.random())
	{
		wfValOld=wfValNew;	//	Accept and update global co-ords
		
		for(unsigned int k=0;k<mcData->nbrToMove;k++)
		{
			*(theta+*(toMove+k))=*(thetaNew+k);
			*(phi+*(toMove+k))=*(phiNew+k);
		}
		//std::cout<<"accepted!"<<std::endl;
		acceptCount++;
	}
	else
	{
		//	Reject and restore old values
		for(unsigned int k=0;k<mcData->nbrToMove;k++)
		{
			*(u+*(toMove+k))=*(uTmp+k);	
			*(v+*(toMove+k))=*(vTmp+k);	
		}
		//std::cout<<"rejected!"<<std::endl;
	}

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief 	This function performs a Metropolis Monte Carlo test for a FQHE
//!	state in the disc geometry.
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::MetropolisTestDisc()
{
	if (exp(2*real(wfValNew)-2*real(wfValOld))>mt.random())
	{
		//	Accept 
		wfValOld=wfValNew;	

		//std::cout<<"accepted!"<<std::endl;
		
		acceptCount++;
	}
	else
	{
		//	Reject and restore old values
		for(unsigned int k=0;k<mcData->nbrToMove;k++)
		{
			*(z+*(toMove+k))=*(zTmp+k);	
		}
		//std::cout<<"rejected!"<<std::endl;
	}

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//!	\brief 	This function performs a Metropolis Monte Carlo test for a FQHE
//!	state in the torus geometry.
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::MetropolisTestTorus()
{
	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_==3
	
	std::cout<<"\n\tSTART Metropolis::MetropolisTestTorus()"<<std::endl;
	
	#endif
	////////////////////////////////////////////////////////////////////////////////

	if (exp(2*real(wfValNew)-2*real(wfValOld))>mt.random())
	{
		//	Accept 
		wfValOld=wfValNew;	
				
		//std::cout<<"accepted!"<<std::endl;getchar();
		acceptCount++;
	}
	else
	{
		//	Reject and restore old values
		for(unsigned int k=0;k<mcData->nbrToMove;k++)
		{
			*(latticeConfig+*(toMove+k))=*(latticeConfigTmp+k);	
		}
		//std::cout<<"rejected!"<<std::endl;
	}
	
	////////////////////////////////////////////////////////////////////////////////
	#if _DEBUG_==3
	
	std::cout<<"\n\tEND Metropolis::MetropolisTestTorus()"<<std::endl;
	
	#endif
	////////////////////////////////////////////////////////////////////////////////
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief Calculates coulomb energy in the sphere geometry using a chord
//! distance measure
//!
//////////////////////////////////////////////////////////////////////////////////

double FQHE::Metropolis::CoulombEnergy(int n, dcmplx *u, dcmplx *v)
{	
	double energy=-n/(2*wfData->discRadius);
	
	for(int i=0;i<n-1;i++)
	{
		for(int j=i+1;j<n;j++)
		{
			energy+=1.0/(n*2*wfData->sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
		}
	}
	
	return energy;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief Calculates coulomb energy in the disc geometry  in the disc geometry 
//!	(see O. Ciftja and C. Wexler PRB 67 075304 (2003))
//!
//////////////////////////////////////////////////////////////////////////////////

double FQHE::Metropolis::CoulombEnergy(int n, dcmplx *z)
{	
	//	Background-background contribution
	double energy=(double)8*n/(3*PI*wfData->discRadius);
	double dist;
	
	for(int i=0;i<n;i++)
	{
		//	Background-electron contribution
		dist=abs(*(z+i))/wfData->discRadius;
		if(dist<5)
		{	
			energy-=*(bgFunc+(int)floor(BG_FUNC_DIM*dist/5.0+0.5))*2/wfData->discRadius;
		}
		else
		{
			energy-=1/(wfData->discRadius*dist);	
		}
		
		for(int j=i+1;j<n;j++)
		{
			//electron-electron contribution
			energy+=1.0/(n*abs(*(z+i)-*(z+j)));
		}
	}
	
	return energy;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Calculates running average then resets the accumulated energy
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::CalculateRunningMean()
{
	meanEnergy=accumEnergy/energyCounter;
	stdDevEnergy=sqrt(accumEnergySquared/energyCounter-meanEnergy*meanEnergy);
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to file
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::ConfigurationToFileSphere()
{

	f_phi.write(reinterpret_cast<char*>(phi),wfData->nbr*SIZE_OF_DOUBLE);
	f_theta.write(reinterpret_cast<char*>(theta),wfData->nbr*SIZE_OF_DOUBLE);
	
	if(mcData->wave)
	{
		double realVal,imagVal;
		
		realVal=real(wfValOld);
		imagVal=imag(wfValOld);
		
		f_wave.write(reinterpret_cast<char*>(&realVal),SIZE_OF_DOUBLE);
		f_wave.write(reinterpret_cast<char*>(&imagVal),SIZE_OF_DOUBLE);
	}

}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to file
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::ConfigurationToFileDisc()
{

	for(int i=0;i<wfData->nbr;i++)
	{
		dcmplx *tmp=z+i;
		*(r+i)=abs(*tmp);
		*(theta+i)=arg(*tmp); 
	}
	
	f_r.write(reinterpret_cast<char*>(r),wfData->nbr*SIZE_OF_DOUBLE);
	f_theta.write(reinterpret_cast<char*>(theta),wfData->nbr*SIZE_OF_DOUBLE);
	
	if(mcData->wave)
	{
		double realVal,imagVal;
		
		realVal=real(wfValOld);
		imagVal=imag(wfValOld);
		
		f_wave.write(reinterpret_cast<char*>(&realVal),SIZE_OF_DOUBLE);
		f_wave.write(reinterpret_cast<char*>(&imagVal),SIZE_OF_DOUBLE);
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to file
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::ConfigurationToFileTorus()
{

	for(int i=0;i<wfData->nbr;i++)
	{
		int realVal,imagVal;
		
		realVal=real(latticeConfig[i]);
		imagVal=imag(latticeConfig[i]);
		
		f_latticeConfig.write(reinterpret_cast<char*>(&realVal),SIZE_OF_INT);
		f_latticeConfig.write(reinterpret_cast<char*>(&imagVal),SIZE_OF_INT);
	}
	
	if(mcData->wave)
	{
		double realVal,imagVal;
		
		realVal=real(wfValOld);
		imagVal=imag(wfValOld);
		
		f_wave.write(reinterpret_cast<char*>(&realVal),SIZE_OF_DOUBLE);
		f_wave.write(reinterpret_cast<char*>(&imagVal),SIZE_OF_DOUBLE);
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to stdout
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::PrintConfigurationSphere()
{
	std::cout<<"\n\tCurrent configuration is:"<<std::endl;
	for(int i=0;i<wfData->nbr;i++)
	{
		std::cout<<"\tu["<<i<<"]="<<u[i]<<"\t"<<"v["<<i<<"]="<<v[i]<<std::endl;
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to stdout
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::PrintConfigurationDisc()
{
	std::cout<<"\n\tCurrent configuration is:"<<std::endl;
	for(int i=0;i<wfData->nbr;i++)
	{
		std::cout<<"\tz["<<i<<"]="<<z[i]<<std::endl;
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief	Print current configuration to stdout
//!
//////////////////////////////////////////////////////////////////////////////////

void FQHE::Metropolis::PrintConfigurationTorus()
{
	std::cout<<"\n\tCurrent configuration is:"<<std::endl;
	for(int i=0;i<wfData->nbr;i++)
	{
		std::cout<<"\tlatticeConfig["<<i<<"]="<<latticeConfig[i]<<std::endl;
	}
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//

//////////////////////////////////////////////////////////////////////////////////
//! \brief Function to check that the maximum occupancy condition is satisfied
//!	for Metropolis moves on a lattice. 
//!
//!	\return true if the condition is satisfied and false if it is not
//!
//////////////////////////////////////////////////////////////////////////////////

bool FQHE::Metropolis::CheckOccupancy(
	std::complex<int>* config)		//!<	A lattice occupancy configuration
{

	//	make a list of all the current sites which are occupied
	
	std::vector<std::complex<int> > list;
	std::vector<std::complex<int> >::iterator k;
	
	int i,j,counter=1,flag;
	dcmplx element;

	list.push_back(config[0]);
	
	for(i=1;i<wfData->nbr;i++)
	{
		element=config[i];
		
		for(j=0;j<counter;j++)
		{
			flag=0;
			
			if(real(element)==real(list[j]) && imag(element)==imag(list[j]))
			{
				flag=1;
				break;
			}
		}
		
		if(flag==0)
		{	
			//	if element is not in list then add element to the list
			list.push_back(element);
			
			//	counter keeps track of the list length
			counter++;
		}
	}
	
	int occupations[counter];
	
	for(j=0;j<counter;j++)
	{
		occupations[j]=0;
	}
	
	for(i=0;i<wfData->nbr;i++)
	{
		element=config[i];
		
		for(j=0;j<counter;j++)
		{
			if(real(element)==real(list[j]) && imag(element)==imag(list[j]))
			{	
				occupations[j]++;
			}
		}
	}
	
	//for(j=0;j<counter;j++)
	//{
	//	std::cout<<list[j]<<" "<<occupations[j]<<std::endl;
	//}
	
	for(j=0;j<counter;j++)
	{
		if(occupations[j]>wfData->maxOccupation)
		{
			return false;
		}
	}

	return true;
}

//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\//
