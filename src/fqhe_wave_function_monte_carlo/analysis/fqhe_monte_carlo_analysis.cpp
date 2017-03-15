////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file 
//!		This is the library file contains a number of possible electrostatic
//!		potentials to consider in a FQHE system
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

#include "fqhe_monte_carlo_analysis.hpp"

//!
//! Default constructor for analysis options - sets default values
//!
AnalysisOptions::AnalysisOptions()
    :
	pairCorrel(false),
	dens(false),
	coulomb(false),
	secllCoulomb(false),
	finiteThickness(false),
	bilayer(false),
	chargePlate(false),
	biChargePlate(false),
	partialCoulomb(false),
	resample(false),		    
	nbrResamples(1),	        
	dataDirIn("./"),	        
	dataDirOut("./"),           
	skip(0),			        
	check(false),
	nPairBins(500),		        
	nDensBins(500),		        
	secllType("d"),
	finiteType("h"),
	minThickness(0.01),
	maxThickness(5.0),
	thicknessStepSize(0.01),
	minBilayer(0.01),
	maxBilayer(5.0),
	bilayerStepSize(0.01),
	minPlate(0.01),
	maxPlate(5.0),
	plateStepSize(0.01),
	partialCoulombType(1),
    getFileName(false)
    {}

//!
//! Initialize the analysis options struct variables from command line
//! arguments
//!
void AnalysisOptions::InitFromCommandLine(
    boost::program_options::variables_map* vm)  //!<    Program options map
{
    this->dataDirIn   = (*vm)["input"].as<std::string>();
    this->dataDirOut  = (*vm)["output"].as<std::string>();
    this->getFileName = (*vm)["get-file-name"].as<bool>();
    this->skip        = (*vm)["skip"].as<int>();
    this->check       = (*vm)["check"].as<bool>();
    this->pairCorrel  = (*vm)["pair"].as<bool>();
    this->nPairBins   = (*vm)["pair-bins"].as<int>();
    this->dens        = (*vm)["dens"].as<bool>();
    this->nDensBins   = (*vm)["dens-bins"].as<int>();
    this->coulomb     = (*vm)["coulomb"].as<bool>();
    this->secllCoulomb= (*vm)["second"].as<bool>();
    this->secllType   = (*vm)["second-type"].as<std::string>();
    this->resample    = (*vm)["re-sample"].as<bool>();
    this->nbrResamples = (*vm)["nbr-re-sample"].as<int>();
    this->finiteThickness = (*vm)["finite"].as<bool>();
    this->finiteType  = (*vm)["finite-type"].as<std::string>();
    this->minThickness= (*vm)["min-thickness"].as<double>();
    this->maxThickness= (*vm)["max-thickness"].as<double>();
    this->thicknessStepSize = (*vm)["step-thickness"].as<double>();
    this->bilayer     = (*vm)["bilayer"].as<bool>();
    this->minBilayer  = (*vm)["min-bilayer"].as<double>();
    this->maxBilayer  = (*vm)["max-bilayer"].as<double>();
    this->bilayerStepSize = (*vm)["step-bilayer"].as<double>();
    this->chargePlate = (*vm)["charge-plate"].as<bool>();
    this->biChargePlate = (*vm)["plate-bilayer"].as<bool>();
    this->minPlate    = (*vm)["min-plate"].as<double>();
    this->maxPlate    = (*vm)["max-plate"].as<double>();
    this->plateStepSize = (*vm)["step-plate"].as<double>();
    this->partialCoulomb = (*vm)["partial-coulomb"].as<bool>();
    this->partialCoulombType = (*vm)["partial-type"].as<int>();
}

//!
//! Default constructor
//!
AnalysisMethods::AnalysisMethods()
    :
    #if _ENABLE_ANALYSIS_METHODS_
    pairFreqUp(0),
    pairFreqDown(0),
    pairFreqDiff(0),
    densityFreq(0),
    f_pair(0), 
    f_dens(0),          
    f_background(0),
    discBgFunc(0),
    thicknessBg(0),
    secondBg(0),        
    biChargePlateBg(0),
    partialBg(0),
    finiteThicknessPotential(0),
    f_potential(0),     
    f_coeffs(0),
    coeffList(0),
	partialType(""),    
	f_energy(0),        
	f_second(0),        
	f_thickness(0),     
	f_bilayer(0),       
	f_plate(0),         
	f_biChargePlate(0), 
	f_partialCoulomb(0),
    coulombEnergy(0),		
	secondEnergy(0),	
	partialEnergy(0),   
	thicknessEnergy(0), 
	bilayerEnergy(0),   
	chargePlateEnergy(0), 	                    
	biChargePlateEnergy(0),                   
	coulombStdev(0),    
	secondStdev(0),     
	partialStdev(0),    
	coulombMeanEnergy(0),                   
	secondMeanEnergy(0),
	partialMeanEnergy(0),
    #endif
    u(0),               
    v(0),               
    z(0),               
    phi(0),             
    theta(0),
    latticeConfig(0),             
    f_theta(0),         
    f_phi(0),           
    f_r(0),
    f_config(0),        
    maxRuns(0),
    nbrSamples(0),
    wf(0)          
{}

//!
//! Constructor from analysis options
//!
AnalysisMethods::AnalysisMethods(
    AnalysisOptions& options)    //!<    Address of analysis options struct
    :
    #if _ENABLE_ANALYSIS_METHODS_
    pairFreqUp(0),
    pairFreqDown(0),
    pairFreqDiff(0),
    densityFreq(0),
    f_pair(0),
    f_dens(0),
    f_background(0),
    discBgFunc(0),
    thicknessBg(0),
    secondBg(0),
    biChargePlateBg(0),
    partialBg(0),
    finiteThicknessPotential(0),
    f_potential(0),
    f_coeffs(0),
    coeffList(0),
	partialType(""),
	f_energy(0),
	f_second(0),
	f_thickness(0),
	f_bilayer(0),
	f_plate(0),
	f_biChargePlate(0),
	f_partialCoulomb(0),
	options(options),
    coulombEnergy(0),
	secondEnergy(0),
	partialEnergy(0),
	thicknessEnergy(0),
	bilayerEnergy(0),
	chargePlateEnergy(0),
	biChargePlateEnergy(0),
	coulombStdev(0),
	secondStdev(0),
	partialStdev(0),
	coulombMeanEnergy(0),
	secondMeanEnergy(0),
	partialMeanEnergy(0),   
    #endif
    u(0),
    v(0),
    z(0),
    phi(0),
    theta(0),
    latticeConfig(0),
    f_theta(0),
    f_phi(0),
    f_r(0),
    f_config(0),
    maxRuns(0),
    nbrSamples(0),
    wf(0)
{}

//!
//! Destructor
//!
AnalysisMethods::~AnalysisMethods()
{
	if(u!=0)				delete[] u;
	if(v!=0)				delete[] v;
	if(z!=0)				delete[] z;
	if(phi!=0)				delete[] phi;
	if(theta!=0)			delete[] theta;
	if(latticeConfig!=0)    delete[] latticeConfig;
	if(wf!=0)               delete wf;
	if(f_theta.is_open()!=0)	f_theta.close();
	if(f_phi.is_open()!=0)		f_phi.close();
	if(f_r.is_open()!=0)		f_r.close();
	if(f_config.is_open()!=0)	f_config.close();
	
}

//!
//! Initialize analysis methods options from command line arguments
//!
void AnalysisMethods::InitFromCommandLine(
    #if _ENABLE_MPI_
    boost::program_options::variables_map* vm,  //!<    Program options map
    const utilities::MpiWrapper& mpi)           //!<    Address of mpi wrapper class
    #else
    boost::program_options::variables_map* vm)  //!<    Program options map
    #endif
{
    #if _ENABLE_MPI_
    wfData.InitFromCommandLine(vm,mpi);
    cfData.InitFromCommandLine(vm,mpi);
    #else
    wfData.InitFromCommandLine(vm);
    cfData.InitFromCommandLine(vm);
    #endif
    options.InitFromCommandLine(vm);
    if(wfData.type==FQHE::_LAUGHLIN_)
	{
		wf = new FQHE::Laughlin(&wfData);
	}
	else if(wfData.type==FQHE::_MOORE_READ_)
	{	
		wf = new FQHE::MooreRead(&wfData);
	}
	else if(wfData.type==FQHE::_NASS_)
	{	
		wf = new  FQHE::NonAbelianSpinSinglet(&wfData);
	}
	else if(wfData.type==FQHE::_COMPOSITE_FERMION_ || wfData.type==FQHE::_BONDERSON_SLINGERLAND_)
	{
		wf = new FQHE::CompositeFermion(&wfData,&cfData);
	}
	else if(wfData.type==FQHE::_PARAFERMION_)
	{
	    wf = new FQHE::Parafermion(&wfData);
	}
	if(wfData.geometry==FQHE::_SPHERE_)
	{
		phi   = new (std::nothrow) double[wfData.nbr];
		theta = new (std::nothrow) double[wfData.nbr];
		u     = new (std::nothrow) std::complex<double>[wfData.nbr];
		v     = new (std::nothrow) std::complex<double>[wfData.nbr];
	}
	else if(wfData.geometry==FQHE::_DISC_)
	{
		z = new (std::nothrow) std::complex<double>[wfData.nbr];
	}
	else if(wfData.geometry==FQHE::_TORUS_)
	{
		latticeConfig = new (std::nothrow) std::complex<int>(wfData.nbr);
	}
}

//!
//! Print the file name and generate a text file containing 
//! that information
//!
void AnalysisMethods::PrintName()
{
	std::cout << std::endl << "\tFile Name:\t\t" << wfData.wfFileName << std::endl;
	std::ofstream f_out;
	f_out.open("fileName.txt");
	f_out << wfData.wfFileName;
	f_out.close();
}

//!
//! A function to deal with generating file names and opening files
//!	file names depend on wf parameters, directories can be specified
//!
int AnalysisMethods::OpenConfigurationDataFile(
    std::string dataDirIn,      //!<    Name of directory to read in from
    const int runNo)            //!<    Integer label of the "run" file to read in 
{
	if(f_theta.is_open()!=0)	f_theta.close();
	if(f_phi.is_open()!=0)		f_phi.close();
	if(f_r.is_open()!=0)		f_r.close();
	if(f_config.is_open()!=0)	f_config.close();
	//	Open new configuration data file
	std::stringstream fileName;
	fileName.str("");
	if(dataDirIn!="")  
	{
		fileName << dataDirIn;
		fileName << "/";
	}
	fileName << wfData.wfFileName << "_run_" << runNo;
	if(wfData.geometry==FQHE::_SPHERE_)
	{
	    std::string phiFile;
		phiFile = fileName.str();
		phiFile.append("_phi.dat");
		f_phi.open(phiFile.c_str(), std::ios::binary);
		if(f_phi.is_open()==0)
		{
			std::cout << "\n\tCannot open file: " << phiFile << std::endl << std::endl;
			return EXIT_FAILURE;
		}
		else
		{
			std::cout << "\n\tReading phi values from file:\t" << phiFile << std::endl;
		}
		std::string thetaFile;
		thetaFile = fileName.str();
		thetaFile.append("_theta.dat");
		f_theta.open(thetaFile.c_str(), std::ios::binary);
		if(f_theta.is_open()==0)
		{
			std::cout << "\n\tCannot open file: " << thetaFile << std::endl << std::endl;
			return EXIT_FAILURE; 
		}
		else
		{
			std::cout << "\tReading theta values from file:\t" << thetaFile << std::endl;
		}
	}
	else if(wfData.geometry==FQHE::_DISC_)
	{
	    std::string rFile;
		rFile = fileName.str();
		rFile.append("_r.dat");
		f_r.open(rFile.c_str(), std::ios::binary);
		if(f_r.is_open()==0)
		{
			std::cout << "\n\tCannot open file: " << rFile << std::endl << std::endl;
			return EXIT_FAILURE; 
		}
		else
		{
			std::cout << "\n\tReading r values from file:\t" << rFile << std::endl;
		}
		std::string thetaFile;
		thetaFile = fileName.str();
		thetaFile.append("_theta.dat");
		f_theta.open(thetaFile.c_str(), std::ios::binary);
		if(f_theta.is_open()==0)
		{
			std::cout << "\n\tCannot open file: " << thetaFile << std::endl << std::endl;
			return EXIT_FAILURE;
		}
		else
		{
			std::cout << "\tReading theta values from file:\t" << thetaFile << std::endl;
		}
	}
	else if(wfData.geometry==FQHE::_TORUS_)
	{
	    std::string configFile;
		configFile = fileName.str();
		configFile.append(".dat");
		f_config.open(configFile.c_str(), std::ios::binary);
		if(f_config.is_open()==0)
		{
			std::cout << "\n\tCannot open file: " << configFile << std::endl << std::endl;
			return EXIT_FAILURE; 
		}
		else
		{
			std::cout << "\n\tReading configuration values from file:\t" << configFile << std::endl;
		}
	}
	return 0;
}

//!
//! Attempt to open all relavant configuration data files and count the
//!	number of sample configurations containined in those files
//!
int AnalysisMethods::LookForConfigurationDataFiles()
{
	size_t size=0;		//	Memory size to be allocated to store data on each sample
	int exitFlag=0;		//	Flag to cause the program to exit if a file cannot be opened
	std::cout << "\n\n\tSearching for configuration data files...\n";
	this->maxRuns = 1;
	exitFlag = this->OpenConfigurationDataFile(this->options.dataDirIn, this->maxRuns);
	if(exitFlag==EXIT_FAILURE) 
	{
	    return EXIT_FAILURE;
	}
	while(exitFlag==0)
	{
		if(wfData.geometry==FQHE::_TORUS_)
		{
			f_config.seekg(0, std::ios::end);
			size += f_config.tellg();
			size -= options.skip*(wfData.nbr*sizeof(std::complex<int>));
		}
		else
		{
			f_phi.seekg(0, std::ios::end);
			size += f_phi.tellg();
			size -= options.skip*(wfData.nbr*sizeof(double));
		}
		exitFlag = OpenConfigurationDataFile(options.dataDirIn, this->maxRuns+1);
		if(exitFlag!=0) 
		{
		    break;
		}
		++maxRuns;
	}
	if(f_theta.is_open()!=0)	f_theta.close();
	if(f_phi.is_open()!=0)		f_phi.close();
	if(f_config.is_open()!=0)	f_config.close();
	if(wfData.geometry==FQHE::_TORUS_)
	{
		nbrSamples = size/(wfData.nbr*sizeof(std::complex<int>));
	}
	else
	{
		nbrSamples = size/(wfData.nbr*sizeof(double));
	}
	std::cout << "\t...detected " << nbrSamples << " configuration samples in total." << std::endl;
	return 0;
}

//!
//! Read in a single configuration from a file (sphere geometry)
//!
void AnalysisMethods::GetConfigFromFileSphere()
{
	f_phi.read(reinterpret_cast<char*>(this->phi), this->wfData.nbr*SIZE_OF_DOUBLE);
	f_theta.read(reinterpret_cast<char*>(this->theta), this->wfData.nbr*SIZE_OF_DOUBLE);
}

//!
//! Read in a single configuration from a file (disc geometry)
//!
void AnalysisMethods::GetConfigFromFileDisc()
{
	double r[wfData.nbr];
	double theta[wfData.nbr];
	f_r.read(reinterpret_cast<char*>(r), this->wfData.nbr*SIZE_OF_DOUBLE);
	f_theta.read(reinterpret_cast<char*>(theta), this->wfData.nbr*SIZE_OF_DOUBLE);
	for(int i=0; i<this->wfData.nbr; ++i)
	{
		this->z[i]=dcmplx(r[i]*cos(theta[i]),r[i]*sin(theta[i]));
	}
}

//!
//! Read in a single configuration from a file (torus geometry)
//!
void AnalysisMethods::GetConfigFromFileTorus()
{
	int realVal,imagVal;
	for(int i=0; i<this->wfData.nbr; ++i)
	{
		f_config.read(reinterpret_cast<char*>(&realVal), SIZE_OF_INT);
		f_config.read(reinterpret_cast<char*>(&imagVal), SIZE_OF_INT);
		latticeConfig[i]=std::complex<int>(realVal, imagVal);
	}
}

//!
//! Convert theta,phi co-ordinates to u,v co-ordinates
//!
void AnalysisMethods::PolarToSpinor()
{
	for (int i=0; i<this->wfData.nbr; ++i)
	{
		this->u[i]=dcmplx(sqrt((1.0+cos(this->theta[i]))/2)*cos(this->phi[i]/2), +sqrt((1.0+cos(this->theta[i]))/2)*sin(this->phi[i]/2));
		this->v[i]=dcmplx(sqrt((1.0-cos(this->theta[i]))/2)*cos(this->phi[i]/2), -sqrt((1.0-cos(this->theta[i]))/2)*sin(this->phi[i]/2));
	}
}

//!
//! Check that particles did not get stuck at the north pole
//!
void AnalysisMethods::CheckConfig(
    const int k)      //!<    MC sample index
{
	int flag = 0;
    int flag2 = 0;
	for(int i=0; i<wfData.nbr; ++i)
	{
		if(cos(theta[i])>0.999999)
		{
			std::cout << "\tWarning: Particle " << k << ":" << i << " close to pole (z/r is " << cos(theta[i]);
			std::cout << ") check that it didn't get stuck!" << std::endl << "\tTracking position:" << std::endl;
			flag += 4;
			flag2 = i;
		}
		if(i==flag2 && cos(theta[i])<=0.999999 && flag!=0)
		{
			std::cout << "\t" << k << ":" << i << " (z/r is " << cos(theta[i]) << ")" << std::endl;				
			flag--;
			if(flag==0)
			{
			    std::cout << "\tOK: Particle is now unstuck!" << std::endl;
			}
		}
	}
}

//!
//! This function will resample the energy value if requested, 
//!	in order to determine A precise estimate of the standard deviation
//!
void AnalysisMethods::Resample(
    double *energyStore,    //!<    A buffer to store resampled energies
    const int nbrResamples, //!<    The number of resamples to perform  
    double &stDev,          //!<    Address of standard deviation value to be 
                            //!     returned
    double &meanEnergy)     //!<    Address of mean energy value to be returned
{
	double timer;				    //	Keep track of the time taken
	double runningResample;		    //	resampled energy
	double resampleav,resampleav2;  // <E> and <E^2> of the resampled energy
	int seed;
	if(getenv("PBS_JOBID")==NULL)	seed=time(0);
	else							seed=atoi(getenv("PBS_JOBID"));
	srand(seed);
	long unsigned int initArray[INIT_ARRAY_SIZE];
	for(int i=0; i<INIT_ARRAY_SIZE; ++i)
	{
		initArray[i] = rand();
	}
	MersenneTwister mt;
	mt.init_by_array(initArray,INIT_ARRAY_SIZE);
	// Re-sample energy data nbrResamples times:
	timer = clock();
	resampleav = 0;
	resampleav2 = 0;
	// Hash table to hold resample indices (memory pre-allocated). 
    // Up to one repeated index is allowed (but not more than one)
	#if _SPEED_OPTIMIZED_MAP_
    google::dense_hash_map<unsigned int,bool,utilities::MurmurHasher64Wrapper<unsigned int>> indices;
    google::dense_hash_map<unsigned int,bool,utilities::MurmurHasher64Wrapper<unsigned int>> duplicatedIndices; 
    #elif _MEMORY_OPTIMIZED_MAP_
    google::sparse_hash_map<unsigned int,bool,utilities::MurmurHasher64Wrapper<unsigned int>> indices;
    google::sparse_hash_map<unsigned int,bool,utilities::MurmurHasher64Wrapper<unsigned int>> duplicatedIndices;
    #else
    std::unordered_map<unsigned int,bool> indices;
	std::unordered_map<unsigned int,bool> duplicatedIndices;
    #endif
	#if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
	indices.set_deleted_key(std::numeric_limits<unsigned int>::max());
	duplicatedIndices.set_deleted_key(std::numeric_limits<unsigned int>::max());
	#endif
	#if _SPEED_OPTIMIZED_MAP_
	indices.set_empty_key(std::numeric_limits<unsigned int>::max()-1);
	duplicatedIndices.set_empty_key(std::numeric_limits<unsigned int>::max()-1);
	#endif
	#if _SPEED_OPTIMIZED_MAP_ || _MEMORY_OPTIMIZED_MAP_
    indices.resize(nbrSamples);
    duplicatedIndices.resize(nbrSamples);
    #else
    indices.reserve(nbrSamples);
    duplicatedIndices.reserve(nbrSamples);
    #endif
	for(int j=0; j<nbrResamples; ++j)
	{
	    indices.clear();
	    duplicatedIndices.clear();
	    while(indices.size()+duplicatedIndices.size() < (unsigned int)nbrSamples)
		{
		    int random = mt.genrand_int31() % nbrSamples;
		    if(!duplicatedIndices[random])
		    {   
		        duplicatedIndices[random] = true;
		    }
		    else
		    {
		        indices[random] = true;
		    }
		}
		runningResample=0.0;
		for(auto& x:indices)
		{
		    runningResample += *(energyStore+x.first);
		}
		
		for(auto& x:duplicatedIndices)
		{
		    runningResample += *(energyStore+x.first);
		}
		runningResample /= nbrSamples;
		resampleav  += runningResample;
		resampleav2 += runningResample*runningResample;	
		if((j+1)%(100)==0)
		{
			std::cout << "\t" << (double)(100*(j+1))/nbrResamples << "% completed\t" << "Time remaining (minutes):\t";
			std::cout << (double)1/60*(((double)nbrResamples/(j+1))-1)*( double ) ( clock() - timer ) / CLOCKS_PER_SEC << std::endl;
		}
	}
	resampleav /= nbrResamples;
	resampleav2 /= nbrResamples;
	stDev = sqrt(resampleav2-resampleav*resampleav);
	meanEnergy = resampleav;
	return;
}

#if _ENABLE_ANALYSIS_METHODS_
//!
//!	A function to deal with the relevant background data, based on the 
//!	type of analysis required e.g. read in finite-thickness background and potential
//!	Also, initializes array structures for data storage and constructs file names
//!
void AnalysisMethods::InitializeMethods()
{
    std::stringstream fileName;
	if(nbrSamples<=0)
	{
		std::cout << "\tERROR: cannot initialize methods if no configuration data have been located!\n";
		exit(EXIT_FAILURE);
	}
	if(options.secllCoulomb)
	{
		if(options.secllType=="d")	options.secllType="disc";
		if(options.secllType=="g")	options.secllType="graphene";
		if(options.secllType=="s")	options.secllType="sphere";
	}
	std::cout << "\n\tCalculating pair correlation function\t " << ((options.pairCorrel) ? ("on") : ("off")) << std::endl;
	if(options.pairCorrel)	std::cout << "\tNo. bins\t\t\t\t" << options.nPairBins << std::endl;
	std::cout << "\tCalculating density function\t\t " << ((options.dens) ? ("on") : ("off")) << std::endl;
	if(options.dens)	std::cout << "\tNo. bins\t\t\t\t" << options.nDensBins << std::endl;
	std::cout << "\n\tCalculating with the following potentials (separately):\n\n";
	std::cout << "\tCoulomb energy\t\t\t\t " << ((options.coulomb) ? ("on") : ("off")) << std::endl;
	std::cout << "\t2nd LL Coulomb energy\t\t\t " << ((options.secllCoulomb) ? ("on") : ("off")) << std::endl;
	if(options.secllCoulomb)	std::cout << "\t2nd LL type:\t\t\t\t" << options.secllType << std::endl;
	std::cout << "\tFinite thickness potential\t\t " << ((options.finiteThickness) ? ("on") : ("off")) << std::endl;
	if(options.finiteThickness)
	{
		if(options.finiteType=="h")		options.finiteType="heterojunction";
		if(options.finiteType=="w")		options.finiteType="well";
		std::cout << "\tPotential type:\t\t\t\t " << options.finiteType << std::endl;
		std::cout << "\tThickness parameter from " << options.minThickness << " to " << options.maxThickness;
		std::cout << " (steps of " << options.thicknessStepSize << ")" << std::endl;
	}
	std::cout << "\tBilayer potential energy\t\t " << ((options.bilayer) ? ("on") : ("off")) << std::endl;
	if(options.bilayer)
	{
		std::cout << "\tBilayer parameter from " << options.minBilayer << " to " << options.maxBilayer;
		std::cout << " (steps of " << options.bilayerStepSize << ")" << std::endl;
	}
	std::cout << "\tCharge plate potential energy\t\t " << ((options.chargePlate) ? ("on") : ("off")) << std::endl;
	if(options.chargePlate)
	{	
		std::cout << "\tBilayer parameter from " << options.minPlate << " to " << options.maxPlate;
		std::cout << " (steps of " << options.plateStepSize << ")" << std::endl;
	}
	std::cout << "\tCombined bilayer+charge plate potential\t " << ((options.biChargePlate) ? ("on") : ("off")) << std::endl;
	if(options.biChargePlate)
	{	
		std::cout << "\tBilayer parameter from 0 to 10 in steps of 0.4" << std::endl;
		std::cout << "\tCharge plate parameter from 0 to 10 in steps of 0.5" << std::endl;
	}
	std::cout << "\t'Partial Coulomb' potential energy\t " << ((options.partialCoulomb) ? ("on") : ("off")) << std::endl;
	if(options.partialCoulomb)
	{
		switch(options.partialCoulombType)
		{
			case 1:
				std::cout << "LLL Coulomb (up-up pairs only)";
				break;
			case 2:
				std::cout << "LLL Coulomb (up-down pairs only)";
				break;
			case 3:
				std::cout << "2nd LL Coulomb (up-up pairs only)";
				break;
			case 4:
				std::cout << "2nd LL Coulomb (up-down pairs only)";
				break;
			case 5:
				std::cout << "(graphene) 2nd LL Coulomb (up-up pairs only)";
				break;
			case 6:
				std::cout << "(graphene) 2nd LL Coulomb (up-down pairs only)";
				break;
		}
		if(options.partialCoulombType>6 && options.partialCoulombType<=13)
		{
			std::cout << "Effective potential (sphereRadius)^(" << options.partialCoulombType-6 
			          << ")*exp(-sphereRadius) (up-up pairs only)";
		}
		else if(options.partialCoulombType>13 && options.partialCoulombType<=20)
		{
			std::cout << "Effective potential (sphereRadius)^(" << options.partialCoulombType-13 
			          << ")*exp(-sphereRadius) (up-down pairs only)";
		}
	}
	std::cout << std::endl << "\tGenerating output data files..." << std::endl << std::endl;
	if(options.pairCorrel)
	{
		//	initialize pair correlation function calculation
		this->pairFreqUp   = new (std::nothrow) double[options.nPairBins+1];
		this->pairFreqDown = new (std::nothrow) double[options.nPairBins+1];
		this->pairFreqDiff = new (std::nothrow) double[options.nPairBins+1];
		for(int i=0; i<options.nPairBins+1; ++i)
		{
			this->pairFreqUp[i]=0;
			if(wfData.nbrDown>0)
			{
				this->pairFreqDiff[i]=0;
				this->pairFreqDown[i]=0;
			}
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_pair.dat";
		this->f_pair.open(fileName.str().c_str());
		if(this->f_pair.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	if(options.dens)
	{
		this->densityFreq = new (std::nothrow) double[options.nDensBins+1];
		for(int i=0; i<options.nDensBins+1; ++i)
		{
			this->densityFreq[i]=0;
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_dens.dat";
		this->f_dens.open(fileName.str().c_str());
		if(this->f_dens.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	if(options.coulomb)
	{
		coulombEnergy = new double[nbrSamples];
		std::cout << "\n\tAllocated " << (double)nbrSamples*sizeof(double)/(1024*1024) 
		          << " mb to store coulomb energy values"<<std::endl<<std::endl;
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_energy.dat";
		f_energy.open(fileName.str().c_str());
		if(f_energy.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	if(options.secllCoulomb)
	{
		//  TODO compute 2nd LL background integral from scratch
		secondEnergy = new double[nbrSamples];
		std::cout << "\n\tAllocated " << (double)nbrSamples*sizeof(double)/(1024*1024)
		          << " mb to store second LL coulomb energy values" << std::endl << std::endl;
		fileName.str("");
		fileName << "second_ll_" << options.secllType << "_background.dat";	
		f_background.open(fileName.str().c_str());
		if(f_background.is_open()==0)
		{
			std::cerr << "\tERROR: cannot open background data file:" << fileName.str() << std::endl;
			exit (EXIT_FAILURE);
		}
		//	The file is arranged so that the background associated with flux q is stored at location q
		//	Cycle through the list to obtain the correct value
		for(int j=0; j<(int)wfData.flux; ++j) 
		{
			f_background>>this->secondBg;
		}
		//	Background values are storred without the factor of n^2: we need to put this back in!
		this->secondBg *= wfData.nbr*wfData.nbr;
		std::cout << "\tUsing 2nd LL Background: " << this->secondBg << std::endl;
		f_background.close();
		coeffList = new double[7];
		//  TODO compute pseudopotential coefficients from scratch
		if(options.secllType == "disc")
		{
			coeffList[0] = -50.36597363;
			coeffList[1] = 87.38179510;
			coeffList[2] = -56.08455086;
			coeffList[3] = 17.76579124;
			coeffList[4] = -2.971636200;
			coeffList[5] = 0.2513169758;
			coeffList[6] = -0.008434843187;
		}
		else if(options.secllType == "graphene")
		{
			//	These values are taken from PRB 75 245440 (2007)
			coeffList[0] = -11.1534;
			coeffList[1] = 24.9078;
			coeffList[2] = -18.6461;
			coeffList[3] = 6.63657;
			coeffList[4] = -1.221097;
			coeffList[5] = 0.1122068;
			coeffList[6] = -0.00404269;
		}
		else if(options.secllType == "sphere")
		{
			//	Need to read in values from external file
			fileName.str("");
			fileName << "2nd_LL_sphere_potential/flux_" << wfData.flux <<"_potential.dat";	
			f_coeffs.open(fileName.str().c_str());
			if(f_coeffs.is_open()==0)
			{
				std::cerr << "\tERROR: cannot open file:" << fileName.str() << std::endl;
				exit(EXIT_FAILURE);
			}
			std::cout << "\n\tList of coefficients in effective potential:" << std::endl;
			for(int i=0; i<7; ++i)		
			{
				f_coeffs >> coeffList[i];
			}
			f_coeffs.close();
			for(int i=0; i<7; ++i)
			{
				std::cout << *(i+coeffList) << std::endl;
			}
		}
		else 
		{
			std::cerr << "\tERROR: unknown second LL type: " << options.secllType << std::endl;
			exit(EXIT_FAILURE);
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_second_" << options.secllType << ".dat";
		f_second.open(fileName.str().c_str());
		if(f_second.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	//	If finite thickness potential is selected then we need to read in the appropriate background values
	if(options.finiteThickness)
	{
	    //  TODO compute finite thickness background potential from scratch
	    //  and for a variable range of parameters
		thicknessBg = new double[THICKNESS_BG_SIZE];
		//	Read in lookup table of values of the finite thickness potential
		finiteThicknessPotential = new double[THICKNESS_POTENTIAL_FILE_SIZE];	
		fileName.str("");
		fileName << "background_data/" << "flux_" << wfData.flux << "_" << options.finiteType << "_background.dat";
		f_background.open(fileName.str().c_str());
		if (f_background.is_open()==0)
		{
			std::cerr << "\tERROR: cannot open finite thickness background energy file: " << fileName.str() << std::endl;
			exit(EXIT_FAILURE);
		}
		fileName.str("");
		fileName << "finite_thickness_potential/finite_thickness_" << options.finiteType << ".dat";
		f_potential.open(fileName.str().c_str());
		if(f_potential.is_open()==0)
		{
			std::cerr << "\tERROR: cannot open finite thickness potential file: " << fileName.str() << std::endl;
			exit(EXIT_FAILURE);
		}
		// the file stores THICKNESS_BG_SIZE=1000 values of the background potential, 
		// for finite thickness parameter = 0.01 to 10 in steps of 0.01
		for(int j=0; j<THICKNESS_BG_SIZE; ++j)
		{
			f_background>>thicknessBg[j];
			//background values are storred without the factor of n^2: we need to put this back in!
			thicknessBg[j] *= wfData.nbr*wfData.nbr;
		}
		for(int j=0; j<THICKNESS_POTENTIAL_FILE_SIZE; ++j)
		{
			f_potential >> finiteThicknessPotential[j];
		}
		f_background.close();
		f_potential.close();
		std::cout << "\n\tUsing " << options.finiteType << " type finite-thickness correction to calculate electrostatic energy.\n" << std::endl;
		//	Put a maximum limit on the finite thickness parameter range (due to limited background data)
		if(options.minThickness<0)
		{
			std::cerr << "\tERROR: options.minThickness " << options.minThickness << " out of range." << std::endl;
			exit(EXIT_FAILURE);
		}		
		if(options.maxThickness>10)
		{
			std::cerr << "\tERROR: options.maxThickness " << options.maxThickness << " out of range\n"
			          << "\t(insufficient backgorund data provided)" << std::endl;
			exit(EXIT_FAILURE);
		}
		if((int)(options.thicknessStepSize*THICKNESS_BG_SIZE)%10!=0)
		{
			std::cerr << "\tERROR: thickness parameter step size " << options.thicknessStepSize << " not possible\n"
			          << "\t(insufficient backgorund data provided)" << std::endl;
			exit(EXIT_FAILURE);
		}
		//	Determine the number of energy values we will need to store 
		//  in the required range
		{
		    int nEnergyValues = 0;
		    double thicknessValue=options.minThickness;
		    while(thicknessValue<=options.maxThickness)
		    {
			    ++nEnergyValues;
			    thicknessValue+=options.thicknessStepSize;
		    }
		    thicknessEnergy = new double[nEnergyValues];
		    for(int i=0; i<nEnergyValues; ++i)
		    {
			    thicknessEnergy[i] = 0;
		    }
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_" << options.finiteType << "_";
		fileName << options.minThickness << "_" << options.maxThickness << ".dat";
		f_thickness.open(fileName.str().c_str());
		if(f_thickness.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	//	If a options.bilayer potential is selected, we need to allocate memory to store energy values
	if(options.bilayer)
	{
		int nEnergyValues = 0;
		double bilayerValue = options.minBilayer;
		while(bilayerValue <= options.maxBilayer)
		{
			++nEnergyValues;
			bilayerValue += options.bilayerStepSize;
		}
		bilayerEnergy = new double[nEnergyValues];
		for(int i=0; i<nEnergyValues; ++i)
		{
			bilayerEnergy[i] = 0;
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_bilayer_" << options.minBilayer;
		fileName << "_" << options.maxBilayer << ".dat";
		f_bilayer.open(fileName.str().c_str());
		if(f_bilayer.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	//	If a charge plate potential is selected, we need to allocate memory to store energy values
	if(options.chargePlate)
	{
		int nEnergyValues = 0;
		double chargePlateValue = options.minPlate;
		while(chargePlateValue <= options.maxPlate)
		{
			++nEnergyValues;
			chargePlateValue += options.plateStepSize;
		}
		chargePlateEnergy = new double[nEnergyValues];
		for(int i=0; i<nEnergyValues; ++i)
		{
			chargePlateEnergy[i] = 0;
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_plate_" << options.minPlate;
		fileName << "_" << options.maxPlate << ".dat";
		f_plate.open(fileName.str().c_str());
		if(f_plate.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	//	If a charge plate + options.bilayer potential is selected, we need to allocate memory to store energy values
	//	and we also need to calcualte the background energy
	if(options.biChargePlate)
	{
		// Generate and store array of background values and array to store energy values
		//	Note that the bounds of max and min charge plate and options.bilayer values are fixed here.
		biChargePlateBg = new double[500];
		biChargePlateEnergy = new double[500];
		double sphereRadius = wfData.sphereRadius;
		int counter=0;
		for(int i=0; i<25; ++i)
		{
			for(int j=0; j<20; ++j)
			{
				biChargePlateEnergy[counter]=0;
				double bilayerValue = (double)(i)/2.5;
				double chargePlateValue = (double)(j+1)*2.0;
				//	The energy is in principle given by an infinite sum over image charges. 
				//  As an approximation I use the first 10 terms of the series (i.e the first 10 pairs of image charges)
				*(biChargePlateBg+counter)=-(wfData.nbrUp*wfData.nbrUp+wfData.nbrDown*wfData.nbrDown)*( -(-13*bilayerValue - 27*chargePlateValue - sphereRadius + sqrt(pow(bilayerValue,2) + 6*bilayerValue*chargePlateValue + 9*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(4*pow(bilayerValue,2) + 12*bilayerValue*chargePlateValue + 9*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(16*pow(bilayerValue,2) + 72*bilayerValue*chargePlateValue + 81*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(25*pow(bilayerValue,2) + 100*bilayerValue*chargePlateValue + 100*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(25*pow(bilayerValue,2) + 110*bilayerValue*chargePlateValue + 121*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(36*pow(bilayerValue,2) + 132*bilayerValue*chargePlateValue + 121*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(36*pow(bilayerValue,2) + 144*bilayerValue*chargePlateValue + 144*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(36*pow(bilayerValue,2) + 156*bilayerValue*chargePlateValue + 169*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(49*pow(bilayerValue,2) + 182*bilayerValue*chargePlateValue + 169*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(49*pow(bilayerValue,2) + 196*bilayerValue*chargePlateValue + 196*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(49*pow(bilayerValue,2) + 210*bilayerValue*chargePlateValue + 225*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(64*pow(bilayerValue,2) + 240*bilayerValue*chargePlateValue + 225*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(64*pow(bilayerValue,2) + 256*bilayerValue*chargePlateValue + 256*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(64*pow(bilayerValue,2) + 272*bilayerValue*chargePlateValue + 289*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(81*pow(bilayerValue,2) + 306*bilayerValue*chargePlateValue + 289*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(81*pow(bilayerValue,2) + 324*bilayerValue*chargePlateValue + 324*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(81*pow(bilayerValue,2) + 342*bilayerValue*chargePlateValue + 361*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(100*pow(bilayerValue,2) + 380*bilayerValue*chargePlateValue + 361*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  2*sqrt(100*pow(bilayerValue,2) + 400*bilayerValue*chargePlateValue + 400*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(100*pow(bilayerValue,2) + 420*bilayerValue*chargePlateValue + 441*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  sqrt(121*pow(bilayerValue,2) + 462*bilayerValue*chargePlateValue + 441*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  (2*(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
													  bilayerValue*sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
													  2*chargePlateValue*sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2))))/
												  sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  (4*pow(bilayerValue,2) + 20*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
												   2*bilayerValue*sqrt(4*pow(bilayerValue,2) + 20*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												   5*chargePlateValue*sqrt(4*pow(bilayerValue,2) + 20*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)))/
												  sqrt(4*pow(bilayerValue,2) + 20*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  (9*pow(bilayerValue,2) + 30*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
												   3*bilayerValue*sqrt(9*pow(bilayerValue,2) + 30*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												   5*chargePlateValue*sqrt(9*pow(bilayerValue,2) + 30*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)))/
												  sqrt(9*pow(bilayerValue,2) + 30*bilayerValue*chargePlateValue + 25*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  (2*(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
													  3*bilayerValue*sqrt(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
													  6*chargePlateValue*sqrt(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + pow(sphereRadius,2))))/
												  sqrt(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  (9*pow(bilayerValue,2) + 42*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
												   3*bilayerValue*sqrt(9*pow(bilayerValue,2) + 42*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												   7*chargePlateValue*sqrt(9*pow(bilayerValue,2) + 42*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)))/
												  sqrt(9*pow(bilayerValue,2) + 42*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  (16*pow(bilayerValue,2) + 56*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
												   4*bilayerValue*sqrt(16*pow(bilayerValue,2) + 56*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												   7*chargePlateValue*sqrt(16*pow(bilayerValue,2) + 56*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)))/
												  sqrt(16*pow(bilayerValue,2) + 56*bilayerValue*chargePlateValue + 49*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  (2*(16*pow(bilayerValue,2) + 64*bilayerValue*chargePlateValue + 64*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
													  4*bilayerValue*sqrt(16*pow(bilayerValue,2) + 64*bilayerValue*chargePlateValue + 64*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
													  8*chargePlateValue*sqrt(16*pow(bilayerValue,2) + 64*bilayerValue*chargePlateValue + 64*pow(chargePlateValue,2) + pow(sphereRadius,2))))/
												  sqrt(16*pow(bilayerValue,2) + 64*bilayerValue*chargePlateValue + 64*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
												  (25*pow(bilayerValue,2) + 90*bilayerValue*chargePlateValue + 81*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
												   5*bilayerValue*sqrt(25*pow(bilayerValue,2) + 90*bilayerValue*chargePlateValue + 81*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												   9*chargePlateValue*sqrt(25*pow(bilayerValue,2) + 90*bilayerValue*chargePlateValue + 81*pow(chargePlateValue,2) + pow(sphereRadius,2)))/
												  sqrt(25*pow(bilayerValue,2) + 90*bilayerValue*chargePlateValue + 81*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
												  (2*(4*pow(bilayerValue,2) + 16*pow(chargePlateValue,2) + pow(sphereRadius,2) - 
													  4*chargePlateValue*sqrt(4*pow(bilayerValue,2) + 16*bilayerValue*chargePlateValue + 16*pow(chargePlateValue,2) + pow(sphereRadius,2)) - 
													  2*bilayerValue*(-8*chargePlateValue + sqrt(4*pow(bilayerValue,2) + 16*bilayerValue*chargePlateValue + 16*pow(chargePlateValue,2) + pow(sphereRadius,2)))
													  ))/sqrt(4*pow(bilayerValue,2) + 16*bilayerValue*chargePlateValue + 16*pow(chargePlateValue,2) + pow(sphereRadius,2)))/(2.*pow(sphereRadius,2)));
				*(biChargePlateBg+counter)-=(2*wfData.nbrDown*wfData.nbrUp)*( (40*chargePlateValue + (16*pow(chargePlateValue,2))/
										 sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 (4*pow(sphereRadius,2))/sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(4*pow(bilayerValue,2) + 16*bilayerValue*chargePlateValue + 16*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(16*pow(bilayerValue,2) + 64*bilayerValue*chargePlateValue + 64*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(25*pow(bilayerValue,2) + 100*bilayerValue*chargePlateValue + 100*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(36*pow(bilayerValue,2) + 144*bilayerValue*chargePlateValue + 144*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(49*pow(bilayerValue,2) + 196*bilayerValue*chargePlateValue + 196*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(64*pow(bilayerValue,2) + 256*bilayerValue*chargePlateValue + 256*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(81*pow(bilayerValue,2) + 324*bilayerValue*chargePlateValue + 324*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 4*sqrt(100*pow(bilayerValue,2) + 400*bilayerValue*chargePlateValue + 400*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
										 (4*pow(sphereRadius,2))/sqrt(pow(bilayerValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(9*pow(bilayerValue,2) + 36*bilayerValue*chargePlateValue + 36*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(25*pow(bilayerValue,2) + 100*bilayerValue*chargePlateValue + 100*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(49*pow(bilayerValue,2) + 196*bilayerValue*chargePlateValue + 196*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(81*pow(bilayerValue,2) + 324*bilayerValue*chargePlateValue + 324*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(121*pow(bilayerValue,2) + 484*bilayerValue*chargePlateValue + 484*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(169*pow(bilayerValue,2) + 676*bilayerValue*chargePlateValue + 676*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(225*pow(bilayerValue,2) + 900*bilayerValue*chargePlateValue + 900*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(289*pow(bilayerValue,2) + 1156*bilayerValue*chargePlateValue + 1156*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(361*pow(bilayerValue,2) + 1444*bilayerValue*chargePlateValue + 1444*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) - 
										 2*sqrt(441*pow(bilayerValue,2) + 1764*bilayerValue*chargePlateValue + 1764*pow(chargePlateValue,2) + 4*pow(sphereRadius,2)) + 
										 bilayerValue*(19 + (16*chargePlateValue)/sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2))) + 
										 pow(bilayerValue,2)*(4/sqrt(pow(bilayerValue,2) + 4*bilayerValue*chargePlateValue + 4*pow(chargePlateValue,2) + pow(sphereRadius,2)) + 
														  1/sqrt(pow(bilayerValue,2) + 4*pow(sphereRadius,2))))/(4.*pow(sphereRadius,2)));
				++counter;
			}
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_bichargeplate.dat";
		f_biChargePlate.open(fileName.str().c_str());
		if(f_biChargePlate.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
	//	If 2nd LL type partial-coulomb interaction is selected, we need to read in the background data files
	if(options.partialCoulomb)
	{
		partialEnergy = new double[nbrSamples];
		std::cout << "\tAllocated " << (double)nbrSamples*sizeof(double)/(1024*1024)
		          << " mb to store partial coulomb energy values" << std::endl;
		if(options.partialCoulombType>2)
		{
			if(options.partialCoulombType==3 || options.partialCoulombType==4)
			{
				options.secllType="disc";
			}
			else if(options.partialCoulombType==5 || options.partialCoulombType==6)
			{
				options.secllType="graphene";
			}
			else if(options.partialCoulombType>6 && options.partialCoulombType<=13)
			{
				std::stringstream intType;
				intType.str("");
				intType << "coeff" << options.partialCoulombType-6;
				partialType = intType.str();
			}
			else
			{
				std::stringstream intType;
				intType.str("");
				intType << "coeff" << options.partialCoulombType-13;
				partialType = intType.str();
			}
			fileName.str("");
			fileName << "background_data/second_LL" << partialType << "_background.dat";	
			f_background.open(fileName.str().c_str());
			if(f_background.is_open()==0)
			{
				std::cerr << "\tERROR: cannot open file:" << fileName.str() << std::endl;
				exit(EXIT_FAILURE);
			}
			//	The file is arranged so that the background associated with flux q is stored at location q
			//	Cycle through the list to obtain the correct value
			for(int j=0; j<(int)wfData.flux; ++j) 
			{
				f_background>>this->partialBg;
			}
			//	Background values are storred without the factor of n^2: we need to put this back in!
			if(options.partialCoulombType==3 || options.partialCoulombType==5 || (options.partialCoulombType>6 && options.partialCoulombType<=13))
			{
				this->partialBg*=wfData.nbrUp*wfData.nbrUp+wfData.nbrDown*wfData.nbrDown;
			}
			else
			{
				this->partialBg*=2*wfData.nbrUp*wfData.nbrDown;
			}
			std::cout << "\n\tUsing partial Coulomb 2nd LL Background:" << this->partialBg << std::endl;
			f_background.close();
			//	Set coefficients for effective potential
			if(options.partialCoulombType==3 || options.partialCoulombType==4)
			{
				coeffList[0] = -50.36597363;
				coeffList[1] = 87.38179510;
				coeffList[2] = -56.08455086;
				coeffList[3] = 17.76579124;
				coeffList[4] = -2.971636200;
				coeffList[5] = 0.2513169758;
				coeffList[6] = -0.008434843187;
			}
			else if(options.partialCoulombType==5 || options.partialCoulombType==6)
			{
				//	These values are taken from PRB 75 245440 (2007)
				coeffList[0] = -11.1534;
				coeffList[1] = 24.9078;
				coeffList[2] = -18.6461;
				coeffList[3] = 6.63657;
				coeffList[4] = -1.221097;
				coeffList[5] = 0.1122068;
				coeffList[6] = -0.00404269;
			}
			else if(options.partialCoulombType>6)
			{
				//	Need to read in values from external file
				fileName.str("");
				fileName << "2nd_LL_sphere_potential/" << partialType << "_potential.dat";
				f_coeffs.open(fileName.str().c_str());	
				if(f_coeffs.is_open()==0)
				{
					std::cerr << "\tERROR: cannot open file:" << fileName.str() << std::endl;
					exit(EXIT_FAILURE);
				}
				for(int i=0; i<7; ++i)		
				{
					f_coeffs >> *(i+coeffList);
				}
				f_coeffs.close();
			}
		}
		fileName.str("");
		if(options.dataDirOut!="")	
		{
			fileName << options.dataDirOut << "/";
		}
		fileName << wfData.wfFileName << "_partial_coulomb_";
		fileName << options.partialCoulombType << ".dat";
		f_partialCoulomb.open(fileName.str().c_str());
		if(f_partialCoulomb.is_open()==0)
		{
			std::cerr << "\tCannot open file: " << fileName.str() << std::endl << std::endl;
			exit(EXIT_FAILURE);
		}
		std::cout << "\tWriting data to file\t" << fileName.str().c_str() << std::endl;
	}
}

//!
//! Calculate coulomb energy in the sphere geometry using a chord distance 
//! measure
//!
double AnalysisMethods::CoulombEnergy(
    int n,      //!<    Number of particles
    dcmplx *u,  //!<    Pointer to u data
    dcmplx *v)  //!<    Pointer to v data
{	
	double energy = -n/(2*wfData.sphereRadius);
	for(int i=0; i<n-1; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			energy+=1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
		}
	}
	return energy;
}

//!
//! \brief  Calculate coulomb energy in the disc geometry (see O. Ciftja and 
//! C. Wexler PRB 67 075304 (2003))
//!
double AnalysisMethods::CoulombEnergy(
    int n,      //!<    Number of particles
    dcmplx *z)  //!<    Pointer to z data
{	
	//	Background-background contribution
	double energy=(double)8*n/(3*PI*wfData.sphereRadius);
	for(int i=0; i<n; ++i)
	{
		//	Background-electron contribution
		double dChord = abs(*(z+i))/wfData.sphereRadius;
		if(dChord<5)
		{	
			energy -= *(discBgFunc+(int)floor(DISC_BG_FUNC_DIM*dChord/5.0+0.5))*2/wfData.sphereRadius;
		}
		else
		{
			energy -= 1/(wfData.sphereRadius*dChord);		
		}
		
		for(int j=i+1; j<n; ++j)
		{
			//electron-electron contribution
			energy+=1.0/(n*abs(*(z+i)-*(z+j)));
		}
	}
	return energy;
}

//!
//! This function calculates the average energy from a bilayer type 
//! interaction (bilayer separation bilayerSep)
//!
double AnalysisMethods::BilayerEnergy(
    int nup,            //!<    Number of up particles
    int ndown,          //!<    Number of down particles
    double bilayerValue,//!<    Value of the bilayer parameter
    dcmplx *u,          //!<    Pointer to u data
    dcmplx *v)          //!<    Pointer to v data
{	
	int n = nup+ndown;
	//	Calculate energy using chord distance measure
	//	1. Background term:
	double energy = -nup*nup-ndown*ndown-2*nup*ndown*(-bilayerValue/(2*wfData.sphereRadius)+sqrt(1+bilayerValue*bilayerValue/(4*wfData.sphereRadius*wfData.sphereRadius)));
	energy /= (2*wfData.sphereRadius*n); 
	//	2. Standard Coulomb interaction in the top and bottom layer
	for(int i=0; i<nup-1; ++i)
	{
		for(int j=i+1; j<nup; ++j)
		{
			energy += 1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
		}
	}
	for(int i=nup; i<n; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			energy += 1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
		}
	}
	//	3. Extended potential for the inter-layer Coulomb interaction
	for(int i=0; i<nup; ++i)
	{
		for(int j=nup; j<n; ++j)
		{
			energy+=1.0/(n*sqrt(pow(2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]),2)+bilayerValue*bilayerValue));
		}
	}
	return energy;
}

//!
//! This function calculates the energy due to a charge plate at distance 
//! chargePlateSep above the 2DEG
//!
double AnalysisMethods::ChargePlateEnergy(
    int n,      //!<    Number of particles
    double d,   //!<    Charge plate parameter
    dcmplx *u,  //!<    Pointer to u data
    dcmplx *v)  //!<    Pointer to v data
{	
	//	Calculate energy using chord distance measure
	//	1. Background term:
	double energy = -n*n*(1.0+2.0*d-1.0*sqrt(1+4*d*d))/(2.0*wfData.sphereRadius);
	//  2. Inter-particle term
	for(int i=0; i<n-1; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			energy += 1.0/rChord-1.0/sqrt(rChord*rChord+4*d*d);
		}
	}
	return energy/n;
}

//!
//! This function evaluates the chord distances for all pairs of particles
//! and updates the total count in pairFreqUp, pairFreqDown and  pairFreqDiff
//!
void  AnalysisMethods::PairCorrelation(
    int nup,    //!<    Number of up particles
    int ndown,  //!<    Number of down particles
    dcmplx *u,  //!<    Pointer to u data
    dcmplx *v)  //!<    Pointer to v data
{
	int n = nup+ndown;
	for(int i=0; i<nup; ++i)
	{
		for(int j=i+1; j<nup; ++j)
		{
			double rUp = abs(u[i]*v[j]-u[j]*v[i]);
			//Put these distances into discrete bins between 0 and max_R = Pi*sphereRadius
			this->pairFreqUp[(int)floor(options.nPairBins*rUp+0.5)]++;			
		}
	}
	for(int i=nup; i<n; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			double rDown = abs(u[i]*v[j]-u[j]*v[i]);
			this->pairFreqDown[(int)floor(options.nPairBins*rDown+0.5)]++;
		}
	}
	for(int i=0; i<nup; ++i)
	{
		for(int j=nup; j<n; ++j)
		{
			double rDiff = abs(u[i]*v[j]-u[j]*v[i]);
			this->pairFreqDiff[(int)floor(options.nPairBins*rDiff+0.5)]++;
		}
	}
	return;
}

//!
//! This function evaluates chord distance to the north pole for all 
//! particles in a given configuration, then updates the total count in 
//! densityFreq
//!
void AnalysisMethods::Density(
    int n,      //!<    Number of particles
    dcmplx *u)  //!<    Pointer to u data
{
	for(int i=0; i<n; ++i)
	{
		int binno = (int)floor(options.nDensBins*abs(u[i])+0.5);
		this->densityFreq[binno]++;
	}
	return;
}

//!
//! This function evaluates the electrocstatic energy due to a
//! modified LLL interaction
//!
double AnalysisMethods::SecondLLEnergy(
    int n,
    dcmplx *u,
    dcmplx *v)
{
	//	Calculate energy (per particle) using chord distance measure
	double energy = -this->secondBg;
	for(int i=0; i<n-1; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			energy+=1.0/rChord+exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
			+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
			+rChord*(coeffList[5]+coeffList[6]*rChord))))));
		}
	}
	return energy/n;
}

//!
//! This function evaluates the electrostatic energy due to a modified LLL
//! interactionas a result of finite-thickness effects. There are two types:
//! Quanutm Well ("well") or Heterostructure ("heterojunction").
//!
double AnalysisMethods::FiniteThicknessEnergy(
    int n,      //!<    Number of particles
    double d,   //!<    Value of thickness parameter
    dcmplx *u,  //!<    Pointer to u data
    dcmplx *v)  //!<    Pointer to v data
{
	//	The background energy values for thicknessValue=0.01--10 in steps of 0.01  
	//	are stored in an array. First we determine which array element we need to pick
	int bg = (int)floor(d*100+0.5);
	if(bg<=0 || bg>1000)
	{
		std::cerr << "ERROR: Finite thckness parameter beyond range of stored background values." << std::endl;
		exit(EXIT_FAILURE);
	}
	//  Otherwise set the energy value to the background value
	double energy = -thicknessBg[bg];
	for(int i=0; i<n-1; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			if((rChord/(d))<=1)
			{
				energy += (1/d)*finiteThicknessPotential[(int)floor((rChord/(d))*10000-1)];
			}
			else if((rChord/(d))<=40)
			{
				energy += (1/d)*finiteThicknessPotential[(int)floor((rChord/(d))*2000+8000-1)];	
			}
			else
			{
			    energy += (1/(rChord));
			}
		}
	}
	return energy/n;
}

//////////////////////////////////////////////////////////////////////////////////
//! \brief  This function evaluates the electrostatic energy due to a
//! partial Coulomb	interaction	partialType parameter = 1,3,5 only uses up-up 
//! and down-down pairs partialType parameter = 2,4,6 only uses up-down pairs	
//!	1,2 work in the LLL, 3,4 work in the standard 2nd LL and 5,6 work in 
//! the graphene 2nd LL.
//!
//!	6-13: calculate 2nd LL partial Coulomb terms for up-up and down-down pairs
//!	13-20:	calculate 2nd LL partial Coulomb terms for up-down pairs
//////////////////////////////////////////////////////////////////////////////////																
double AnalysisMethods::PartialCoulombEnergy(
    int nup,        //!<    Number of spin-up particles
    int ndown,      //!<    Number of spin-down particles
    int partialType,//!<    Integer identifier of the type of partial
                    //!<    coulomb interaction
    dcmplx*u,       //!<    Pointer to u data
    dcmplx*v)       //!<    Pointer to v data
{
    double energy = 0.0;
	int n = nup+ndown;
	//	LLL partial Coulomb interaction
	if(partialType==1)
	{
		//	background term:
		energy = -(nup*nup+ndown*ndown)/(2*n*wfData.sphereRadius);
		//	Coulomb energy between up-up and down-down pairs only
		for(int i=0; i<nup-1; ++i)
		{
			for(int j=i+1; j<nup; ++j)
			{
				energy += 1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
			}
		}
		for(int i=nup; i<n-1; ++i)
		{
			for(int j=i+1; j<n; ++j)
			{
				energy += 1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
			}
		}
	}
	else if(partialType==2)
	{
		//	background term:
		energy = -(2*nup*ndown)/(2*n*wfData.sphereRadius);
		//	Coulomb energy between up-down pairs only
		for(int i=0; i<nup; ++i)
		{
			for(int j=nup; j<n; ++j)
			{
				energy += 1.0/(n*2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]));
			}
		}
	}
	//	2nd LL partial Coulomb interaction
	else if(partialType==3 || partialType==5)
	{
		//	background term:
		energy = -this->partialBg;
		//	Coulomb energy between up-up and down-down pairs only
		for(int i=0; i<nup-1; ++i)
		{
			for(int j=i+1; j<nup; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy+=(1/(rChord))+exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		for(int i=nup; i<n; ++i)
		{
			for(int j=i+1; j<n; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy+=(1/(rChord))+exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		energy /= n;
	}
	else if(partialType==4 || partialType==6)
	{
		//	background term:
		energy = -this->partialBg;
		//	Coulomb energy between up-down pairs only
		for(int i=0; i<nup; ++i)
		{
			for(int j=nup; j<n; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy+=(1/(rChord))+exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		energy /= n;
	}
	else if(partialType>6 && partialType<=13)
	{
		//	background term:
		energy = -this->partialBg;
		//	Coulomb energy between up-up and down-down pairs only
		for(int i=0; i<nup-1; ++i)
		{
			for(int j=i+1; j<nup; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy += exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		for(int i=nup; i<n-1; ++i)
		{
			for(int j=i+1; j<n; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy += exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		energy /= n;
	}
	else if(partialType>13)
	{
		//	background term:
		energy = -this->partialBg;
		//	Coulomb energy between up-down pairs only
		for(int i=0; i<nup; ++i)
		{
			for(int j=nup; j<n; ++j)
			{
				double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
				energy += exp(-rChord)*(coeffList[0]+rChord*(coeffList[1]
				+rChord*(coeffList[2]+rChord*(coeffList[3]+rChord*(coeffList[4]
				+rChord*(coeffList[5]+coeffList[6]*rChord))))));
			}
		}
		energy/=n;
	}
	
	return energy;
}

//!
//! This function calculates the energy due to options.bilayer of 
//! layer separation bilayerSep with a charge plate at distance chargePlateSep
//! from each options.bilayer (above or below)
//!																
double AnalysisMethods::BiChargePlateEnergy(
    int nup,            //!<    Number of spin-up particles
    int ndown,          //!<    Number of spin-down particles
    int chargePlateSep, //!<    Charge plate separation parameter
    int bilayerSep,     //!<    Bilayer separation parameter
    int counterBg,      //!<    Counter to point to backgroun energy value
    dcmplx *u,          //!<    Pointer to u data
    dcmplx *v)          //!<    Pointer to v data
{
	int n = nup+ndown;
	double bilayerValue = (double)(bilayerSep)/2.5;
	double chargePlateValue = (double)(chargePlateSep+1)*2.0;
	//	background term:
	double energy = -biChargePlateBg[counterBg];
	//	Interaction in the top and bottom layer
	for(int i=0; i<nup-1; ++i)
	{
		for(int j=i+1; j<nup; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			energy +=  1/sqrt(pow(rChord,2)) - 1/sqrt(36*pow(bilayerValue,2) + 120*bilayerValue*chargePlateValue + 100*pow(chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(64*pow(bilayerValue,2) + 224*bilayerValue*chargePlateValue + 196*pow(chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(4*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 2/sqrt(16*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(36*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(64*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(100*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(144*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(196*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(256*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(324*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(400*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(bilayerValue + 3*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(2*bilayerValue + 3*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(2*bilayerValue + 5*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(3*bilayerValue + 7*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(4*bilayerValue + 9*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(5*bilayerValue + 9*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(5*bilayerValue + 11*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(6*bilayerValue + 11*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(6*bilayerValue + 13*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(7*bilayerValue + 13*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(7*bilayerValue + 15*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(8*bilayerValue + 15*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(8*bilayerValue + 17*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(9*bilayerValue + 17*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(9*bilayerValue + 19*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(10*bilayerValue + 19*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(10*bilayerValue + 21*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(11*bilayerValue + 21*chargePlateValue,2) + pow(rChord,2));
		}
	}
	for(int i=nup; i<n; ++i)
	{
		for(int j=i+1; j<n; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			energy +=  1/sqrt(pow(rChord,2)) - 1/sqrt(36*pow(bilayerValue,2) + 120*bilayerValue*chargePlateValue + 100*pow(chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(64*pow(bilayerValue,2) + 224*bilayerValue*chargePlateValue + 196*pow(chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(4*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 2/sqrt(16*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(36*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(64*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(100*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(144*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(196*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(256*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(324*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(400*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(bilayerValue + 3*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(2*bilayerValue + 3*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(2*bilayerValue + 5*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(3*bilayerValue + 7*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(4*bilayerValue + 9*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(5*bilayerValue + 9*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(5*bilayerValue + 11*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(6*bilayerValue + 11*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(6*bilayerValue + 13*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(7*bilayerValue + 13*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(7*bilayerValue + 15*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(8*bilayerValue + 15*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(8*bilayerValue + 17*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(9*bilayerValue + 17*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(9*bilayerValue + 19*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(10*bilayerValue + 19*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(10*bilayerValue + 21*chargePlateValue,2) + pow(rChord,2)) - 
			1/sqrt(4*pow(11*bilayerValue + 21*chargePlateValue,2) + pow(rChord,2));
		}
	}
	//	Inter-layer interaction
	for(int i=0; i<nup; ++i)
	{
		for(int j=nup; j<n; ++j)
		{
			double rChord = 2*wfData.sphereRadius*abs(u[i]*v[j]-u[j]*v[i]);
			energy +=  1/sqrt(pow(bilayerValue,2) + pow(rChord,2)) + 2/sqrt(4*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(9*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 2/sqrt(16*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(25*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(36*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(49*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(64*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(81*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(100*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(121*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(144*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(169*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(196*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(225*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(256*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(289*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(324*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 
			2/sqrt(361*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) + 
			2/sqrt(400*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2)) - 2/sqrt(441*pow(bilayerValue + 2*chargePlateValue,2) + pow(rChord,2));
		}
	}
	return energy / n;
}

//!
//! This function calculates only the mean and (estimated) 
//! standard deviation energy (without re-sampling procedure)
//!
void AnalysisMethods::MeanEnergy(
    double *energyStore,    //!<    Pointer to list of Coulomb energy values
    double &stDev,          //!<    Address to write standard deviation
    double &meanEnergy)     //!<    Address to write mean
{
	double energy;					//	current value of energy
	double energyav=0,energyav2=0; // <E> and <E^2> of the energy
	for(int k=0; k<nbrSamples; ++k)
	{
		energy = *(energyStore+k);
		energyav += energy;
		energyav2 += energy*energy;
	}
	energyav /= nbrSamples;
	energyav2 /= nbrSamples;
	stDev = sqrt(energyav2-energyav*energyav);
	meanEnergy = energyav;
	return;
}

//!
//! This function tallys up all of the energy values, determines 
//! the average values and then writes data to files. Also normalises 
//! pair correlation and density functions
//!
void AnalysisMethods::FinaliseMethods()
{
	if(options.pairCorrel)
	{
		//	Normalise the pair correaltion functions
		//	up-up is normalised to the value (nup-1)/2*wfData.sphereRadius
		{
		    double normalisation = 0.0;
		    for(int i=0; i<options.nPairBins; ++i)
		    {
			    normalisation += (double)this->pairFreqUp[i]/(i+1);
		    }
		    normalisation /= (double)options.nPairBins*(wfData.nbrUp-1)/(2*wfData.sphereRadius);
		    for(int i=0; i<options.nPairBins; ++i)
		    {
			    this->pairFreqUp[i] /= normalisation;
		    }
		    this->f_pair << "Up-Up" << std::endl;
		    for(int i=0; i<options.nPairBins; ++i)
		    {
			    this->f_pair << this->pairFreqUp[i] << std::endl;
		    }
		}
		if(wfData.nbrDown>0)
		{
			//	down-down is normalised to the value (ndown-1)/2*wfData.sphereRadius
			{
			    double normalisation = 0.0;
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    normalisation += (double)this->pairFreqDown[i]/(i+1);
			    }
			    normalisation /= (double)options.nPairBins*(wfData.nbrDown-1)/(2*wfData.sphereRadius);
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    this->pairFreqDown[i] /= normalisation;
			    }
			    this->f_pair << "Down-Down" << std::endl;
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    this->f_pair << this->pairFreqDown[i] << std::endl;
			    }
			}
			//	up-down is normalised to the value min(nup,ndown)/2*wfData.sphereRadius
			{
			    double normalisation = 0.0;
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    normalisation+=(double)this->pairFreqDiff[i]/(i+1);
			    }
			    normalisation /= (double)options.nPairBins*(std::min(wfData.nbrUp,wfData.nbrDown))/(2*wfData.sphereRadius);
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    this->pairFreqDiff[i] /= normalisation;
			    }
			    this->f_pair << "Up-Down" << std::endl;
			    for(int i=0; i<options.nPairBins; ++i)
			    {
				    this->f_pair << this->pairFreqDiff[i] << std::endl;
			    }
			}
		}
		this->f_pair.close();
	}
	if(options.dens)
	{
		double normalisation = 0.0;
		for(int i=0; i<options.nDensBins; ++i)
		{
			normalisation += (double)this->densityFreq[i]/(i+1);
		}
		normalisation /= options.nDensBins;
		for(int i=0; i<options.nDensBins; ++i)
		{
			this->densityFreq[i] /= normalisation;
		}
		for(int i=0; i<options.nDensBins; ++i)
		{
			this->f_dens << this->densityFreq[i] << std::endl;
		}
		this->f_dens.close();
	}
	if(options.coulomb)
	{
		delete[] coulombEnergy;
		//	Output file containing mean and standard deviation of resampled energy data
		f_energy.precision(15);
		f_energy << coulombMeanEnergy << std::endl << coulombStdev << std::endl;
		f_energy.close();
		std::cout << "\tMean Coulomb energy per electron:\n\n\t" << coulombMeanEnergy 
		          << " +- " << coulombStdev << std::endl << std::endl;
		std::cout << "\tMean Coulomb energy per electron [included 1/sqrt(n/(n-1)) factor]\n\n\t";
		std::cout << coulombMeanEnergy*sqrt((double)(wfData.nbr-1.0)/(wfData.nbr)) << " +- " 
		          << coulombStdev*sqrt((double)(wfData.nbr-1.0)/wfData.nbr) << std::endl << std::endl;
	}
	if(options.secllCoulomb)
	{
		delete[] secondEnergy;
		f_second.precision(15);
		f_second << secondMeanEnergy << std::endl << secondStdev << std::endl;
		f_second.close();
		std::cout << "\tMean second LL (" << options.secllType << " type Potential)";
		std::cout << " energy per electron:\n\n\t" << secondMeanEnergy 
		          << " +- " << secondStdev << std::endl << std::endl;
	}
	if(options.finiteThickness)
	{
		std::cout << "\t\tThickness Parameter\t\tEnergy per electron" << std::endl;	
		std::cout << utilities::cout.HyphenLine() << std::endl;
		double thicknessValue = options.minThickness;
		int counter = 0;
		while(thicknessValue<=options.maxThickness)
		{
			f_thickness.precision(4);
			std::cout.precision(4);
			f_thickness << thicknessValue << "\t";
			std::cout << "\t\t" << thicknessValue << "\t\t\t\t";
			f_thickness.precision(15);
			std::cout.precision(8);
			thicknessEnergy[counter] /= nbrSamples;
			f_thickness << thicknessEnergy[counter] << std::endl;
			std::cout << thicknessEnergy[counter] << std::endl;
			thicknessValue += options.thicknessStepSize;
			++counter;
		}
		std::cout << utilities::cout.HyphenLine() << std::endl;
		f_thickness.close();
		delete[] thicknessEnergy;
	}
	if(options.bilayer)
	{
		std::cout << "\t\tBilayer Separation\t\tEnergy per electron" << std::endl;
		std::cout << utilities::cout.HyphenLine() << std::endl;
		double bilayerValue = options.minBilayer;
		int counter = 0;
		while(bilayerValue<=options.maxBilayer)
		{
			f_bilayer.precision(4);
			std::cout.precision(4);
			f_bilayer << bilayerValue << "\t";
			std::cout << "\t\t" << bilayerValue << "\t\t\t\t";
			f_bilayer.precision(15);
			std::cout.precision(8);
			bilayerEnergy[counter] /= nbrSamples;
			f_bilayer << bilayerEnergy[counter] << std::endl;
			std::cout << bilayerEnergy[counter] << std::endl;
			bilayerValue += options.bilayerStepSize;
			++counter;
		}
		std::cout << utilities::cout.HyphenLine() << std::endl;
		f_bilayer.close();
		delete[] bilayerEnergy;
	}
	if(options.chargePlate)
	{
		std::cout << "\t\tCharge Plate Separation\t\tEnergy per electron" << std::endl;
		std::cout << utilities::cout.HyphenLine() << std::endl;
		double chargePlateValue = options.minPlate;
		int counter = 0;
		while(chargePlateValue<=options.maxPlate)
		{
			f_plate.precision(4);
			std::cout.precision(4);
			f_plate << chargePlateValue << "\t";
			std::cout << "\t\t" << chargePlateValue << "\t\t\t\t";
			f_plate.precision(15);
			std::cout.precision(8);
			chargePlateEnergy[counter] /= nbrSamples;
			f_plate << chargePlateEnergy[counter] << std::endl;
			std::cout << chargePlateEnergy[counter] << std::endl;
			chargePlateValue += options.plateStepSize;
			++counter;
		}
		std::cout << utilities::cout.HyphenLine() << std::endl;
		f_plate.close();
		delete[] chargePlateEnergy;
	}
	if(options.biChargePlate)
	{
		std::cout << "\t\tBilayer Separation\tCharge Plate Separation\t\tEnergy per electron" << std::endl;
		std::cout << utilities::cout.HyphenLine() << std::endl;
		int counter=0;
		for(int d=0; d<25; ++d)
		{
			for(int D=0; D<20; ++D)
			{
				f_biChargePlate.precision(4);
				std::cout.precision(4);
				f_biChargePlate << d << "\t" << D << "\t";
				std::cout << "\t\t" << d << "\t\t\t\t" << D << "\t\t\t\t";
				f_biChargePlate.precision(15);
				std::cout.precision(8);
				biChargePlateEnergy[counter] /= nbrSamples;
				f_biChargePlate << biChargePlateEnergy[counter] << std::endl;
				std::cout << biChargePlateEnergy[counter] << std::endl;
				++counter;
			}
		}
		std::cout << utilities::cout.HyphenLine() << std::endl;
		f_biChargePlate.close();
		delete[] biChargePlateEnergy;
	}
	if(options.partialCoulomb)
	{
		delete[] partialEnergy;
		//	Output file containing mean and standard deviation of resampled energy data
		f_partialCoulomb.precision(15);
		f_partialCoulomb << partialMeanEnergy << std::endl << partialStdev << std::endl;
		f_partialCoulomb.close();
		switch(options.partialCoulombType)
		{
			case 1:
				std::cout << "\tMean LLL Coulomb Potential (up-up pairs only)";
				break;
			case 2:
				std::cout << "\tMean LLL Coulomb Potential (up-down pairs only)";
				break;
			case 3:
				std::cout << "\tMean 2nd LL Coulomb Potential (up-up pairs only)";
				break;
			case 4:
				std::cout << "\tMean 2nd LL Coulomb Potential (up-down pairs only)";
				break;
			case 5:
				std::cout << "\tMean (graphene) 2nd LL Coulomb Potential (up-up pairs only)";
				break;
			case 6:
				std::cout << "\tMean (graphene) 2nd LL Coulomb Potential (up-down pairs only)";
				break;
		}
		if(options.partialCoulombType>6 && options.partialCoulombType<=13)
		{
			std::cout << "\tMean Effective potential (r)^(" << options.partialCoulombType-6 
			          << ")*e^-r (up-up pairs only)"<<std::endl;
		}
			
		if(options.partialCoulombType>13 && options.partialCoulombType<=20)	
		{
			std::cout << "\tMean Effective potential (r)^(" << options.partialCoulombType-13 
			          << ")*e^-r (up-down pairs only)"<<std::endl;
		}
		std::cout << " energy per electron:\n\n\t" << partialMeanEnergy << " +- " << partialStdev << std::endl;
	}
	if(this->pairFreqUp!=0)     delete[] this->pairFreqUp;
	if(this->pairFreqDown!=0)   delete[] this->pairFreqDown;
	if(this->pairFreqDiff!=0)   delete[] this->pairFreqDiff;
	if(this->densityFreq!=0)    delete[] this->densityFreq;
}
#endif
