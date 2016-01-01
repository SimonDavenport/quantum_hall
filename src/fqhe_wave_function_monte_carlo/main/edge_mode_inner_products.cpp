
//  mpicpc edge_mode_inner_products.cpp -o edgeModes

#include <iostream>
#include <fstream>
#include <cmath>
#include "../../utilities/general/dcmplx_type_def.hpp"

dcmplx E1(dcmplx* z,const int n,double r);
dcmplx E2(dcmplx* z,const int n,double r);

int main()
{
    //  Read in data file
    
    std::ifstream f_r,f_theta;
    
    long int dim = 400000;
    int n = 160;
    double p = 3;
    double discR = sqrt(2.0*n*p);
    
    double* rData = new double[dim*n];
    double* tData = new double[dim*n];
    dcmplx* zData = new dcmplx[dim*n];
    
    f_r.open("/scratch/scd51/physics_programs/fqhe_monte_carlo/fermions_disc_n_160_flux_477_run_1_r.dat",std::ios::binary);

    if(!f_r.is_open())   std::cout<<" ERROR 1"<<std::endl;
    
    f_r.read(reinterpret_cast<char*>(rData),dim*n*sizeof(double));

    f_r.close();

    f_theta.open("/scratch/scd51/physics_programs/fqhe_monte_carlo/fermions_disc_n_160_flux_477_run_1_theta.dat",std::ios::binary);
   
    if(!f_theta.is_open())   std::cout<<" ERROR 2"<<std::endl;
    
    f_theta.read(reinterpret_cast<char*>(tData),dim*n*sizeof(double));

    f_theta.close();
    
    //  Convert to z co-ords
    
    for(long int i=0;i<n*dim;i++)
    {
        zData[i] = rData[i]*dcmplx(sin(tData[i]),cos(tData[i]));
    }
    
    delete[] rData;
    delete[] tData;
   
    //  Evaluate function
    
    dcmplx mean  = 0.0;
    dcmplx mean1 = 0.0;
    dcmplx mean2 = 0.0;
    dcmplx mean3 = 0.0;
    dcmplx mean4 = 0.0;
    dcmplx mean5 = 0.0;
    dcmplx* p_z  = zData;
    
    for(long int i=0;i<dim;i++,p_z += n)
    {
        dcmplx e1 = E1(p_z,n,discR);
        dcmplx e2 = E2(p_z,n,discR);
    
        mean  += e1*std::conj(e1);
        mean1 += e2*std::conj(e2);
        mean2 += std::conj(e1)*std::conj(e1)*e2;
        mean3 += std::conj(e1)*std::conj(e2)*e1*e2;
        mean4 += std::conj(e2)*std::conj(e2)*e2*e2;
        mean5 += std::conj(e1)*std::conj(e1)*e1*e1;
    }
    
    mean  /= dim; 
    mean1 /= dim;
    mean2 /= dim;
    mean3 /= dim;
    mean4 /= dim;
    mean5 /= dim;
    
    mean  *= p;
    mean1 *= p;
    mean2 *= p*sqrt(p);
    mean3 *= p*p;
    mean4 *= p*p;
    mean5 *= p*p;
    
    delete[] zData;
    
    std::cout<<"<J_1 J_-1> = "<<mean<<std::endl;
    std::cout<<"<J_2 J_-2> = "<<mean1<<std::endl;
    std::cout<<"<J_1 J_1 J_-2> = "<<mean2<<std::endl;
    std::cout<<"<J_2 J_1 J_-1 J_-2> = "<<mean3<<std::endl;
    std::cout<<"<J_2 J_2 J_-2 J_-2> = "<<mean4<<std::endl;
    std::cout<<"<J_1 J_1 J_-1 J_-1> = "<<mean5<<std::endl;

    return 0;
}

dcmplx E1(dcmplx* z,const int n,const double r)
{
    dcmplx temp = 0.0;
    
    for(int i=0;i<n;i++)
    {
        temp += z[i];
    }

    return temp/r;
}

dcmplx E2(dcmplx* z,const int n,const double r)
{
    dcmplx temp = 0.0;
    
    for(int i=0;i<n;i++)
    {
        temp += z[i]*z[i];
    }

    return temp/(r*r);
}
