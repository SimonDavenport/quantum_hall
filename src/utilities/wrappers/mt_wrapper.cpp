////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!                         \date Last Modified: 28/04/2015
//!
//!  \file
//!		Contains wrappers around the mt random number generator class
//!
//!                    Copyright (C) 2015 Simon C Davenport
//!                                                                             
//!     This program is free software: you can redistribute it and/or modify
//!     it under the terms of the GNU General Public License as published by
//!     the Free Software Foundation, either version 3 of the License,
//!     or (at your option) any later version.
//!                                                                             
//!     This program is distributed in the hope that it will be useful, but
//!     WITHOUT ANY WARRANTY; without even the implied warranty of
//!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//!     General Public License for more details.
//!                                                                             
//!     You should have received a copy of the GNU General Public License
//!     along with this program. If not, see <http://www.gnu.org/licenses/>.
//!                                                                             
////////////////////////////////////////////////////////////////////////////////

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "mt_wrapper.hpp"

namespace utilities
{
    //!
    //! Generate a random seed and use that to seed the Mersenne Twister
    //! random number generator
    //!
    void Random::Seed(
        int seedOffset) //!<    Offset in the random seed compared to the default
                        //!      obtained from the job ID or system time
    {
        int seed;
        if(getenv("PBS_JOBID")!=NULL)
        {
            seed = atoi(getenv("PBS_JOBID"));
        }
        else if(getenv("SLURM_JOB_ID")!=NULL)
        {
            seed = atoi(getenv("SLURM_JOB_ID"));
        }
        else
        {
            seed = time(0);
        }		
        seed += seedOffset;
        srand(seed);
        //	Now seed the Mersenne Twister random number generator using
        //	an array of these standard random numbers
        int initArraySize = 100;
        long unsigned int initArray[initArraySize];
        for(int i=0; i<initArraySize; ++i)
        {
	        initArray[i] = rand();
        }
        mt.init_by_array(initArray, initArraySize);
    }
    
    //!
    //! Generate a random double using the Merseene Twister method
    //!
    double Random::GenerateDouble()
    {
        return mt.random();
    }

    //!
    //! Generate a random complex double using the Merseene Twister method
    //!
    dcmplx Random::GenerateComplex()
    {
        double arg = mt.random();
        double phase = mt.random()*2.0*PI;
        return dcmplx(arg*cos(phase), arg*sin(phase));
    }

    //!
    //! Generate a random integer using the Merseene Twister method
    //!
    int Random::GenerateInt()
    {
        return mt.genrand_int31();
    }
}   //  End namespace utilities
