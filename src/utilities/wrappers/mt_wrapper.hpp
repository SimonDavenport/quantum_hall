////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport 
//!
//!  \file
//!		Contains wrappers around the mt random number generator class
//!
//!                    Copyright (C) Simon C Davenport
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

#ifndef _MT_WRAPPER_HPP_INCLUDED_
#define _MT_WRAPPER_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../mathematics/mt.hpp"
#include "../general/template_tools.hpp"
#include "../general/dcmplx_type_def.hpp"
#include "../general/pi_const_def.hpp"

namespace utilities
{
    class Random
    {
        private:
        MersenneTwister mt;
        void Seed(int seedOffset);

        public:
        
        Random()
        {
            this->Seed(0);
        }
        
        Random(int seed)
        {
            this->Seed(seed);
        }

        double GenerateDouble();
        dcmplx GenerateComplex();
        int GenerateInt();
    };
}  //  End namespace utilities
#endif
