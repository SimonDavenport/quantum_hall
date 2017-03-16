////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport     
//!                                                                             
//!	 \file
//!     This file declares the program options associated with calculation of 
//!     the RSES
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

#ifndef _RSES_OPTIONS_HPP_INCLUDED_
#define _RSES_OPTIONS_HPP_INCLUDED_

#include <boost/program_options.hpp>

namespace myOptions
{
    namespace po = boost::program_options;
    inline po::options_description GetEntanglementSpectrumOptions()
    {
        po::options_description rsesOpt("RSES options");
        rsesOpt.add_options()
        ("nbr-n,n", po::value<int>()->default_value(12),
        "Number of particles\n")
        ("nbr-a,a", po::value<int>()->default_value(6),
        "Number of particles in the A cut\n")
        ("sector,s", po::value<int>()->default_value(0),
        "Specify the sector to obtain OR\n")
        ("multi-sector,m", po::value<bool>()->default_value(false),
        "Set to true to use a range of sectors rather than just one, or alternatively specify a maximum and minimum sector to diagonalize.\n")
        ("max-sector", po::value<int>()->default_value(0),
        "Specify the maximum sector to obtain\n")
        ("min-sector", po::value<int>()->default_value(0),
        "Specify the minimum sector to obtain\n")
        ("nbr-ll,l", po::value<int>()->default_value(1),
        "specify the number of Landau levels (maximum 3 levels currently programmed)\n");
        return rsesOpt;
    };
}   //  End myOptions namespace
#endif
