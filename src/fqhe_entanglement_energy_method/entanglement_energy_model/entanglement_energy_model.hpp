////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport        
//!                                                                             
//!	 \file
//!     This file implements an entanglement energy type fitting function
//!     as an alternative to the CFT form of the entanglement Hamiltonian           
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

#ifndef _ENTANGLEMENT_ENERGY_MODEL_HPP_INCLUDED_
#define _ENTANGLEMENT_ENERGY_MODEL_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../occupation_basis/occupation_basis.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/program_options.hpp>
#include "../../utilities/algorithms/quick_sort.hpp"
#include <algorithm>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif
#if _ENABLE_HIGH_PRECISION_
#include <mpfr.h>
#endif

////////////////////////////////////////////////////////////////////////////////
//! \brief This class implements an ansatz expression for the entanglement
//! energy assigned to FQHE RSES. The form is based on some small perturbations
//! to the integer quantum Hall form of the entanglement energy, implemented
//! in iqheRses.h
////////////////////////////////////////////////////////////////////////////////
class EntanglementEnergyModel
{
    private:
    void GenerateEntanglementEnergies(
        double* energyLevels, const double q, double& normalization, const int nbrLevels,
        const int nbrStates, std::vector<double>& modelParameters) const;
    void ParametersFromFile(
        boost::program_options::variables_map* optionList, 
        std::vector<double>& modelParameters, std::vector<double>& qModelParameters) const;
    void SpectrumToFile(
        boost::program_options::variables_map* optionList, double* spectrum,
        std::vector<int>& spectrumLabels, const int spectrumDim, const int sector) const;
    void ApplyPerturbation(
        std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,
        double* energyLevels, double* spectrum, const int nbrStates, 
        std::vector<double>& modelParameters) const;
    public:
    void BuildEntanglementEnergyModel(
        boost::program_options::variables_map* optionList, const int lzaMin, 
        const int lzaMax) const;
};
#endif
