////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!		This file defines a class to calculate the numerically exact real
//!     space entanglement spectrum for integer quantum Hall states
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

#ifndef _IQHE_RSES_HPP_INCLUDED_
#define _IQHE_RSES_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include "../../occupation_basis/occupation_basis.hpp"
#include "../../utilities/wrappers/lapack_wrapper.hpp"
#include "../../utilities/general/pi_const_def.hpp"
#include "../../utilities/mathematics/binomials.hpp"
#include "../../utilities/general/cout_tools.hpp"
#include <boost/program_options.hpp>
#include <gsl/gsl_integration.h>
#include <vector>
#include <algorithm>
#if _DEBUG_
#include "../../utilities/general/debug.hpp"
#endif
//////      Predeclare integrand function compatible with gsl_integrate
double Integrand(double x, void* params);

//!
//! Define a class to calculate the numerically exact real space 
//! entanglement spectrum of integer quantum Hall states
//!
class IqheRses
{
    private:
    friend class EntanglementEnergy;
    double EvaluateMonopoleHarmonic(const double q, const int l, const double m, const double x) const;
    double MonopoleHarmonicNorm(const double q, const int l, const double m) const;
    friend double Integrand(double x, void* params);
    void GenerateOverlapMatrix(double* overlapMatrixBuffer, const double q,
                               const int nbrLevels, const int nbrStates) const;
    void GenerateNewLevels(double* overlapMatrixBuffer, double* energyLevels, double& normalization,
                           const int nbrLevels, const int nbrStates) const;
    #if _DEBUG_
    void CheckAngularMomentum(std::vector<std::vector<std::bitset<_MAX_NUMBER_> > >& occupationsList,
                              const double monopoleStrength, const int nbrA, const int nbrLevels, 
                              const int sector) const;
    #endif
    std::string GenFileName(boost::program_options::variables_map* optionList, const int sector) const;
    void SpectrumToFile(boost::program_options::variables_map* optionList, double* spectrum,
                        const int spectrumDim, const int sector) const;
    public:
    IqheRses();
    ~IqheRses();
    void GenerateIqheSpectrum(boost::program_options::variables_map* optionList,
                              const int lzaMin, const int lzaMax);
};
#endif
