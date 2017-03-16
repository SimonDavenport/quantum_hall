////////////////////////////////////////////////////////////////////////////////
//!
//!                         \author Simon C. Davenport
//!
//!  \file 
//!		This file contains functions to build file names for input/output
//!     files of the RSES program
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

#ifndef _FILE_NAME_GENERATOR_HPP_INCLUDED_
#define _FILE_NAME_GENERATOR_HPP_INCLUDED_

///////		LIBRARY INCLUSIONS		////////////////////////////////////////////
#include "../fqhe_wave_function_algorithms/fqhe_wave_function.hpp"
#include "../fqhe_wave_function_algorithms/composite_fermion.hpp"

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the wave function label for a given file
//!	
//!	\return Wave function label
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateWaveFunctionName(
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
    if(wfData->type == FQHE::_LAUGHLIN_)
	{
		return "laughlin_";
	}
	else if(wfData->type == FQHE::_COMPOSITE_FERMION_)
	{
	    if(cfData->NEF)
	    {
	        return "nef_";
	    }
	    else
	    {
		    return "pef_";
		}
	}
	else if(wfData->type == FQHE::_BONDERSON_SLINGERLAND_)
	{
		return "bonderson_slingerland_";
	}
	else if(wfData->type == FQHE::_MOORE_READ_)
	{
		return "moore_read_";
	}
	else
	{
	    return "";
	}
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the base file name labelling the output data
//!	
//!	\return base file name
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateBaseFileName(
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
    std::stringstream name;
	name.str("");
	name << parameters->path << "/" << GenerateWaveFunctionName(cfData, wfData);
	name << wfData->fillNumerator << "_" << wfData->fillDenominator << "_n_" << wfData->nbr 
	     << "_size_" << parameters->rows << "x" << parameters->columns;
    return name.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate index file name
//!	
//!	\return indexes file name
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateIndexesFileName(
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
    std::stringstream name;
	name.str("");
    name << parameters->path << "/" << GenerateWaveFunctionName(cfData, wfData);
	name << wfData->fillNumerator << "_" << wfData->fillDenominator << "_n_" << wfData->nbr << "_indexes.dat";
    return name.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate plot file name
//!	
//!	\return plot file name
////////////////////////////////////////////////////////////////////////////////
inline std::string GeneratePlotFileName(
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
    std::stringstream name;
	name.str("");
    name << GenerateBaseFileName(parameters, cfData, wfData) << "_plot.pdf";
    return name.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the file name labelling the output data produced by each
//!	parallel process.
//!	
//!	\return File name
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateNodeFileName(
    const int id,	                            //!<	Node label
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
	std::stringstream name;
	name.str("");
	name << parameters->path << "/" << "process_" << id << "_" << GenerateWaveFunctionName(cfData, wfData);
	name << wfData->fillNumerator << "_" << wfData->fillDenominator << "_n_" << wfData->nbr << "_realcut_" 
	     << parameters->realCut << "_nbr_sectors_" << parameters->nbrSect << "_lza_" << parameters->lzA 
	     << "_size_" << parameters->rows << "x" << parameters->columns << ".tmp";
	return name.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the file name labelling the output data associated with
//!	each angular momentum sector.
//!	
//!	\return File name
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateSectorFileName(
    const int sector,	                        //!<	Sector label
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
	std::stringstream name;
	name.str("");
	name << GenerateBaseFileName(parameters, cfData, wfData) 
	     << "_sector_" << sector+parameters->lzA << "_matrix.dat";
	return name.str();
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Generate the file name labelling the output data associated with
//!	each angular momentum sector.
//!	
//!	\return File name
////////////////////////////////////////////////////////////////////////////////
inline std::string GenerateEigenvaluesFileName(
    const int sector,	                        //!<	Sector label
    const GeneralOptions* parameters,           //!<    Pointer to set of program parameters
    const FQHE::CompositeFermionData* cfData,   //!<    Pointer to composite fermion wave function data
    const FQHE::WaveFunctionData* wfData)       //!<    Pointer to wave function data
{
	std::stringstream name;
	name.str("");
	name << GenerateBaseFileName(parameters, cfData, wfData) 
	     << "_sector_" << sector+parameters->lzA << "_eigenvalues.dat";
	return name.str();
}
#endif
