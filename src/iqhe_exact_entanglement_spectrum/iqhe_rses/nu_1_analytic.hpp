////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!                        \author Simon C. Davenport                        
//!                                                                             
//!	 \file
//!     This file implements the analytic form of the nu=1 entanglement
//!     energy function          
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

#ifndef _NU_1_ANALYTIC_HPP_INCLUDED_
#define _NU_1_ANALYTIC_HPP_INCLUDED_

///////     LIBRARY INCLUSIONS     /////////////////////////////////////////////
#include <gsl/gsl_sf_gamma.h> 
#include <gsl/gsl_sf_hyperg.h>
inline double AnalyticEntanglementEnergy(
    double q,
    double m)
{
    double gslFactorialResult;
    double gslFactorialResult2;
    double gslFactorialResult3;
    double gslHyperResult;
    gslFactorialResult = gsl_sf_gamma(2.0*q+1.0);
    std::cout << "\tgslFactorialResult = " << gslFactorialResult << std::endl;
    gslFactorialResult2 = gsl_sf_gamma(q-m+1.0);
    std::cout << "\tgslFactorialResult2 = " << gslFactorialResult2 << std::endl;
    gslFactorialResult3 = gsl_sf_gamma(q+m+2.0);
    std::cout << "\tgslFactorialResult3 = " << gslFactorialResult3 << std::endl;
    gslHyperResult = gsl_sf_hyperg_2F1(1.0, m-q, 2.0+q+m, -1.0);
    std::cout << "\tgslHyperResult = " << gslHyperResult << std::endl;
    double alpha = (1.0+2.0*q)*gslFactorialResult;
    alpha /= std::pow(2.0, 1+2.0*q);
    alpha *= gslHyperResult;
    alpha /= gslFactorialResult2;
    alpha /= gslFactorialResult3;
    alpha = fabs(alpha);
    std::cout << "\talpha = " << alpha << std::endl;
    if(alpha>1.0)
    {
        gslFactorialResult = gsl_sf_gamma(2.0*q+1.0);
        std::cout << "\tgslFactorialResult = " << gslFactorialResult << std::endl;
        gslFactorialResult2 = gsl_sf_gamma(m+1.0);
        std::cout << "\tgslFactorialResult2 = " << gslFactorialResult2 << std::endl;
        gslFactorialResult3 = gsl_sf_gamma(2*q-m+2.0);
        std::cout << "\tgslFactorialResult3 = " << gslFactorialResult3 << std::endl;
        gslHyperResult = gsl_sf_hyperg_2F1(1.0, -m, 2.0+2*q-m, -1.0);
        std::cout << "\tgslHyperResult = " << gslHyperResult << std::endl;
        double alpha = (1.0+2.0*q)*gslFactorialResult;
        alpha /= std::pow(2.0, 1+2.0*q);
        alpha *= gslHyperResult;
        alpha /= gslFactorialResult2;
        alpha /= gslFactorialResult3;
        alpha = fabs(alpha);
        std::cout << "\t 1 - alpha = " << alpha << std::endl;
        #if _ENABLE_HIGH_PRECISION_
        const int P = 512;    //  Set precision
        mpfr_t temp;
        mpfr_t temp1;
        mpfr_init2(temp, P);
        mpfr_init2(temp1, P);
        mpfr_set_d(temp, alpha, MPFR_RNDN);
        mpfr_ui_sub(temp1, 1, temp, MPFR_RNDN);
        mpfr_div(temp1, temp1, temp, MPFR_RNDN);
        mpfr_log(temp1, temp1, MPFR_RNDN);
        double returnVal = mpfr_get_d(temp1, MPFR_RNDN);
        mpfr_clear(temp);
        mpfr_clear(temp1);
        return returnVal;
        #else
        return log((1.0-alpha)/alpha);
        #endif
    }
    else
    {
        #if _ENABLE_HIGH_PRECISION_
        const int P = 512;    //  Set precision  
        mpfr_t temp;
        mpfr_t temp1;
        mpfr_init2(temp, P);
        mpfr_init2(temp1, P);
        mpfr_set_d(temp, alpha, MPFR_RNDN);
        mpfr_ui_sub(temp1, 1, temp, MPFR_RNDN);
        mpfr_div(temp1, temp, temp1, MPFR_RNDN);
        mpfr_log(temp1, temp1, MPFR_RNDN);
        double returnVal = mpfr_get_d(temp1, MPFR_RNDN);
        mpfr_clear(temp);
        mpfr_clear(temp1);
        return returnVal;
        #else
        return log(alpha/(1.0-alpha));
        #endif
    }
}
#endif
