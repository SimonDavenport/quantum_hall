////////////////////////////////////////////////////////////////////////////////
//!                                                                             
//!         \author Simon C Davenport and Ivan D Rodriguez
//!                                                                             
//!                      \date Last Modified: 03/06/2015
//!                                                                             
//!	 \mainpage
//!     This suite of programs contains a number of codes used to calculate 
//!     the real space entanglement spectrum (RSES) of model integer and
//!     fractional quantum Hall wave functions.
//!     
//!     INSTALLATION:
//!
//!     To install this code you will first need to download and install 
//!     the following libraries:
//!
//!     - https://gmplib.org/ (high precision library)
//!     - http://www.mpfr.org/ (high precision library)
//!     - http://www.multiprecision.org/ (high precision library)
//!     - http://www.boost.org/ (we only need to the program options)
//!     - https://www.gnu.org/software/gsl/ (for numerical quadrature)
//!     - http://fftw.org/ (for fast Fourier transforms)
//!
//!     Once this is completed, you can install our codes using a 
//!     standard GNU build system.
//!     
//!     - to see the configuration options run ./configure --help 
//!     - when you are happy run ./configure [options]
//!     (note that the installation directory is set with --prefix)
//!     - modify the options/check library installations until
//!     all the configure tests are passed. 
//!     - run make -j to compile, and make -j install to move
//!     the compiled codes to the installation directory. Run make
//!     clean to tidy up the build (and leave the installation 
//!     unchanged). Run make uninstall to remove the installation. 
//!
//!     Now that the code is installed, several methods are available 
//!     to you for RSES calculations:
//!
//!     - fqhe entanglement energy method: 
//!     See entanglement_energy_model.cpp. Treats the RSES for 
//!     fractional quantum Hall states in a non-interacting approximation 
//!     in terms of a set of single-particle entanglement energy levels.
//!     Author: Simon C Davenport
//!
//!     - fqhe entanglement hamiltonian method: 
//!     See fqhe_entanglement_hamiltonian.cpp. Treats the RSES for integer
//!     or fractional quantum Hall states in the formalism of Dubail-Read
//!     and Rezayi's paper Phys. Rev. B 86, 245310 (2012). The formalism
//!     has also been extended in an attempt to treat composite fermion
//!     fractional quantum Hall states.
//!     Author: Simon C Davenport
//!
//!     - fqhe entanglement wave function method: 
//!     See entanglement_wave_function.cpp. A specialist method
//!     based on the "entanglement wave functions" proposed by Rodriguez 
//!     et al. in Phys. Rev. B 88, 155307. This method can be used to 
//!     calculate the RSES of the Laughlin and positive effective field
//!     Jain states very very large system sizes. 
//!     Author: Ivan D Rodriguez
//!
//!     - fqhe pseudo reduced density matrix method: 
//!     See fqhe_entanglement_spectrum.cpp. A more general method
//!     to calculate the approximate RSES associated with quantum Hall 
//!     states with known trial wave functions. Based on the 
//!     observations of Rodriguez et al. in Phys. Rev. Lett. 108, 256806.
//!     Authors: Ivan D Rodriguez and Simon C Davenport
//!
//!     - iqhe exact entanglement spectrum: 
//!     See iqhe_entanglement_spectrum.cpp. Calculates the RSES of the
//!     integer quantum Hall states according to the analytical 
//!     result (but doing the requisite integrals numerically).
//!     Author: Simon C Davenport
//!
//!     There are also several other sets of algorithms and utilities
//!     bundled with the code:
//!
//!     - fqhe wave function algorithms: 
//!     See the FQHE namespace. A set of fairly state-of-the-art
//!     algorithms to evaluate a large number of different types of 
//!     fractional quantum Hall ground state wave function. Mostly
//!     done in the sphere geometry case, but disk and torus geometry
//!     is sometimes available. 
//!     Author: Simon C Davenport
//!
//!     - fqhe wave function monte carlo: 
//!     See fqhe_monte_carlo.cpp. An older code containing 
//!     some basic routines to evaluate the ground state energy and
//!     ground state correlation functions of fractional quantum 
//!     Hall states using monte carlo sampling. 
//!     Author: Simon C Davenport
//!
//!     - occupation basis: 
//!     See occupation_basis.cpp. Some shared utilities for generating the 
//!     angular momentum basis used by the other RSES calculations
//!     in the sphere geometry case. 
//!     Author: Simon C Davenport
//!
//!     - program_options: 
//!     See e.g. general_options.h Some shared utilities for program option
//!     handling.
//!     Author: Simon C Davenport
//!
//!     - python codes:
//!     A set of additional python codes that can be used
//!     to analyse the RSES data generated by the other programs. In 
//!     particular "fit_entanglement_spectrum.py" can be used to 
//!     fit one set of RSES data to one of the models. 
//!     Author: Simon C Davenport
//!
//!     - utilities: 
//!     See the utilities namespace. A set of general purpose utility  
//!     functions required for low level mathematical functions e.g.
//!     linear algebra, integer partitions and some C++ hacks. 
//!                                                                                                         
////////////////////////////////////////////////////////////////////////////////

