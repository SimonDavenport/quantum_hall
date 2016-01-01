bin_PROGRAMS = \
fqhe_entanglement_energy \
fqhe_entanglement_hamiltonian \
fqhe_entanglement_wave_function \
fqhe_monte_carlo \
analyse_monte_carlo \
edge_mode_inner_products \
fqhe_entanglement_spectrum \
iqhe_entanglement_spectrum

##  FQHE entanglement energy model

fqhe_entanglement_energy_SOURCES = \
fqhe_entanglement_energy_method/main/fqhe_entanglement_energy.cpp \
occupation_basis/occupation_basis.cpp \
fqhe_entanglement_energy_method/entanglement_energy_model/entanglement_energy_model.cpp \
utilities/mathematics/integer_partitions.cpp

fqhe_entanglement_energy_LINK = $(FC) $(LDFLAGS) -o fqhe_entanglement_energy

##  FQHE entanglement Hamiltonian

fqhe_entanglement_hamiltonian_SOURCES = \
fqhe_entanglement_hamiltonian_method/main/fqhe_entanglement_hamiltonian.cpp \
fqhe_entanglement_hamiltonian_method/entanglement_hamiltonian/list_of_terms.cpp \
fqhe_entanglement_hamiltonian_method/entanglement_hamiltonian/hilbert.cpp \
fqhe_entanglement_hamiltonian_method/entanglement_hamiltonian/term.cpp \
occupation_basis/occupation_basis.cpp \
utilities/mathematics/integer_partitions.cpp

fqhe_entanglement_hamiltonian_LINK = $(FC) $(LDFLAGS) -o fqhe_entanglement_hamiltonian 

##  FQHE entanglement wave function method

fqhe_entanglement_wave_function_SOURCES = \
fqhe_entanglement_wave_function_method/algebra/algebra.cpp \
fqhe_entanglement_wave_function_method/algebra/beta.cpp \
fqhe_entanglement_wave_function_method/algebra/deter.cpp \
fqhe_entanglement_wave_function_method/algebra/gamma.cpp \
fqhe_entanglement_wave_function_method/gnuplot/gnuplot_i.cpp \
fqhe_entanglement_wave_function_method/hilbert_ewf/hilbert_ewf.cpp \
fqhe_entanglement_wave_function_method/main/entanglement_wave_function.cpp \
fqhe_entanglement_wave_function_method/metropolis/metropolis.cpp \
fqhe_entanglement_wave_function_method/wave_functions/wave_functions.cpp \
utilities/mathematics/mt.cpp

fqhe_entanglement_wave_function_LINK = $(FC) $(LDFLAGS) -o fqhe_entanglement_wave_function

##  FQHE wave function monte carlo

fqhe_monte_carlo_SOURCES = \
fqhe_wave_function_monte_carlo/main/fqhe_monte_carlo.cpp \
fqhe_wave_function_monte_carlo/metropolis/fqhe_metropolis.cpp \
fqhe_wave_function_algorithms/fqhe_wave_function.cpp \
fqhe_wave_function_algorithms/composite_fermion.cpp \
fqhe_wave_function_algorithms/laughlin.cpp \
fqhe_wave_function_algorithms/moore_read.cpp \
fqhe_wave_function_algorithms/nass.cpp \
fqhe_wave_function_algorithms/parafermion.cpp \
utilities/mathematics/mt.cpp \
utilities/mathematics/binomial_table.cpp \
utilities/mathematics/theta_functions.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/wrappers/mpi_wrapper.cpp

analyse_monte_carlo_SOURCES = \
fqhe_wave_function_monte_carlo/main/analyse_monte_carlo.cpp \
fqhe_wave_function_monte_carlo/analysis/fqhe_monte_carlo_analysis.cpp  \
fqhe_wave_function_algorithms/fqhe_wave_function.cpp \
fqhe_wave_function_algorithms/composite_fermion.cpp \
fqhe_wave_function_algorithms/laughlin.cpp \
fqhe_wave_function_algorithms/moore_read.cpp \
fqhe_wave_function_algorithms/nass.cpp \
fqhe_wave_function_algorithms/parafermion.cpp \
utilities/mathematics/mt.cpp \
utilities/mathematics/binomial_table.cpp \
utilities/mathematics/theta_functions.cpp \
utilities/mathematics/binary_number_tools.cpp \
utilities/mathematics/murmur_hash2.cpp \
utilities/mathematics/murmur_hash3.cpp \
utilities/wrappers/mpi_wrapper.cpp

analyse_monte_carlo_LINK = $(FC) $(LDFLAGS) -o analyse_monte_carlo 

edge_mode_inner_products_SOURCES = \
fqhe_wave_function_monte_carlo/main/edge_mode_inner_products.cpp

edge_mode_inner_products_LINK = $(FC) $(LDFLAGS) -o edge_mode_inner_products

##  FQHE pseudo reduced density matrix method

fqhe_entanglement_spectrum_SOURCES = \
fqhe_pseudo_reduced_density_matrix_method/main/fqhe_entanglement_spectrum.cpp \
utilities/mathematics/mt.cpp \
utilities/mathematics/theta_functions.cpp \
fqhe_wave_function_algorithms/fqhe_wave_function.cpp \
fqhe_wave_function_algorithms/composite_fermion.cpp \
fqhe_wave_function_algorithms/laughlin.cpp \
fqhe_wave_function_algorithms/moore_read.cpp \
utilities/wrappers/mpi_wrapper.cpp

fqhe_entanglement_spectrum_LINK = $(FC) $(LDFLAGS) -o fqhe_entanglement_spectrum 

##  IQHE exact entanglement spectrum

iqhe_entanglement_spectrum_SOURCES = \
iqhe_exact_entanglement_spectrum/main/iqhe_entanglement_spectrum.cpp \
iqhe_exact_entanglement_spectrum/iqhe_rses/iqhe_rses.cpp \
utilities/mathematics/binomial_table.cpp \
occupation_basis/occupation_basis.cpp \
utilities/mathematics/integer_partitions.cpp

iqhe_entanglement_spectrum_LINK = $(FC) $(LDFLAGS) -o iqhe_entanglement_spectrum

