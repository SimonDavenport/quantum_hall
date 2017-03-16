A collection of programs for the analysis of fractional quanutm Hall wave functions. 

Applications of the code are described in more detail in these research papers:

https://arxiv.org/abs/1203.0004  (Monte Carlo computations of energies for multicomponent Composite Fermion states)

https://arxiv.org/abs/1305.5834  (A Monte Carlo apporach to approximate computation of the entanglement spectrum of
quantum Hall states, expecially useful for application to Composite Fermion states)

https://arxiv.org/abs/1506.06741 (An alternative approach to computation of entanglement spectra)

INSTALL: autoconf automake chmod +x configure

./configure --enable-blas --enable-lapack --enable-speed-optimization

DEPENDENCIES: boost program options

Google sparse-hash (if using speed optimization)
