#! /usr/bin/env python2.7
################################################################################
##                                                                            ##
##						  Author: Simon C. Davenport						  ##
##                                                                            ##
##							Last Modified: 10/03/2015						  ##
##                                                                            ##
##                                                                            ##
##		Perform fitting of the real space entanglement spectrum of a          ##
##      specified model. The fitting function can be either an                ##
##		entanglement Hamiltonian model of the Dubail et al. form, a           ##
##      a generalization of their entanglement Hamiltonian to multiple        ##
##      Landau levels, or alternatively a single particle entanglement        ##
##      energy model.                                                         ##
##                                                                            ##
##					Copyright (C) 2013-2015 Simon C Davenport			      ##
##                                                                            ##
################################################################################

################################################################################
################################################################################
##  Import modules

##  Set module import paths for utility libraries
import sys
sys.path.append(r'/scratch/scd51/Dropbox/physics_programs/python')
sys.path.append(r'/rscratch/scd51/python')

##  Import utility libraries
from utilities import plot_utilities
from utilities import curve_fitting

##  Import RSES data processing libraries from rses_data_processing.py
from rses_data_processing import RsesData
from rses_data_processing import EwfRses
from rses_data_processing import RdmRses
from rses_data_processing import ExactRses

##  Import fitting model processing libraries from fitting_model_processing.py
from fitting_model_processing import FittingModel
from fitting_model_processing import SingleParticleEnergy
from fitting_model_processing import DRR431Hamiltonian
from fitting_model_processing import ExtendedDRRHamiltonian

##  Import additional standard python libraries
import math
import numpy as np
import scipy
from scipy.optimize import minimize  ##version 0.11 and up
##from scipy.optimize import fmin    ##version 0.10
import os
import argparse

#	these functions enable zoomed insets
#	see http://matplotlib.org/mpl_toolkits/axes_grid/users/overview.html
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

################################################################################
################################################################################
##  Define a class to import and contain program options
class ProgramOptions:
    
    ############################################################################
    ##  Define a full list of command line arguments to be parsed
    ##
    def __init__(self):

        parser = argparse.ArgumentParser(description='This program fits a specified entanglement Hamiltonian or single particle energy spectrum to a given RSES')
        
        ##  Set running options
        
        parser.add_argument('--generate-minima-file',action='store_true',default=False,help='Specify this option to read in exact RSES data and generate files contianing the minimum entanglement energy values')
        parser.add_argument('--get-counting',action='store_true',default=False,help='Specify this option to count the number of RSES states below a certian entanglement energy, then quit')
        parser.add_argument('--fit-rses',action='store_true',default=False,help='Specify this option to fit the specified entanglement Hamiltonian to a given set of RSES data')
        parser.add_argument('--make-plots',action='store_true',default=False,help='Specify this option to generate plots of the RSES data alongside the fitting data')
        parser.add_argument('--plot-rses',action='store_true',default=False,help='Specify this option to plot the numerical RSES only')
        parser.add_argument('--plot-rses-convergence',action='store_true',default=False,help='Specify this option to plot the convergence of the RSES eigenvalues')
        parser.add_argument('--plot-fitting-convergence',action='store_true',default=False,help='Specify this option to plot the convergence of the fitting model parameters as a function of the sector fit')
        parser.add_argument('--plot-thermodynamic-convergence',action='store_true',default=False,help='Specify this option to plot the convergence of the fitting model parameters as a function of the system size')
        parser.add_argument('--calculate-entanglement-entropy',action='store_true',default=False,help='Specify this option to recalculate the sum of exp(-E_i) and the Von Neumann, Renyi and minimum entanglement entropy for a given model')
        
        ##  Set system parameters
        
        parser.add_argument('--wf-type',default='laughlin',help='Set the wave function type to look for RSES data e.g. laughlin,pef,iqhe')
        parser.add_argument('--nbr',default=12,type=int,help='Set the number of particles in the spectrum ')
        parser.add_argument('--nbr-a',default=6,type=int,help='Set the overall particle cut by specifying the number of particles in the A subsystem.')
        parser.add_argument('--extra-nbr-a',default=0,type=int,help='Set to EXTRA_NBR_A>0 to additionally fit the RSES with all particle cuts within nbr-a +/-EXTRA_NBR_A')
        parser.add_argument('--nbr-ll',default=1,type=int,help='Number of CF LLs in the construction')
        parser.add_argument('--jastrow',default=1,type=int,help='Power in Jastrow factor - NOT including Slater determinant part (set to zero for IQHE)')
        parser.add_argument('--fitting-model',default='DRR431',help='Specify fitting model name e.g. DRR431, extended_DRR, or entanglement_energy')
        parser.add_argument('--statistics',default="fermions",help='Specify fermions or bosons (value will be set consistently with --jastrow)')
        parser.add_argument('--fit-offset',action='store_true',default=False,help='Include offset parameters depending on the value of NA in the fitting procedure')
        parser.add_argument('--use-perturbation',action='store_true',default=False,help='Specify this option to include a perturbation to the single particle energy levels.\n')

        ##  RSES data options
        
        parser.add_argument('--rses-method',default='ewf',help='Set method used to calculate rses: "ewf" (entanglement wave function overlaps),"rdm" (approximated reduced density matrix) or "exact" (the IQHE RSES can be calculated exactly). This option sets the format for data files to be read in, and performs no calcualtions.')
        parser.add_argument('--rses-steps',default=6400,type=int,help='For "ewf" type RSES: sets the number of Monte Carlo steps in the RSES calcualtion (in 1000s)')
        parser.add_argument('--rses-size',default=4096,type=int,help='For "rdm" type RSES: dimension RSES_SIZE of pseudo reduced density matrix for the RSES data (pseudo reduced density matrix dimensions will be RSES_SIZE by 2*RSES_SIZE)')
        parser.add_argument('--rses-skip',default=1,type=int,help='For "rdm" type RSES: skip the first RSES_SKIP sectors files whilst reading RSES data.')
        parser.add_argument('--rses-norm-zero',action='store_true',default=False,help='Normalize the RSES data to zero before fitting or plotting.')
        parser.add_argument('--rses-norm-file',action='store_true',default=False,help='Normalize the RSES data to values stored in a file before fitting or plotting')
        parser.add_argument('--rses-norm-file-label',default='rses',help='Set the file label of the normalization data file. E.g. rses sets the exact normalizations, DRR431 sets the normalizations calculated from the entanglement Hamiltonian model.')
        
        ##  Counting options
        parser.add_argument('--counting-cut',default=0,type=float,help='Set the entanglement energy cut below which to cout the RSES states')
        parser.add_argument('--counting-max-sector',default=0,type=int,help='Set the highest sector to get the coutning for')
        
        ##  Fitting and plotting options
        
        parser.add_argument('--min-fit-sector',default=6,type=int,help='Minimum number of sectors to use for fitting')
        parser.add_argument('--max-fit-sector',default=6,type=int,help='Maximum number of sectors to use for fitting')
        parser.add_argument('--plot-sector',default=8,type=int,help='Number of sectors to use for plotting')
        parser.add_argument('--highlight-branches',default=1,type=int,help='Specify this option to produce a plot with the branches highlighted (also turns on block diagonalization in entanglement Hamiltonian calculations)')
        parser.add_argument('--convergence-min',default=200,type=int,help='Specify minimum number of Monte Carlo steps in RSES calculation (in 1000s)')
        parser.add_argument('--convergence-max',default=6400,type=int,help='Specify maximum number of Monte Carlo steps in RSES calculation (in 1000s)')
        parser.add_argument('--plot-counting',action='store_true',default=False,help='Set this option to plot the full conting of entanglement energy states')
        parser.add_argument('--plot-degeneracy',action='store_true',default=False,help='Set this option to plot numbers counting only degenerate states')
        parser.add_argument('--plot-legend',action='store_true',default=False,help='Set this option in order to turn on the plot legend')
        
        
        ##  Set data i/o directories and external program names and directories
        
        parser.add_argument('--bin-path',default='//rscratch//scd51//bin',help='Set the path to the entanglement Hamiltonian and entanglement energy binaries directory.')
        parser.add_argument('--entanglement-hamiltonian-bin',default='fqhe_entanglement_hamiltonian',help='Set the program name for calcualtion of the entanglement Hamiltonian.')
        parser.add_argument('--entanglement-energy-bin',default='fqhe_entanglement_energy',help='Set the program name for calcualtion of the single particle entanglement energy.')
        parser.add_argument('--base-data-path',default='//scratch//scd51//physics_data',help='Set the first part of the common path to the input/output data directories.')
        parser.add_argument('--ewf-rses-path',default='ewf_spectrum_unnormalized',help='Set the name of the directory containing RSES data calculated by the EWF method. To be appended to --base-data-path.')
        parser.add_argument('--rdm-rses-path',default='entanglement_spectrum',help='Set the name of the directory containing RSES data calculated by the pseudo reduced density matrix method. To be appended to --base-data-path.')
        parser.add_argument('--exact-rses-path',default='exact_entanglement_spectrum',help='Set the name of the directory containing RSES data calculated by a numerically exact method. To be appended to --base-data-path.')
        parser.add_argument('--rses-norm-path',default='entanglement_energy_minima',help='Set the name of the directory containing the exact RSES normalizaiton data. To be appended to --base-data-path.')
        parser.add_argument('--output-path',default='rses_fitting',help='Set the name of the directory containing fitted RSES output data. To be appended to --base-data-path.')
        
        ##  Process the arguments

        args = parser.parse_args(namespace=self)
        
        ##  Check consistency of options
        
        if self.wf_type=='pef' and self.nbr_ll==1:
            print('\n\nERROR - set nbr_ll=1 with pef type')
            quit()
            
        if self.wf_type=='pef' and self.jastrow==1 and self.statistics=='fermions':
            print('\n\nERROR - Set fermion statistics with wrong Jastrow power')
            quit() 
        
        if self.wf_type=='laughlin' and self.jastrow % 2 == 1 and self.statistics=='fermions':
            self.statistics='bosons'
        
        if self.wf_type=='laughlin' and self.jastrow % 2 == 0 and self.statistics=='bosons':
            self.statistics='fermions'
        
        if self.wf_type=='pef' and self.jastrow==2 and self.statistics=='fermions':
            print('\n\nERROR - Set fermion statistics with wrong Jastrow power')
            quit()   
            
        if self.fit_rses and self.nbr_ll>1 and self.fitting_model=='DRR431':
            print('\n\nERROR - Using DRR431 Hamiltonian to fit multi-Landau level case!')
            quit()

        if self.fit_offset and self.fitting_model=='entanglement_energy':
            self.rses_norm_zero = False
            print('\n\nWARNING - disabling rses_norm_zero in order to fit offset')
        
        if self.rses_norm_zero and self.rses_norm_file:
            self.rses_norm_zero = False
            print('\n\nWARNING - disabling rses_norm_zero in order to use rses_norm_file option')

        self.GeneratePossibleNa()
        
        ##  Print a summary

        self.Print()

    ############################################################################
    ##  A function to generate a list of possible NA values to be used in the
    ##  combined fitting
    ##
    def GeneratePossibleNa(self):
        
        ##  Define one more argument and populate it
        self.possibleNa = []
        
        self.possibleNa.append(self.nbr_a)
        
        if self.extra_nbr_a > 0:
            
            self.possibleNa = self.possibleNa + range(self.nbr_a+1,min(self.nbr+1,self.nbr_a+self.extra_nbr_a+1))
            
            self.possibleNa = range(max(0,self.nbr_a-self.extra_nbr_a),self.nbr_a) + self.possibleNa
            
    ############################################################################
    ##  Determine which index in possibleNa corresponds to the given nbr_a
    ##
    def GetPossibleNaIndex(self,nbr_a):
    
        k = 0 
        
        for i in self.possibleNa:
        
            if i == nbr_a:
            
                return k
                
            k += 1
        
        ##  If not found in the list return 1 + the length of possibleNa
        return k
        
    ############################################################################
    ##  Function to print all program option parameters
    ##
    def Print(self):
        print ('================================================================================')
        print ('================================================================================')
        print ('\tUsing Scipy version: '+scipy.__version__)
        print ('\tWave function type: '+str(self.wf_type))
        print ('\tNumber of electrons: '+str(self.nbr))
        print ('\tNumber of particles in the A cut: '+str(self.possibleNa))
        print ('\tNumber of effective Landau levels: '+str(self.nbr_ll))
        print ('\tPower in the Jastrow factor: '+str(self.jastrow))
        print ('\tFilling factor: '+str(self.nbr_ll)+'/'+str(self.jastrow*self.nbr_ll+1))
        print ('\tStatistics: '+str(self.statistics))
        print ('\tNormalizing RSES data to zero? '+str(self.rses_norm_zero))
        print ('\tNormalizing RSES data to file value? '+str(self.rses_norm_file))
        if self.rses_norm_file:
            print ('\tRSES norm data path:\n\t\t'+self.base_data_path+'//'+self.rses_norm_path)
            print ('\tRSES norm data label: '+str(self.rses_norm_file_label))
        print ('\tFitting RSES offsets? '+str(self.fit_offset))
        if self.rses_method == 'ewf':
            print ('\tMonte Carlo steps used in RSES calculation: '+str(self.rses_steps)+'000')
            print ('\tRSES data path:\n\t\t'+self.base_data_path+'//'+self.ewf_rses_path)
        elif self.rses_method == 'rdm':
            print ('\tDimension of pseudo-reduced density matrix in RSES calculation: '+str(self.rses_size))
            print ('\tNumber of RSES sectors skipped: '+str(self.skip_sector))
            print ('\tRSES data path:\n\t\t'+self.base_data_path+'//'+self.rdm_rses_path)
        elif self.rses_method == 'exact':
            print ('\tPerforming exact calculation of RSES.')
            print ('\tRSES data path:\n\t\t'+self.base_data_path+'//'+self.exact_rses_path)
        if self.fitting_model=='DRR431' or self.fitting_model=='extended_DRR':
            print ('\tFitting entanglement Hamiltonian: '+str(self.fitting_model))
            print ('\tLooking for entanglement Hamiltonian binary in:\n\t\t'+str(self.bin_path)+'//'+str(self.entanglement_hamiltonian_bin))
        elif self.fitting_model=='entanglement_energy':
            print ('\tFitting single particle entanglement energy model')
            print ('\tLooking for entanglement energy binary in:\n\t\t'+str(self.bin_path)+'//'+str(self.entanglement_energy_bin))
        print ('\tMinimum sector used for fitting: '+str(self.min_fit_sector))
        print ('\tMaximum sector used for fitting: '+str(self.max_fit_sector))
        if self.make_plots:
            print ('\tMaximum sector used for plotting: '+str(self.plot_sector))
            print ('\tHighlight branches in the plot: '+str(self.highlight_branches))
        print ('\tUse block diagonalization: '+str(self.highlight_branches))
        print ('\tOutput directory:\n\t\t'+self.BuildOutDirectoryName())
        print ('================================================================================')
        print ('================================================================================')
    
    ############################################################################
    ##  Generate the output directory name
    ##
    def BuildOutDirectoryName(self):
        
        dir = self.base_data_path+'//'+self.output_path+'//'+self.wf_type+'_'+str(self.nbr_ll)+'_'+str(self.jastrow*self.nbr_ll+1)+'_n_'+str(self.nbr)+'//'
        
        ##  Check for the existence of the directory
        
        if not os.path.isdir(dir):
            print("ERROR: could not find directory "+dir)
            
            quit()

        return dir
        
    ############################################################################
    ##  Generate the RSES input file name
    ##
    def BuildRsesFileName(self,sector,nbr_a):
        
        ##  Set the file name for rdm type data
        if self.rses_method == 'rdm':
        
            return self.base_data_path+'//'+self.rdm_rses_path+'/'+self.wf_type+'_'+str(self.nbr_ll)+'_'+str(self.jastrow*self.nbr_ll+1)+'_n_'+str(self.nbr)+'_size_'+str(self.rses_size)+'x'+str(2*self.rses_size)+'_sector_'+str(sector)+'_eigenvalues.dat'
        
        ##  Set the file name for exact type data
        if self.rses_method == 'exact':
        
            return self.base_data_path+'//'+self.exact_rses_path+'/Exact_IQHE_Eigenvalues_n_'+str(self.nbr)+'_na_'+str(nbr_a)+'_LL_'+str(self.nbr_ll)+'.dat'
        
        ##  Set the file name for ewf type data
        if self.rses_method == 'ewf':
        
            if self.wf_type=='laughlin':
                extraFlux = 1
            else:
                extraFlux = 0
                
            ##  Remove a kludge in file names
            #   jastrow=0 indicates IQHE case in general, but if LL=1 also then we need jastrow=1
            
            if self.jastrow==0 and self.nbr_ll==1:
            
                extraFlux = 1
        
            return self.base_data_path+'//'+self.ewf_rses_path+'/Eigenvalues_n_'+str(self.nbr)+'_na_'+str(nbr_a)+'_flux_'+str(self.jastrow+extraFlux)+'_LL_'+str(self.nbr_ll)+'_iter_'+str(self.rses_steps)+'000.dat'

    ############################################################################
    ##  Generate the RSES normalization file name
    ##
    def BuildRsesNormFileName(self,label,nbr_a):

        return self.base_data_path+'//'+self.rses_norm_path+"/minimum_"+str(label)+'_n_'+str(programOptions.nbr)+'_na_'+str(nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'.dat'

    ############################################################################
    ##  Generate the fit model data output file name for a given nbr_a
    ##
    def BuildFitModelFileName(self,sector,nbr_a):
        
        fileName = self.BuildOutDirectoryName()+'/eigenvalues_'+self.fitting_model+'_model_n_'+str(self.nbr)+'_na_'+str(nbr_a)+'_ll_'+str(self.nbr_ll)+'_sector_'+str(sector)+'.dat'

        ##  Fitting for an entanglement Hamiltonian as in Eq 4.31 of Dubail et al.   
        ##  (see PRB 86 245310)    
        if self.fitting_model=='DRR431':

            return fileName

        ##  Fitting for an entanglement Hamiltonian that is two copies of the 
        ##  DRR431 (one for each CF LL)
        elif self.fitting_model=='extended_DRR':

            if self.nbr_ll == 1:
                
                print ('ERROR: extended_DRR type Hamiltonian cannot be used with nbr_lls set to 1')
                quit()

            return fileName

      
        ##  Fitting a polynomial single particle entanglement energy
        elif self.fitting_model=='entanglement_energy':

            return fileName
        
        else:
        
            print ('Unknown fit model type specified: '+str(self.fitting_model))
            quit()
           
################################################################################
################################################################################
##  A class to contain plot layout specifications
##
class PlotSpecification:

    ############################################################################
    ##  Generate the default plot axis format
    ##
    def __init__(self):

        fontSize = 28
 
        self.insetOn = False
        self.legendOn = False

        self.main_axes = plot_utilities.plt.axes([0.01,0.01,0.75,0.6])	# ([left, bottom, width, height])

        ##  Remove the ticks from the top and left axes
        self.main_axes.get_xaxis().tick_bottom()
        self.main_axes.get_yaxis().tick_left()
        self.main_axes.spines['right'].set_visible(False)
        self.main_axes.spines['top'].set_visible(False)
        
        if programOptions.wf_type=='laughlin':
            
            self.insetOn = True
        
            self.inset_axes = plot_utilities.plt.axes([0.59,0.15,0.154,0.24])	# ([left, bottom, width, height])
        
            self.inset_axes.get_xaxis().tick_bottom()
            ##self.inset_axes.get_yaxis().tick_left()
            ##self.inset_axes.spines['right'].set_visible(False)
            ##self.inset_axes.spines['top'].set_visible(False)
        
        if programOptions.wf_type=='iqhe' and programOptions.nbr_ll==1:
            
            self.insetOn = False
        
            ##self.inset_axes = plot_utilities.plt.axes([0.59,0.2,0.154,0.42])	# ([left, bottom, width, height])
        
            ##self.inset_axes.get_xaxis().tick_bottom()
            ##self.inset_axes.get_yaxis().tick_left()
            ##self.inset_axes.spines['right'].set_visible(False)
            ##self.inset_axes.spines['top'].set_visible(False)

        ############    SPECIFY FORMATTING FOR SPECIAL CASES        ############

        ### N=12 1/2 Laughlin

        if programOptions.wf_type=='laughlin' and programOptions.nbr==12:
            self.main_axes.set_ylim(-2,22)             ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(0,25,5))
            self.main_axes.set_yticklabels([r'\textbf{0}',r'\textbf{5}',r'\textbf{10}',r'\textbf{15}',r'\textbf{20}'],fontsize=fontSize-10)
        
        ### N=50 1/2 LAUGHLIN
        
        if programOptions.wf_type=='laughlin' and programOptions.nbr==50:
            self.main_axes.set_ylim(-0.6,7)             ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(0,8))
            self.main_axes.set_yticklabels([r'\textbf{0}',r'\textbf{1}',r'\textbf{2}',r'\textbf{3}',r'\textbf{4}',r'\textbf{5}',r'\textbf{6}',r'\textbf{7}'],fontsize=fontSize-10)
            
            self.inset_axes.set_ylim(0.6,1.95)             ##  SET PLOT Y-LIMITS
            self.inset_axes.set_yticks([0.75,1.0,1.25,1.5,1.75])
            self.inset_axes.set_yticklabels([r'\textbf{0.75}',r'\textbf{1.0}',r'\textbf{1.25}',r'\textbf{1.5}',r'\textbf{1.75}'],fontsize=fontSize-14)
            
        ### N=50 1 IQHE
        
        if programOptions.wf_type=='iqhe' and programOptions.nbr==50:
            '''
            self.main_axes.set_ylim(-0.6,9)                ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(0,10))
            self.main_axes.set_yticklabels([r'\textbf{0}',r'\textbf{1}',r'\textbf{2}',r'\textbf{3}',r'\textbf{4}',r'\textbf{5}',r'\textbf{6}',r'\textbf{7}',r'\textbf{8}',r'\textbf{9}'],fontsize=fontSize-10)
            '''
            self.main_axes.set_ylim(-0.8,8)                ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(0,9))
            self.main_axes.set_yticklabels([r'\textbf{0}',r'\textbf{1}',r'\textbf{2}',r'\textbf{3}',r'\textbf{4}',r'\textbf{5}',r'\textbf{6}',r'\textbf{7}',r'\textbf{8}'],fontsize=fontSize-10)
            
            '''
            self.inset_axes.set_ylim(1.4,3.3)             ##  SET PLOT Y-LIMITS
            self.inset_axes.set_yticks([1.5,1.75,2.0,2.25,2.5,2.75,3.00,3.25])
            self.inset_axes.set_yticklabels([r'\textbf{1.5}',r'\textbf{1.75}',r'\textbf{2.0}',r'\textbf{2.25}',r'\textbf{2.5}',r'\textbf{2.75}',r'\textbf{3.0}',r'\textbf{3.25}'],fontsize=fontSize-14)
            '''
            '''
            self.inset_axes.set_ylim(1.6,3.1)             ##  SET PLOT Y-LIMITS
            self.inset_axes.set_yticks([1.75,2.0,2.25,2.5,2.75,3.00])
            self.inset_axes.set_yticklabels([r'\textbf{1.75}',r'\textbf{2.0}',r'\textbf{2.25}',r'\textbf{2.5}',r'\textbf{2.75}',r'\textbf{3.0}'],fontsize=fontSize-14)
            '''

        ### N=48 2 IQHE
        
        if programOptions.wf_type=='iqhe' and programOptions.nbr==48:
            self.main_axes.set_ylim(-8,23)                ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(-5,25,5))
            self.main_axes.set_yticklabels([r'$-$\textbf{5}',r'\textbf{0}',r'\textbf{5}',r'\textbf{10}',r'\textbf{15}',r'\textbf{20}'],fontsize=fontSize-10)

        ### N=46 2/3 JAIN
        
        if programOptions.wf_type=='pef' and programOptions.nbr==46:
        
            self.main_axes.set_ylim(-9.9,30)             ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(-5,35,5))
            self.main_axes.set_yticklabels([r'\textbf{$-$5}',r'\textbf{0}',r'\textbf{5}',r'\textbf{10}',r'\textbf{15}',r'\textbf{20}',r'\textbf{25}',r'\textbf{30}'],fontsize=fontSize-10)

        ### N=48 2/3 JAIN
            
        if programOptions.wf_type=='pef' and programOptions.nbr==48:
        
            self.main_axes.set_ylim(-8,35)             ##  SET PLOT Y-LIMITS
            self.main_axes.set_yticks(range(-5,40,5))
            self.main_axes.set_yticklabels([r'\textbf{$-$5}',r'\textbf{0}',r'\textbf{5}',r'\textbf{10}',r'\textbf{15}',r'\textbf{20}',r'\textbf{25}',r'\textbf{30}',r'\textbf{35}'],fontsize=fontSize-10)
        
        ### N=60 3/5 JAIN
        
        ##self.main_axes.set_ylim(-19,36)             ##  SET PLOT Y-LIMITS
        ##self.main_axes.set_yticks(range(-15,40,5),(r'\textbf{$-$15}',r'\textbf{$-$10}',r'\textbf{$-$5}',r'\textbf{0}',r'\textbf{5}',r'\textbf{10}',r'\textbf{15}',r'\textbf{20}',r'\textbf{25}',r'\textbf{30}',r'\textbf{35}'))
        
        ########################################################################

        self.main_axes.set_ylabel(r'$\Delta\xi$',rotation='horizontal',fontsize=fontSize-10)
        
        self.main_axes.yaxis.set_label_coords(0.02,1.01)
        self.main_axes.xaxis.set_label_coords(1.03,-0.02)

        self.main_axes.set_xlim(-programOptions.plot_sector-0.5,0.8)
        self.main_axes.set_xlabel(r'$\Delta L_z^A$',fontsize=fontSize-10)
        self.main_axes.set_xticks(range(0,-programOptions.plot_sector-1,-1))
        ##self.main_axes.set_xticklabels([r'\textbf{0}',r'\textbf{$-$1}',r'\textbf{$-$2}',r'\textbf{$-$3}',r'\textbf{$-$4}',
        ##r'\textbf{$-$5}',r'\textbf{$-$6}',r'\textbf{$-$7}',r'\textbf{$-$8}',r'\textbf{$-$9}',r'\textbf{$-$10}',r'\textbf{$-$11}',r'\textbf{$-$12}',r'\textbf{$-$13}',r'\textbf{$-$14}'],fontsize=fontSize-10)
        self.main_axes.set_xticklabels([r'\textbf{0}',r'\textbf{1}',r'\textbf{2}',r'\textbf{3}',r'\textbf{4}',
        r'\textbf{5}',r'\textbf{6}',r'\textbf{7}',r'\textbf{8}',r'\textbf{9}',r'\textbf{10}',r'\textbf{11}',r'\textbf{12}',r'\textbf{13}',r'\textbf{14}'],fontsize=fontSize-10)
        
        if self.insetOn:
            
            '''
            if programOptions.wf_type=='iqhe':
                self.inset_axes.set_xlim(-6.6,-3.25)
            else:
            '''
            self.inset_axes.set_xlim(-6.4,-3.2)
            self.inset_axes.set_xticks([-4,-5,-6])
            self.inset_axes.set_xticklabels([r'\textbf{4}',r'\textbf{5}',r'\textbf{6}'],fontsize=fontSize-14)
            
            self.inset_axes.yaxis.set_label_coords(0.02,1.01)
            self.inset_axes.xaxis.set_label_coords(1.03,-0.02)
    
    ############################################################################
    ##  Finalize the plots and save the output
    ##
    def Finalize(self,fileLabel,fitSector):
        
        fontSize = 28
        
        if self.legendOn:
            leg = self.main_axes.legend(loc='best',scatterpoints=1,markerscale=2,fontsize=fontSize-10)       

            ##  Don't colour the legend
            
            if programOptions.highlight_branches and programOptions.nbr_ll>1:
            
                for l in leg.legendHandles:
                    l.set_color('black')
        
        if self.insetOn:
            # 	draw a box of the region of the inset axes in the parent axes and
            #	connecting lines between the bbox and the inset axes area
            #	ec sets transparency of connecting lines
            #	loc 1 and loc 2 specify where they are connected
            mark_inset(self.main_axes,self.inset_axes, loc1=2, loc2=4, fc="none", ec="0.5")
        
        if fitSector==-1:
        
            fileName=programOptions.BuildOutDirectoryName()+'plot_'+str(fileLabel)+'_n_'+str(programOptions.nbr)+'_nA_'+str(nbr_a)+'.pdf'
        
        else:
        
            fileName=programOptions.BuildOutDirectoryName()+'plot_'+str(fileLabel)+'_'+programOptions.fitting_model+'_model_n_'+str(programOptions.nbr)+'_nA_'+str(nbr_a)+'_fit_'+str(fitSector)+'.pdf'
        
        self.main_axes.get_figure().savefig(fileName,bbox_inches='tight')

        print ('Generated plot '+fileName)
     
    ############################################################################
    ##  Clean up the object when it goes out of scope
    ##
    def __del__(self):
        plot_utilities.plt.close()
        
################################################################################
################################################################################
######      MAIN        ########################################################

if __name__=='__main__':

    ##  Declare an instance of the programOptions class and initialize from 
    ##  command line arguments

    programOptions = ProgramOptions()
    
    ##  Initialize the RSES data
        
    if programOptions.rses_method=="ewf":
        rses = EwfRses()
    elif programOptions.rses_method=="rdm":
        rses = RdmRses()
    elif programOptions.rses_method=="exact":
        rses = ExactRses() 
    else:
        print("ERROR: unknown RSES type: "+programOptions.rses_method)
        quit()
    
    if programOptions.get_counting:
        
        ##  Read in RSES data
        
        rses.GetSpectrumFromFile(programOptions,programOptions.counting_max_sector)
        
        print(rses.GetCounting(programOptions.counting_cut))
        
        quit()
        
    ############################################################################
    ######    GENERATE FILES OF EXACT RSES MINIMA DATA        ##################

    if programOptions.generate_minima_file:
        
        ##  Don't normalize the spectrum to zero or with values from a file
        programOptions.rses_norm_zero=False
        programOptions.rses_norm_file=False
        
        ##  Read in RSES data from the 0th sector
        
        rses.GetSpectrumFromFile(programOptions,0)
        
        ##  Output the minium entanglement energy of sector 0 in a new file
        
        rses.EntanglementEnergyMinimaToFile(programOptions)
    
    ############################################################################
        
    ##  Initialize the fitting model

    if programOptions.fitting_model=='entanglement_energy':
        model = SingleParticleEnergy(programOptions)
    elif programOptions.fitting_model=='DRR431':
        model = DRR431Hamiltonian(programOptions)
    elif programOptions.fitting_model=='extended_DRR':
        model = ExtendedDRRHamiltonian(programOptions)
    else:
        print("ERROR: unknown fitting model type: "+programOptions.fitting_model)
        quit()
    
    ############################################################################
    ######    FIT MODEL TO RSES DATA        ####################################
    if programOptions.fit_rses:
            
        print ("\n\tFITTING "+programOptions.fitting_model+" MODEL TO RSES DATA FOR SECTORS "+str(programOptions.min_fit_sector)+" TO "+str(programOptions.max_fit_sector))
        
        for fitSector in range(programOptions.min_fit_sector,programOptions.max_fit_sector+1):
        
            ##  Read in RSES data to be fitted

            rses.GetSpectrumFromFile(programOptions,fitSector)

            ##  Generate the counting of states in the spectrum and use that
            ##  to remove the continuum of states in the RSES. The RSES
            ##  is calculated for an initial random set of fitting parameters

            model.GenerateModelSpectrum(programOptions,fitSector)
            
            model.GetSpectrumFromFile(programOptions,fitSector,False,{})
            
            rses.RemoveContinuum(model.allSectors)
            
            ##print(rses.allSectors)
            ##print(rses.allEigenvalues)
            ##raw_input("Press Enter to continue...")
            
            ##  Optimize the fitting parameter values

            fittingResult = minimize(model.MinFunction,model.initialFittingParameters,args=(rses,programOptions,fitSector,),method='powell',options={'xtol': 1e-6, 'disp': True})
            
            model.StoreFinalFitParameters(fittingResult,programOptions,fitSector)
            
            ##  Store a copy of the final model parameters
            
            ##  In case of entanglement Hamitlonian models, set the offset
            ##  parameter so that the spectrum is normalized
            
            if programOptions.fitting_model=='DRR431' or programOptions.fitting_model=='extended_DRR':
            
                model.RenormalizeSpectrum(fitSector,programOptions)
            
                ##  Store an updated copy of the final model parameters
            
                model.StoreFinalFitParameters(fittingResult,programOptions,fitSector)
            
            ##  Output normalization data to a file
            
            model.RenormalizeSpectrum(fitSector,programOptions)
            
            model.ClearWorkingFiles(programOptions)

    ############################################################################
    ######    GENERATE PLOTS        ############################################
    if programOptions.make_plots or programOptions.plot_rses:

        print ("\n\tGENERATING PLOTS FOR SECTORS "+str(programOptions.min_fit_sector)+" TO "+str(programOptions.max_fit_sector))
        
        ##  Read in RSES data to be plotted

        rses.GetSpectrumFromFile(programOptions,programOptions.plot_sector)
        
        for fitSector in range(programOptions.min_fit_sector,programOptions.max_fit_sector+1):
        
            if not programOptions.plot_rses:
        
                ##  Obtain fitting parameter values for the given fitSector value
                
                model.GetFittingParametersFromFile(programOptions,fitSector)
                
                ##  Re-calculate the final spectrum for plotting
                
                print ("\n\tRecalculating spectrum...")

                model.GenerateModelSpectrum(programOptions,programOptions.plot_sector)
                
                model.GetSpectrumFromFile(programOptions,programOptions.plot_sector,False,model.CalculateOffsets(True,programOptions))
                
                model.ClearWorkingFiles(programOptions)
            
            #####   MAKE RSES PLOTS  ######
          
            ##	Generate a full plot of the fitted data, superimposed on 
            ##  a full plot of the RSES (including the upper continuum)
            '''
            for nbr_a in programOptions.possibleNa:
                
                plotSpec = PlotSpecification()
                
                ##  Add RSES data
                
                rses.Plot(plotSpec.main_axes,nbr_a,False,programOptions)
                
                ##  Add fitting data
                
                model.Plot(plotSpec.main_axes,nbr_a,fitSector,False,programOptions)
                
                ##  Finalize the plot
                
                plotSpec.Finalize('full_fitted_spectrum',fitSector)
                
                del plotSpec
            '''
            ##  Plot the RSES alone
            
            for nbr_a in programOptions.possibleNa:
            
                plotSpec = PlotSpecification()
                
                ##  Add RSES data
                
                rses.Plot(plotSpec.main_axes,nbr_a,False,programOptions)
                
                if plotSpec.insetOn:
                    rses.Plot(plotSpec.inset_axes,nbr_a,True,programOptions)
                
                ##  Finalize the plot
                
                plotSpec.Finalize('rses',-1)
                
                del plotSpec
            
            ##	Generate a full plot of the fitted data, superimposed on 
            ##  a full plot of the RSES (NOT including the upper continuum)
            
            if not programOptions.plot_rses:
            
                rses.RemoveContinuum(model.allSectors)
                
                for nbr_a in programOptions.possibleNa:
                    
                    plotSpec = PlotSpecification()
                    
                    if programOptions.plot_legend:
                        plotSpec.legendOn = True
                    
                    ##  Add RSES data
                    
                    rses.Plot(plotSpec.main_axes,nbr_a,False,programOptions)
                    if plotSpec.insetOn:
                        rses.Plot(plotSpec.inset_axes,nbr_a,True,programOptions)
                    
                    ##  Add fitting data
                    
                    model.Plot(plotSpec.main_axes,nbr_a,fitSector,False,programOptions)
                    if plotSpec.insetOn:
                        model.Plot(plotSpec.inset_axes,nbr_a,fitSector,True,programOptions)
                    
                    ##  Finalize the plot
                    
                    plotSpec.Finalize('fitted_spectrum',fitSector)
        
                #####   MAKE MINIMUM ENTANGLEMENT ENERGY PLOTS  ######
                
                ##  If applicable, check the spectrum normalization and output
                ##  the lowest entanglement eigenvalue to a file
                
                if programOptions.extra_nbr_a > 0:
                
                    ##  Obtain mimimum fitted entanglement energy values from an existing file
                    
                    fitNorms = []
                    
                    for nbr_a in programOptions.possibleNa:

                        print ('\n\tReading in RSES norm data from '+programOptions.BuildRsesNormFileName(programOptions.fitting_model,nbr_a)+'\n')

                        fin=open(programOptions.BuildRsesNormFileName(programOptions.fitting_model,nbr_a))
                        
                        ##  Skip the comment line        
                        fin.readline()
                        
                        ##  Get the norm value from the file
                        fitNorms.append(float(fin.readline()))
                            
                        fin.close()
                
                    programOptions.rses_norm_zero = False
                    print('WARNING - disabling rses_norm_zero in order to plot entanglement energy minima')

                    rsesNorms = rses.CalculateNormalizations(12,programOptions)
                    
                    main_axes = plot_utilities.plt.axes([0.1,0.1,0.75,0.75])	# ([left, bottom, width, height])
                    
                    main_axes.scatter(programOptions.possibleNa,fitNorms,label="Fitted minimum",marker='_',lw=1.0,s=100,c=plot_utilities.colour_list[3])
                    
                    main_axes.scatter(programOptions.possibleNa,rsesNorms,label="Actual minimum",marker='_',lw=1.0,s=100,c=plot_utilities.colour_list[4])
                    
                    main_axes.legend(loc='best',scatterpoints=1,markerscale=2)

                    main_axes.set_xlabel(r'$N_A$',**plot_utilities.font)
                    main_axes.set_ylabel(r'Minimum entanglement energy',**plot_utilities.font)

                    fileName = programOptions.BuildOutDirectoryName()+'plot_minimum_comparison_'+programOptions.fitting_model+'_model_n_'+str(programOptions.nbr)+'_fit_'+str(fitSector)+'.pdf'

                    plot_utilities.plt.savefig(fileName,bbox_inches='tight')
            
                    plot_utilities.plt.close()
                    
                    print ('\n\tGenerated plot '+fileName)    
        
    ############################################################################
    ######    ANALYSE AND PLOT RSES CONVERGENCE      ###########################
    if programOptions.plot_rses_convergence:

        print ('\tAnalysing RSES convergence ')

        ##    Generate the counting of states in the spectrum
        ##    and use that to remove the continuum of states in the RSES
        
        model.GenerateModelSpectrum(programOptions,programOptions.max_fit_sector)
        
        model.GetSpectrumFromFile(programOptions,programOptions.max_fit_sector,False,{})

        model.ClearWorkingFiles(programOptions)
        
        ##  Read in the RSES spectrum for each number of MC samples specified
        ##  Remove the continuum of each data set before appending to an overall 
        ##  list
        
        nbrSizes = int(math.log(programOptions.convergence_max/programOptions.convergence_min,2))
        
        programOptions.rses_steps = programOptions.convergence_min
        
        eigenvalues=[]
        
        print ('\tReading in RSES files... ')
        
        for i in range(0,nbrSizes+1):
        
            print ('\t#MC steps in RSES calculation: '+str(programOptions.rses_steps)+'000')

            rses.GetSpectrumFromFile(programOptions,programOptions.plot_sector)
            
            rses.RemoveContinuum(model.allSectors)
            
            eigenvalues.append(rses.allEigenvalues[0])
            
            programOptions.rses_steps *= 2
            
        ##    Calculate average differences in mean eigenvalue for each sector

        meanEigs = np.empty((programOptions.plot_sector,nbrSizes+1))

        fitSectors = model.allSectors[0]

        j=0

        print ('\tCalculating mean eigenvalues in each sector ')

        for row in eigenvalues:
            
            k=0
            currSector = 0
            counter    = 0
            currMean   = 0.0
            
            ##print(row)
            
            for eig in row :
            
                ##print (eig,k,fitSectors[k],currSector)
                
                if fitSectors[k]==currSector :

                    ##  add up all values for that sector  

                    currMean+=eig
                        
                    counter+=1    
                        
                else :
                   
                    ##print (counter)

                    currMean/=counter
                        
                    ##print(currDiff)
                   
                    meanEigs[-currSector,j] = currMean

                    ##    go to the next sector
                    
                    currSector-=1
                    currMean=0.0
                    counter=0

                    currMean+=eig
                        
                    counter+=1
                k+=1
            j+=1    

        print (meanEigs)

        print ('\n\tMaking eigenvalue tracking plot ')

        ##    Generate plot tracking only the mean values

        plotX=[]

        curr = programOptions.convergence_min

        for i in range(0,nbrSizes+1):

            plotX.append(math.log(curr*1000,2))
            
            curr *= 2

        main_axes = plot_utilities.plt.axes([0.12,0.12,0.85,0.85])    # ([left, bottom, width, height])

        #    remove ticks on the top and left axes
        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()
        
        main_axes.spines['right'].set_visible(False)
        main_axes.spines['top'].set_visible(False)

        counter=0

        for plotY in meanEigs:

            print(plotX,plotY)

            main_axes.plot(plotX,plotY,marker=plot_utilities.marker_list[counter],c=plot_utilities.colour_list[counter],label=str(counter))    
            #    (xdata, ydata, marker width, marker type, marker height)

            counter+=1
            
        main_axes.set_xlabel(r'Log (base 2) of no. MC steps',**plot_utilities.font)
        main_axes.set_ylabel(r'Mean $\xi$ per sector',**plot_utilities.font)
        main_axes.legend(loc='best',scatterpoints=1,numpoints=1,prop={'size':16},title=r'$|\Delta L_z^A|$')
        main_axes.set_xlim(math.floor(math.log(programOptions.convergence_min*1000,2)),math.ceil(math.log(programOptions.convergence_max*1000,2)))

        fileName = programOptions.BuildOutDirectoryName()+'track_rses_means_'+str(programOptions.nbr_ll)+'_'+str(programOptions.nbr_ll*programOptions.jastrow+1)+'_n_'+str(programOptions.nbr)+'_nA_'+str(programOptions.nbr_a)+'.pdf'

        plot_utilities.plt.savefig(fileName,bbox_inches='tight')
        
        plot_utilities.plt.close()
        
        print ('\n\tGenerated plot '+fileName)    
        
        print ('\n\tCalcualting and plotting mean differences between eigenvalues')
        
        ##  Now determine the differences in the mean eigenvalues for each sector
            
        diffs=np.empty((programOptions.plot_sector,nbrSizes))

        for i in range(0,programOptions.plot_sector):

            for j in range(nbrSizes):

                diffs[i,j] = meanEigs[i,j]-meanEigs[i,j+1]
                
        ##print (diffs)      

        ##    Generate plot tracking the average differences

        plotX=[]

        curr = programOptions.convergence_min

        for i in range(0,nbrSizes):

            plotX.append(math.log(curr*1000,2))
            
            curr *= 2 

        main_axes = plot_utilities.plt.axes([0.18,0.12,0.80,0.85])    # ([left, bottom, width, height])

        #    remove ticks on the top and left axes
        main_axes.get_xaxis().tick_bottom()
        main_axes.get_yaxis().tick_left()

        main_axes.spines['right'].set_visible(False)
        main_axes.spines['top'].set_visible(False)
        
        counter=0

        for plotY in diffs:

            ##print (plotY,plotX)

            main_axes.plot(plotX,plotY,marker=plot_utilities.marker_list[counter],c=plot_utilities.colour_list[counter],label=str(counter))    
            #    (xdata, ydata, marker width, marker type, marker height)

            counter+=1

        main_axes.set_xlabel(r'Log (base 2) of no. MC steps',**plot_utilities.font)
        main_axes.set_ylabel(r'Mean difference in $\xi$ per sector',**plot_utilities.font)
        main_axes.legend(loc='best',scatterpoints=1,numpoints=1,prop={'size':16},title=r'$|\Delta L_z^A|$')
        main_axes.set_xlim(math.floor(math.log(programOptions.convergence_min*1000,2)),math.ceil(math.log(programOptions.convergence_max*1000,2)))
        
        main_axes.set_xticks(range(18,23,1),(r'\textbf{18}',r'\textbf{19}',r'\textbf{20}',r'\textbf{21}',r'\textbf{22}',r'\textbf{23}'))
        
        if programOptions.nbr==12 and programOptions.wf_type=='laughlin':
        
            main_axes.set_yticks([-0.002,0.00,0.002,0.004,0.006],(r'\textbf{$-$0.002}',r'\textbf{0}',r'\textbf{0.002}',r'\textbf{0.004}',r'\textbf{0.006}'))
            
        if programOptions.nbr==50 and programOptions.wf_type=='laughlin':
        
           main_axes.set_yticks([-0.06,-0.04,-0.02,0,0.02,0.04,0.06],(r'\textbf{$-$0.06}',r'\textbf{$-$0.04}',r'\textbf{$-$0.02}',r'\textbf{0}',r'\textbf{0.02}',r'\textbf{0.04}',r'\textbf{0.06}'))
           
        if programOptions.nbr==46 and programOptions.wf_type=='pef':
        
           main_axes.set_yticks([-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.20],(r'\textbf{$-$0.15}',r'\textbf{$-$0.10}',r'\textbf{$-$0.05}',r'\textbf{0}',r'\textbf{0.05}',r'\textbf{0.10}',r'\textbf{0.15}',r'\textbf{0.20}'))
           
        if programOptions.nbr==48 and programOptions.wf_type=='pef':
        
           main_axes.set_yticks([-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5],(r'\textbf{$-$0.4}',r'\textbf{$-$0.3}',r'\textbf{$-$0.2}',r'\textbf{$-$0.1}',r'\textbf{0}',r'\textbf{0.1}',r'\textbf{0.2}',r'\textbf{0.3}',r'\textbf{0.4}',r'\textbf{0.5}'))
            

        fileName = programOptions.BuildOutDirectoryName()+'rses_convergence_'+str(programOptions.nbr_ll)+'_'+str(programOptions.nbr_ll*programOptions.jastrow+1)+'_n_'+str(programOptions.nbr)+'_nA_'+str(programOptions.nbr_a)+'.pdf'

        plot_utilities.plt.savefig(fileName,bbox_inches='tight')

        print ('\n\tGenerated plot '+fileName)

    ############################################################################
    ######    ANALYSE AND PLOT FITTING CONVERGENCE      ########################
    
    if programOptions.plot_fitting_convergence:
    
        print ("\n\tANALYSING FITTING CONVERGENCE FOR SECTORS "+str(programOptions.min_fit_sector)+" TO "+str(programOptions.max_fit_sector))
        
        extrapList = []
        
        sectorList = range(programOptions.min_fit_sector,programOptions.max_fit_sector+1)
        
        for fitSector in sectorList:
            
            ##  Obtain fitting parameter values for the given fitSector value
            
            model.GetFittingParametersFromFile(programOptions,fitSector)
            
            ##  Store fitting parameters in a set of lists
            
            extrapList.append(model.modelParameterList.copy())
            
        ##  Make plots of the fitting parameter convergence, including the
        ##  Extrapolated parameter values
        
        for parameter in extrapList[0].keys():

            valList = []
            
            for vals in extrapList:
                valList.append(vals[parameter])
        
            ##  Ignore if all values are zero (since that means the parameter
            ##  was not included in the fitting procedure)
            
            if sum(valList) == 0.0:
                pass
            else:

                print("Plotting "+parameter+" convergence")
                
                main_axes = plot_utilities.plt.axes([0.15,0.09,0.75,0.85])	# ([left, bottom, width, height])
                
                fitData,error,extrap = curve_fitting.ExponentialFit(sectorList,valList,100000)

                message = "Extrapolated $"+parameter+" = "+str(extrap)+" \\pm "+str(error)+"$"

                ##print("\n\t"+message)

                main_axes.set_title(message)

                main_axes.scatter(sectorList,valList,s=50,c=plot_utilities.colour_list[1],marker='o',label=parameter,edgecolor='none')
                
                main_axes.plot(sectorList,fitData,c=plot_utilities.colour_list[1])

                main_axes.set_ylabel(parameter,**plot_utilities.font)
                main_axes.set_xlabel(r'Maximum $\Delta L_z^A$ fitted',**plot_utilities.font)
                main_axes.set_xticks(sectorList)
                
                fileName = programOptions.BuildOutDirectoryName()+programOptions.fitting_model+'_model_'+parameter+'_convergence_'+str(programOptions.nbr_ll)+'_'+str(programOptions.nbr_ll*programOptions.jastrow+1)+'_n_'+str(programOptions.nbr)+'_nA_'+str(programOptions.nbr_a)+'.pdf'

                plot_utilities.plt.savefig(fileName,bbox_inches='tight')

                plot_utilities.plt.close()

                print ('\n\tGenerated plot '+fileName)
                   
    ############################################################################
    ######    ANALYSE AND PLOT THERMODYNAMIC CONVERGENCE      ##################
    
    if programOptions.plot_thermodynamic_convergence:
    
        print ("\n\tANALYSING THERMODYNAMIC FITTING CONVERGENCE")
    
    ############################################################################
    ######    CALCULATE VON NEUMANN/RENYI ENTROPY       ########################
    
    if programOptions.calculate_entanglement_entropy:

        programOptions.rses_norm_zero = False
        print('WARNING - disabling rses_norm_zero in order to calculate entanglement entropy (Minimum, Von Neumann and Renyi)')
        
        ##  Fix the sector and extra nbr a values used in the norm calculation
        test_sector = 20
        programOptions.extra_nbr_a = 4
        programOptions.GeneratePossibleNa()
        
        model.RenormalizeSpectrum(programOptions.max_fit_sector,programOptions)

        print ("\n\tCALCULATING MINIMUM ENTANGLEMENT ENTROPY")

        minima,total = model.CalculateNormalizations(test_sector,programOptions)

        model.ClearWorkingFiles(programOptions)
    
        print ("\n\tCALCULATING VON NEUMANN/RENYI ENTROPY")
        
        for q in range(1,10):
        
            model.CalculateRenyiEntropy(q,test_sector,programOptions.max_fit_sector,programOptions)

        model.ClearWorkingFiles(programOptions)

    ############################################################################

    print ('\n\tALL TASKS COMPLETED') 
