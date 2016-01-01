#! /usr/bin/env python2.7
################################################################################
##                                                                            ##
##						  Author: Simon C. Davenport						  ##
##                                                                            ##
##							Last Modified: 22/01/2015						  ##
##                                                                            ##
##                                                                            ##
##		This file contains a class to handle reading in and processing        ##
##      RSES fitting model data, using a single particle entanglement         ##
##      energy model or an entanglement Hamiltonian model                     ##
##                                                                            ##
################################################################################

################################################################################
################################################################################
##  Import modules
import os
import random
import math
from utilities import plot_utilities
from mpmath import mp

################################################################################
################################################################################
##  Define a base class for calculating and importing fitted RSES data of
##  various types
##
class FittingModel():

    ############################################################################
    ##  Base class constructor
    ##
    def __init__(self):
        self.allEigenvalues = []
        self.allLabels = []
        self.allSectors = []
        ##  Number of actual fitting parameters used
        self.nbrFittingParameters = 0
        self.initialFittingParameters = []
        ##  List of all possible parameters
        self.modelParameterList = {}
    
    ############################################################################
    ##  Abstract GenerateModelSpectrum method
    ##
    def GenerateModelSpectrum(self,programOptions,fitSector):
        raise NotImplementedError("Subclass must implement abstract method")
    
    ############################################################################
    ##  Abstract MinFunction method
    ##
    def MinFunction(self,rawFittingParameters,rsesData,programOptions,fitSector):
        raise NotImplementedError("Subclass must implement abstract method")
    
    ############################################################################
    ##  Abstract StoreFinalFitParameters method
    ##
    def StoreFinalFitParameters(self,fittingResult,programOptions,fitSector):
        raise NotImplementedError("Subclass must implement abstract method")
    
    ############################################################################
    ##  Abstract GetFittingParametersFromFile method
    ##
    def GetFittingParametersFromFile(self,programOptions,fitSector):
        raise NotImplementedError("Subclass must implement abstract method")
    
    ############################################################################
    ##  Abstract ClearWorkingFiles method
    ##
    def ClearWorkingFiles(self,programOptions):
        raise NotImplementedError("Subclass must implement abstract method")
        
    ############################################################################
    ##  CalculateOffsets default method 
    ##
    def CalculateOffsets(self,includeLinear,programOptions):
        return {}
    
    ############################################################################
    ##	A function to evaluate a least-squares fitting parameter. 
    ##	This is done by taking the lowest eigenvalues in each sector 
    ##	and matching these up. The fitting parameter is the sum of 
    ##	the squares of the differences in equivalent eigenvalues.	
    ##
    ##  This method has been modified on 07/01/2014 such that the parameter
    ##  value is weighted by the number of eigenvalues in a given sector.
    ##  This avoids bias towards higher sectors, especially in cases where
    ##  there are a large number of eigenvalues in higher sectors. 
    ##	
    ##	NB this method will not work if the counting of the fitting 
    ##	spectrum is not the same as the spectrum we're trying to fit
    ##
    ##  This method has been updated on 23/01/2014 such that there is an
    ##  additional weight of exp(-eigenvalue). This makes it so that the 
    ##  lowest lying states in  the spectrum are given more weight
    ##
    def FittingParameter(self,rsesData,programOptions):
	
        value = 0

        if programOptions.wf_type == 'laughlin' or programOptions.nbr_ll==1:
            for i in range(len(self.allEigenvalues)):
                for j in range(len(self.allEigenvalues[i])):
                    value+=((rsesData.allEigenvalues[i][j]-self.allEigenvalues[i][j])**2)/(rsesData.allFitFunctionWeights[i][j]**2)*math.exp(-rsesData.allEigenvalues[i][j])
        else:
            for i in range(len(self.allEigenvalues)):
                for j in range(len(self.allEigenvalues[i])):
                    value+=((rsesData.allEigenvalues[i][j]-self.allEigenvalues[i][j])**2)/(rsesData.allFitFunctionWeights[i][j]**2)
	
        print ('=======================================\n\tFitting parameter '+str(value)+'\n=======================================')

        return value

    ############################################################################
    ##	A function to read in the sector data from the entanglement 
    ##	Hamiltonian eigenvalue files. 
    ##
    def GetSpectrumFromFile(self,programOptions,maxSector,normalize,offsetList):
	
        self.allEigenvalues = []
        self.allLabels      = []
        self.allSectors     = []

        ##  Turn off offsets if the given list is empty

        for nbr_a in programOptions.possibleNa:

            #	find first sector file

            condition=1
            sector=0

            while condition==1:
                if os.path.isfile(programOptions.BuildFitModelFileName(sector,nbr_a)):
                    condition=0
                    break
                else:
                    ##print(programOptions.BuildFitModelFileName(sector,nbr_a))

                    condition=1
                    sector+=1	

	            if sector>400:
	
		            print ('\tERROR: CANNOT FIND FILES IN  '+programOptions.BuildOutDirectoryName())
	
		            quit()

            eigenvalues=[]
            sectors=[]
            labels=[]

            while condition!=1:

                fin=open(programOptions.BuildFitModelFileName(sector,nbr_a))

                for line in fin:

                    columns = line.split('\t')

                    eigenvalues.append(float(columns[0]))
                    labels.append(float(columns[1]))
                    sectors.append(-sector)

                fin.close()

                ##print ('read sector '+str(sector))

                sector+=1

                #   check for next file

                if os.path.isfile(programOptions.BuildFitModelFileName(sector,nbr_a)):
                    condition=0
                else:
                    condition=1  

                if sector>maxSector:
                    condition = 1

            ##  Append to list of eigenvalues.
        
		    self.allEigenvalues.append(eigenvalues)
		    self.allLabels.append(labels)
		    self.allSectors.append(sectors)

        ##  Optionally normalize the whole spectrum to zero 
        
        if programOptions.rses_norm_zero:
        
            self.NormalizeToZero(programOptions)

        ##  Include offsets
        
        if offsetList != {}:

            k=0
            
            for nbr_a in programOptions.possibleNa:
                i=0
                for eig in self.allEigenvalues[k]:
                    self.allEigenvalues[k][i] = eig + offsetList[nbr_a]
                    i+=1
                k+=1
        
        ##  Optionally renormalize the spectrum such that it satisfies the
        ##  constraint that sum Exp(-E_i) = 1
        
        ##  Note that high precision arithmetic is used to add the logs
        
        if normalize and not programOptions.rses_norm_zero:

            total = 0

            k=0
            for nbr_a in programOptions.possibleNa:
                for eig in self.allEigenvalues[k]:

                    total += mp.exp(-eig)
                k+=1
            
            total = mp.log(total)
            
            k=0 
            for nbr_a in programOptions.possibleNa:
                i=0
                for eig in self.allEigenvalues[k]:
                    self.allEigenvalues[k][i] = eig + total
                    i+=1
                k+=1
        
            if offsetList != {}:
            
                self.modelParameterList['LinearOffset'] = total   
        
        ##print(self.allSectors)
        ##print(self.allEigenvalues)
        ##raw_input("Press Enter to continue...")        
          
    ############################################################################
    ##	Normalize the fitting data so that the lowest lying entanglement 
    ##  eigenvalue in the lowest sector is set to 0
    ##
    def NormalizeToZero(self,programOptions):
        
        ##  First obtain the lowest lying eigenvalue in sector 0
        ##  and for each na value

        k=0
        
        for nbr_a in programOptions.possibleNa:

            minEigenvalue=100000
            i=0

            for sector in self.allSectors[k]:
        
                if sector == 0:
                    
                    if self.allEigenvalues[k][i]<minEigenvalue:
                    
                        minEigenvalue = self.allEigenvalues[k][i]
                    
                i+=1
            
            ##  Subtract that value from the whole spectrum
        
            i=0
            for eig in self.allEigenvalues[k]:
            
                self.allEigenvalues[k][i] = eig - minEigenvalue
                i+=1
                
            k+=1

        ##raw_input("Press Enter to continue...")

    #############################################################################
    ##  Calculate the sum of exp(-eigenvalue) in order to check the overall
    ##  rses normalization
    ##
    ##  Also make a file containing the final value of the lowest entanglement
    ##  energy, which is the normalization of the state
    ##
    def CalculateNormalizations(self,normSector,programOptions):

        print("\n\tCHECKING SUM OF EXP(-E) AND LOWEST ENTANGLEMENT ENERGY STATE\n")

        ##  First calculate the spectrum for normsector sectors
        ##  and do not normalize it to zero
        
        self.GenerateModelSpectrum(programOptions,normSector)
        
        self.GetSpectrumFromFile(programOptions,normSector,False,self.CalculateOffsets(True,programOptions))

        currTotal = 0.0
        k = 0
        norms = []

        for nbr_a in programOptions.possibleNa:

            counter = 0
            sectorCounter = 0
            currSector = 0
            
            sectorTotals = []
           
            i = 0

            for eig in self.allEigenvalues[k]:

                if i==0:
                    
                    norm = self.allEigenvalues[k][i]
                
                    norms.append(norm)

                    fout=open(programOptions.BuildRsesNormFileName(programOptions.fitting_model,nbr_a),'w')
                
                    fout.write('##  THIS FILE CONTAINS THE MINIMUM ENTANGLEMENT ENERGY VALUE IN SECTOR 0:\n')

                    fout.write(str(norm)+'\n')
            
                    fout.close()

                if self.allSectors[k][i] == currSector:
                
                    currTotal += math.exp(-eig)

                else:
                
                    sectorTotals.append(currTotal)
                    
                    currTotal += math.exp(-eig)

                    currSector -= 1
                    
                i += 1
            
            k += 1

            print("\n\tSector totals of exp(-E):",sectorTotals)
       
        print("\n\tOverall total of exp(-E):",currTotal)
                
        print ("\n\tLowest entanglement energy states:",norms)
        
        print ("\n")
        
        return norms,currTotal
    
    ############################################################################
    ##  Plot the RSES data stored within the class. Note that matplotlib
    ##  and plot_utilities need to be imported before this function can be 
    ##  called
    ##
    def Plot(self,axes,nbr_a,fitSector,isInset,programOptions):
        
        sectorOffset = 0.4
        arrowOffset  = 0.12
        markerSize   = 90
        
        if isInset:
            markerSize   = 80
            arrowOffset  = 0.12
            sectorOffset = 0.4
            
        numberOffsetY = 0.5
        markerSize   = 100
        numberOffsetX = 0
        annotateFontSize = 16
        
        if isInset:
            numberOffsetY = 0.15
            markerSize   = 100
        
        if programOptions.plot_rses:
            markerSize = 140
            numberOffsetX = 0.15
        
        if programOptions.nbr_ll == 2:
            
            numberOffsetY=1.75       ##  1.5 for nu=2
        
        if programOptions.plot_degeneracy:
        
            numberOffsetX = 0.5
            numberOffsetY = 0.01
            annotateFontSize = 6
            
            if programOptions.nbr_ll == 2:
                
                numberOffsetX = 0.25
                annotateFontSize = 6
        
        if programOptions.fitting_model == 'entanglement_energy':
            
            ##fitModelTypeLabel=r"EE MODEL"
            if programOptions.use_perturbation:
            
                ##fitModelTypeLabel=r"Single-Particle Model $+$ Perturbation"
                fitModelTypeLabel=r"Model"
            else:
                ##fitModelTypeLabel=r"Single-Particle Model"
                fitModelTypeLabel=r"Model"
            
        elif programOptions.fitting_model == 'DRR431':
        
            fitModelTypeLabel=r"DRR 4.31"
            
        elif programOptions.fitting_model == 'extended_DRR':
        
            fitModelTypeLabel=r"extended DRR"
            
        ##  Determine which index corresponds to nbr_a
        
        nbrIndex = programOptions.GetPossibleNaIndex(nbr_a)

        if programOptions.highlight_branches and programOptions.nbr_ll>1:
            
            ##  Divide up the plot data so that the different branches can be highlighted

            valList,myColourMap,bounds,myNorm = plot_utilities.DiscreteColourMap(self.allLabels[nbrIndex])
        
            fig = axes.scatter([x+sectorOffset for x in self.allSectors[nbrIndex]],self.allEigenvalues[nbrIndex],s=markerSize,c=valList,marker='_',lw=0.4,cmap=myColourMap,norm=myNorm)
            
            fig = axes.scatter([x+sectorOffset+arrowOffset for x in self.allSectors[nbrIndex]],self.allEigenvalues[nbrIndex],s=15,c=valList,marker='<',lw=0.2,cmap=myColourMap,norm=myNorm,edgecolor='none',label=fitModelTypeLabel)
            
            cbar = axes.get_figure().colorbar(fig,cmap=myColourMap,norm=myNorm,spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i',anchor=(-0.45,0.0))
            
            ##axes.get_figure().text(0.95,1.02,r"$E'_{\mbox{\small Cyclotron}}$",transform=axes.transAxes,**plot_utilities.font)

            if programOptions.nbr_ll == 2:
                
                cbar.set_ticks([x+0.5 for x in range(0,6)])
                
                if programOptions.nbr==60 or programOptions.nbr==62:
                    cbar.ax.set_yticklabels([r'\textbf{13}',r'\textbf{14}',r'\textbf{15}',r'\textbf{16}',r'\textbf{17}',r'\textbf{18}'])
                
                if programOptions.nbr==46 or programOptions.nbr==48:
                    cbar.ax.set_yticklabels([r'\textbf{10}',r'\textbf{11}',r'\textbf{12}',r'\textbf{13}',r'\textbf{14}',r'\textbf{15}'])
                    
                if programOptions.nbr==36 or programOptions.nbr==38:
                    cbar.ax.set_yticklabels([r'\textbf{6}',r'\textbf{7}',r'\textbf{8}',r'\textbf{9}',r'\textbf{10}',r'\textbf{11}'])
            
        else:
        
            axes.scatter([x+sectorOffset for x in self.allSectors[nbrIndex]],self.allEigenvalues[nbrIndex],s=markerSize,c=plot_utilities.colour_list[3],marker='_',lw=0.4)
            axes.scatter([x+sectorOffset+arrowOffset for x in self.allSectors[nbrIndex]],self.allEigenvalues[nbrIndex],s=15,c=plot_utilities.colour_list[3],marker='<',lw=0.2,edgecolor='none',label=fitModelTypeLabel)
    
        if programOptions.nbr_ll == 2:

            axes.get_figure().text(0.62,0.62,r"$\sigma'_{\mbox{\small tot.}}$",transform=axes.transAxes,**plot_utilities.font)

        ##  Highlight the fitted region
            
        if fitSector != programOptions.plot_sector:
            
            axes.axvspan(-fitSector-0.3,0.8, facecolor='black', alpha=0.1,label=r"Fitted Region")
            
        ##  Plot degeneracy number label next to each set of levels

        if programOptions.plot_counting or programOptions.plot_degeneracy:
        
            if programOptions.plot_counting:
                branchGap = 1.8     ##  1.0 for nu=2, 1.8 for 2/3
            elif programOptions.plot_degeneracy:
                branchGap = 0.00001     ##  Set small to test for degeneracies
            
            minCount = 1
            maxSector = 8
            currLza  = 0
            count    = 0
            countList  = []
            sectorList = []
            energyList = []
            innerEnergyList = []
            colourLabels = []

            for i in range(0,len(self.allSectors[nbrIndex])):
            
                x = self.allSectors[nbrIndex][i]
                e = self.allEigenvalues[nbrIndex][i]

                if x<currLza or i == len(self.allSectors[nbrIndex])-1:
                    ##  Skip to next sector
                    if count>minCount:    
                        countList.append(count)
                    else:
                        countList.append('')
                        
                    sectorList.append(currLza)
                    energyList.append(sorted(innerEnergyList)[0])
                    colourLabels.append(label)
                        
                    innerEnergyList = []
                    currLza -= 1
                    count = 0
                    
                    if currLza < -maxSector:
                        break
                
                if len(innerEnergyList)>0:
                    if abs(innerEnergyList[-1]-e)>branchGap:
                        ##  Skip to next branch
                        ##if count <= 100:
                        
                        if count>minCount:    
                            countList.append(count)
                        else:
                            countList.append('')
                            
                        sectorList.append(currLza)
                        energyList.append(sorted(innerEnergyList)[0])
                        colourLabels.append(label)
                        
                        innerEnergyList = []
                        count = 0
                
                count += 1
                innerEnergyList.append(e)
                
                label = self.allLabels[nbrIndex][i]
                
            print("Count of plotted eigenvalues")
            print(countList)

            if programOptions.highlight_branches and programOptions.nbr_ll>1:
                
                for i in range(0,len(countList)):
                
                    ##(-1)**i*
                    axes.annotate(str(countList[i]),xy=(sectorList[i]+numberOffsetX+sectorOffset,energyList[i]-numberOffsetY),size=annotateFontSize,color=plot_utilities.colour_list[int(colourLabels[i]-10)])
                
            else:

                for i in range(0,len(countList)):
                        
                    axes.annotate(str(countList[i]),xy=(sectorList[i]+(-1)**i*numberOffsetX+sectorOffset,energyList[i]-numberOffsetY),size=annotateFontSize)
   
    ################################################################################
    ##  This function calculates the Renyi entropy:
    ##  H_q = 1/(1-q) Log[Sum exp(-q E_i)]
    ##
    ##  Setting q=1 calculates the Von Neumann entropy:
    ##
    ##  S = Sum E_i Exp(-E_i)
    ##
    def CalculateRenyiEntropy(self,q,normSector,fitSector,programOptions):
 
        ##  First calculate the spectrum for normsector sectors
        ##  and do not normalize it to zero
        
        self.GenerateModelSpectrum(programOptions,normSector)
        
        self.GetSpectrumFromFile(programOptions,normSector,False,self.CalculateOffsets(True,programOptions))

        currTotal = 0.0
        k = 0

        for nbr_a in programOptions.possibleNa:

            counter = 0
            sectorCounter = 0
            currSector = 0
            
            sectorTotals = []
           
            i = 0

            for eig in self.allEigenvalues[k]:

                if self.allSectors[k][i] == currSector:
                
                    if q==1:
                        currTotal += eig * mp.exp(-eig)
                    else:
                        currTotal += mp.exp(-q*eig)

                else:
                
                    sectorTotals.append(currTotal)
                    
                    if q==1:
                        currTotal += eig * mp.exp(-eig)
                    else:
                        currTotal += mp.exp(-q*eig)

                    currSector -= 1
                    
                i += 1
            
            k += 1

            ##if q==1:
            ##    print("\n\tSector totals of Von Neumann entropy:",sectorTotals)
            ##else:
            ##    print("\n\tSector totals of "+str(q)+"th Renyi entropy:",(1/(1-q))*mp.log(sectorTotals))
        
        if q==1:
            print("\n\tOverall total of Von Neumann entropy:",currTotal)
        else:
            print("\n\tOverall total of "+str(q)+"th Renyi entropy:",(1/(1-q))*mp.log(currTotal)) 
    
    ################################################################################
    ##  This function reajusts the linear offset parameter such that the spectrum
    ##  is normalized
    ##
    def RenormalizeSpectrum(self,fitSector,programOptions):

        tempExtraNa = programOptions.extra_nbr_a
        
        ##  Fix the sector and extra nbr a values used in the norm calculation
        test_sector = 14
        programOptions.extra_nbr_a = 4
        programOptions.GeneratePossibleNa()
        
        ##  Obtain fitting parameter values for the given fitSector value
        
        self.GetFittingParametersFromFile(programOptions,fitSector)

        tempOption = programOptions.rses_norm_zero

        if programOptions.rses_norm_zero:
            programOptions.rses_norm_zero = False
            print('WARNING - disabling rses_norm_zero in order to calculate entanglement entropy')

        minima,total = self.CalculateNormalizations(test_sector,programOptions)
        
        programOptions.rses_norm_zero = tempOption
        
        #  Update the linear offset such that the spectrum is correctly
        ##  normalized
        
        print ("\tSHIFTING SPECTRUM TO ACHIEVE NORM OF 1")

        if total > 0.0:
            self.modelParameterList['LinearOffset'] += math.log(total)
        
        ##  Reset original parameters
        programOptions.extra_nbr_a = tempExtraNa
        programOptions.GeneratePossibleNa()
                   
################################################################################
##  Define a class to import and contain single particle entanglement energy 
##  model parameters and data
##
class SingleParticleEnergy(FittingModel):
    
    ############################################################################
    ##  This function initializes a set of random fitting parameters, depending
    ##  on the type of model to be fit. See MinFunction for the definitions of
    ##  the fitting models
    ##
    def __init__(self,programOptions):
    
        FittingModel.__init__(self)
        
        if programOptions.wf_type == 'iqhe' and programOptions.nbr_ll == 1:
            self.nbrFittingParameters = 2
        elif programOptions.wf_type == 'laughlin':      
            self.nbrFittingParameters = 3
        else:
            self.nbrFittingParameters = 4
        
        if programOptions.fit_offset:
            self.nbrFittingParameters += 1
        
        if programOptions.use_perturbation:
            self.nbrFittingParameters = 2
            
        self.initialFittingParameters = []
        
        random.seed()
        
        for i in range(0,self.nbrFittingParameters):
            self.initialFittingParameters.append(random.random()*0.2)

        ##self.initialFittingParameters[0]=-0.243
        ##self.initialFittingParameters[1]=0.00130
        ##self.initialFittingParameters[2]=-0.00214
        
        self.initialFittingParameters[0]=0.394581859561
        self.initialFittingParameters[1]=0.9279121962
        self.initialFittingParameters[2]=3.37095904351
        self.initialFittingParameters[3]=0.642918301597

        ##  Initialize the full list of model parameters

        self.modelParameterList = {'LinearOffset':0,'QuadraticOffset':0,'a_0':0,'a_1':0,'a_2':0,'a_3':0,'a_5':0,'b_0':0,'b_1':0, \
        'b_2':0,'b_3':0,'V1':0,'V3':0,'V5':0,'V7':0,'c_0':0,'c_1':0,'c_2':0,'c_3':0,'c_5':0,'d_0':0,'d_1':0,'d_2':0,'d_3':0}

    ############################################################################
    ##	This function calls the single particle entanglement energy calculation
    ##  program. Fitting parameters are passed to the program via a text file
    ##
    def GenerateModelSpectrum(self,programOptions,fitSector):

        for nbr_a in programOptions.possibleNa:

            ##  First generate parameter.dat file to contain current
            ##  fitting parameters

            ## epsilon_{m,sigma} = a_0 + a_1 m + a_2 m^2 + a_3 m^3  + a_5 m^5
            ## + b_0 sigma + b_1 sigma m + b_2 sigma m^2 + b_3 sigma m^3

            fileName = 'parameters_'+str(programOptions.fitting_model)+'_model_n_'+str(programOptions.nbr)+'_na_'+str(nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'.dat'
                
            fout=open(programOptions.BuildOutDirectoryName()+fileName,'w')

            fout.write('##  THIS FILE CONTAINS ENTANGLEMENT ENERGY FITTING PARAMETER DATA FOR:\n')

            if programOptions.use_perturbation:

                fout.write('##  Perturbed entanglement energy model\n')
                fout.write('##  a_0, a_1,a_2,a_3,a_5, b_0,b_1,b_2,b_3, Perturbation parameter, c_0,c_1,c_2,c_3,c_5,d_0,d_1,d_2,d_3\n\n')
            
            else:

                fout.write('##  Entanglement energy model\n')
                fout.write('##  a_0, a_1,a_2,a_3,a_5, b_0,b_1,b_2,b_3,c_0,c_1,c_2,c_3,c_5,d_0,d_1,d_2,d_3\n\n')

            fout.write(str(self.modelParameterList['a_0'])+'\n')
            fout.write(str(self.modelParameterList['a_1'])+'\n')
            fout.write(str(self.modelParameterList['a_2'])+'\n')
            fout.write(str(self.modelParameterList['a_3'])+'\n')
            fout.write(str(self.modelParameterList['a_5'])+'\n')
            fout.write(str(self.modelParameterList['b_0'])+'\n')
            fout.write(str(self.modelParameterList['b_1'])+'\n')
            fout.write(str(self.modelParameterList['b_2'])+'\n')
            fout.write(str(self.modelParameterList['b_3'])+'\n')     
            fout.write(str(self.modelParameterList['V1'])+'\n')
            fout.write(str(self.modelParameterList['V3'])+'\n')
            fout.write(str(self.modelParameterList['V5'])+'\n')
            fout.write(str(self.modelParameterList['V7'])+'\n')
            fout.write(str(self.modelParameterList['c_0'])+'\n')
            fout.write(str(self.modelParameterList['c_1'])+'\n')
            fout.write(str(self.modelParameterList['c_2'])+'\n')
            fout.write(str(self.modelParameterList['c_3'])+'\n')
            fout.write(str(self.modelParameterList['c_5'])+'\n')
            fout.write(str(self.modelParameterList['d_0'])+'\n')
            fout.write(str(self.modelParameterList['d_1'])+'\n')
            fout.write(str(self.modelParameterList['d_2'])+'\n')
            fout.write(str(self.modelParameterList['d_3'])+'\n')

            fout.close()

            ##	Call a C++ program to evaluate the spectrum for these parameters

            ##print(fitSector)
            
            ##print(str(programOptions.bin_path)+'//'+str(programOptions.entanglement_energy_bin)+' -v 0 -n '+str(programOptions.nbr)+' -a '+str(nbr_a)+' --nbr-ll '+str(programOptions.nbr_ll)+' --path '+programOptions.BuildOutDirectoryName()+' --max-sector '+str(fitSector)+' --min-sector '+str(0)+' --use-perturbation '+str(int(programOptions.use_perturbation)))

            os.system(str(programOptions.bin_path)+'//'+str(programOptions.entanglement_energy_bin)+' -v 0 -n '+str(programOptions.nbr)+' -a '+str(nbr_a)+' --nbr-ll '+str(programOptions.nbr_ll)+' --in-path '+programOptions.BuildOutDirectoryName()+' --out-path '+programOptions.BuildOutDirectoryName()+' --max-sector '+str(fitSector)+' --min-sector '+str(0)+' --use-perturbation '+str(int(programOptions.use_perturbation)))
            
            ##raw_input("Press Enter to continue...")
    
    ############################################################################
    ##  CalculateOffsets method. Note that the offset function must be 
    ##  normalized so that sum exp(-E_i,NA) = 1.
    ##  
    ##  To specify that, we will set the quadratic coefficient in the
    ##  fitting function, then impose the summation constraint to 
    ##  set the linear coefficient once the specturm is calculated
    ##
    def CalculateOffsets(self,includeLinear,programOptions):

        offsetList = {}

        if programOptions.fit_offset:
            
            naMin = int(programOptions.nbr/2.0)

            for na in programOptions.possibleNa:

                offset = 0

                if includeLinear:
                    
                    offset += self.modelParameterList['LinearOffset']
                
                offset += self.modelParameterList['QuadraticOffset']*(na-naMin)**2
                    
                offsetList[na] = offset
                
                #print(offsetList)
                
        return offsetList
        
    ############################################################################
    ##  This function calls external programs to generate the fitting model 
    ##  spectrum, then it determines the current value of the fitting 
    ##  parameter.
    ##
    ##  The particular models that are used for fitting are defined here
    ##
    def MinFunction(self,rawFittingParameters,rsesData,programOptions,fitSector):
        
        ## epsilon_{m,sigma} = a_0 + a_1 m + a_2 m^2 + a_3 m^3  + a_5 m^5
        ## + b_0 sigma + b_1 sigma m + b_2 sigma m^2 + b_3 sigma m^3

        message = '\tcurrent iteration: '

        for val in rawFittingParameters:
        
            message += str(val)+'\t'

        print (message)

        i = 0
        
        if programOptions.fit_offset:
            
            ##  The linear offset is set by the overall normalization constraint
            
            self.modelParameterList['QuadraticOffset'] = rawFittingParameters[0]
            
            i += 1
        
        self.modelParameterList['c_0'] = 0##rawFittingParameters[i]
        self.modelParameterList['c_1'] = 0##rawFittingParameters[i+1]
        self.modelParameterList['c_2'] = 0##rawFittingParameters[i+2]
        self.modelParameterList['c_3'] = 0
        self.modelParameterList['c_5'] = 0
        
        ##  use d_0 fit for the Jain state
        self.modelParameterList['d_0'] = rawFittingParameters[i]
        self.modelParameterList['d_1'] = 0##rawFittingParameters[i+2]
        self.modelParameterList['d_2'] = 0
        self.modelParameterList['d_3'] = 0
        
        if programOptions.wf_type == 'iqhe' and programOptions.nbr_ll == 1  and not programOptions.double_energy:

            self.modelParameterList['a_0'] = 0
            self.modelParameterList['a_1'] = rawFittingParameters[i+1]
            self.modelParameterList['a_2'] = 0
            self.modelParameterList['a_3'] = rawFittingParameters[i+2]
            self.modelParameterList['a_5'] = 0
            
            self.modelParameterList['b_0'] = 0
            self.modelParameterList['b_1'] = 0
            self.modelParameterList['b_2'] = 0
            self.modelParameterList['b_3'] = 0
            
            self.modelParameterList['V1'] = 0
            self.modelParameterList['V3'] = 0
            self.modelParameterList['V5'] = 0
            self.modelParameterList['V7'] = 0

        elif programOptions.wf_type == 'laughlin':

            self.modelParameterList['a_0'] = 0
            self.modelParameterList['a_1'] = rawFittingParameters[i]
            self.modelParameterList['a_2'] = rawFittingParameters[i+1]
            self.modelParameterList['a_3'] = rawFittingParameters[i+2]
            self.modelParameterList['a_5'] = 0
            
            self.modelParameterList['b_0'] = 0
            self.modelParameterList['b_1'] = 0
            self.modelParameterList['b_2'] = 0
            self.modelParameterList['b_3'] = 0
            
            if programOptions.use_perturbation:
                self.modelParameterList['V1'] = rawFittingParameters[i]
                self.modelParameterList['V3'] = rawFittingParameters[i+1]
                self.modelParameterList['V5'] = 0##rawFittingParameters[i+2]
                self.modelParameterList['V7'] = 0
            else:
                self.modelParameterList['V1'] = 0
                self.modelParameterList['V3'] = 0
                self.modelParameterList['V5'] = 0
                self.modelParameterList['V7'] = 0    
        else:

            self.modelParameterList['a_0'] = 0##rawFittingParameters[i]
            self.modelParameterList['a_1'] = rawFittingParameters[i+1]
            self.modelParameterList['a_2'] = 0##rawFittingParameters[i+5]
            self.modelParameterList['a_3'] = 0##rawFittingParameters[i+5]
            self.modelParameterList['a_5'] = 0##rawFittingParameters[i+7]
            
            self.modelParameterList['b_0'] = rawFittingParameters[i+2]
            self.modelParameterList['b_1'] = rawFittingParameters[i+3]
            self.modelParameterList['b_2'] = 0##rawFittingParameters[i+2]
            self.modelParameterList['b_3'] = 0##rawFittingParameters[i+5]
            
            if programOptions.use_perturbation:
                self.modelParameterList['V1'] = rawFittingParameters[i]
                self.modelParameterList['V3'] = rawFittingParameters[i+1]
                self.modelParameterList['V5'] = 0##rawFittingParameters[i+2]
                self.modelParameterList['V7'] = 0##rawFittingParameters[i+3]
            else:
                self.modelParameterList['V1'] = 0
                self.modelParameterList['V3'] = 0
                self.modelParameterList['V5'] = 0
                self.modelParameterList['V7'] = 0

        self.GenerateModelSpectrum(programOptions,fitSector)
          
        self.GetSpectrumFromFile(programOptions,fitSector,True,self.CalculateOffsets(False,programOptions))

        ##raw_input("Press Enter to continue...")
        
        ##	Determine least-squares fitting parameter
        
        return self.FittingParameter(rsesData,programOptions)

    ############################################################################
    ##  Function to store final fitting parameter values in a file
    ##
    def StoreFinalFitParameters(self,fittingResult,programOptions,fitSector):

        print ("Final entanglement energy Parameters:")
        	
        print (self.modelParameterList)

        fout=open(programOptions.BuildOutDirectoryName()+'parameters_entanglement_energy_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'_fit_'+str(fitSector)+'.dat','w')
        
        fout.write('##  THIS FILE CONTAINS ENTANGLEMENT ENERGY FITTING PARAMETER DATA FOR:\n')
        
        if programOptions.use_perturbation:

            fout.write('##  Perturbed entanglement energy model\n')
            fout.write('##  Stored: epsilon_{m,sigma} = a_0 + a_1 m + a_2 m^2 + a_3 m^3  + a_5 m^5 + b_0 sigma + b_1 sigma m + b_2 sigma m^2 + b_3 sigma m^3, 4 Perturbation parameters then c_0 to c_5 and d_0 to d_3\n\n')
        else: 

            fout.write('##  Entanglement energy model\n')
            fout.write('##  Stored: epsilon_{m,sigma} = epsilon_{m,sigma} = a_0 + a_1 m + a_2 m^2 + a_3 m^3  + a_5 m^5 + b_0 sigma + b_1 sigma m + b_2 sigma m^2 + b_3 sigma m^3 then c_0 to c_5 and d_0 to d_3\n\n')

        fout.write(str(self.modelParameterList['a_0'])+'\n')
        fout.write(str(self.modelParameterList['a_1'])+'\n')
        fout.write(str(self.modelParameterList['a_2'])+'\n')
        fout.write(str(self.modelParameterList['a_3'])+'\n')
        fout.write(str(self.modelParameterList['a_5'])+'\n')
        fout.write(str(self.modelParameterList['b_0'])+'\n')
        fout.write(str(self.modelParameterList['b_1'])+'\n')
        fout.write(str(self.modelParameterList['b_2'])+'\n')
        fout.write(str(self.modelParameterList['b_3'])+'\n')  
        fout.write(str(self.modelParameterList['V1'])+'\n')
        fout.write(str(self.modelParameterList['V3'])+'\n')
        fout.write(str(self.modelParameterList['V5'])+'\n')
        fout.write(str(self.modelParameterList['V7'])+'\n')
        fout.write(str(self.modelParameterList['c_0'])+'\n')
        fout.write(str(self.modelParameterList['c_1'])+'\n')
        fout.write(str(self.modelParameterList['c_2'])+'\n')
        fout.write(str(self.modelParameterList['c_3'])+'\n')
        fout.write(str(self.modelParameterList['c_5'])+'\n')
        fout.write(str(self.modelParameterList['d_0'])+'\n')
        fout.write(str(self.modelParameterList['d_1'])+'\n')
        fout.write(str(self.modelParameterList['d_2'])+'\n')
        fout.write(str(self.modelParameterList['d_3'])+'\n')

        if programOptions.fit_offset:
            fout.write('\n## OFFSETS:\n\n')
            fout.write(str(self.modelParameterList['LinearOffset'])+'\n')
            fout.write(str(self.modelParameterList['QuadraticOffset'])+'\n')

        fout.write('\n## FITTING RESULT:\n\n'+str(fittingResult)+'\n')
        
        fout.close()
    
    ############################################################################
    ##  Function to obtain the previously stored final fitting parameter values
    ##
    def GetFittingParametersFromFile(self,programOptions,fitSector):

        fileName = programOptions.BuildOutDirectoryName()+'parameters_entanglement_energy_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'_fit_'+str(fitSector)+'.dat'
        
        print ("\n\tReading in entanglement energy model parameters from file "+str(fileName))
        
        try:       
            fin=open(fileName)
            fin.close()
        except: 
            print("ERROR: file "+str(fileName)+" not found")
            quit()
            
        fin=open(fileName)
            
        ##  Skip comment lines
        fin.readline()
        fin.readline()
        fin.readline()
        fin.readline()
        
        ##  Read in parameter values
        self.modelParameterList['a_0'] = float(fin.readline())
        self.modelParameterList['a_1'] = float(fin.readline())
        self.modelParameterList['a_2'] = float(fin.readline())
        self.modelParameterList['a_3'] = float(fin.readline())
        self.modelParameterList['a_5'] = float(fin.readline())
        self.modelParameterList['b_0'] = float(fin.readline())
        self.modelParameterList['b_1'] = float(fin.readline())
        self.modelParameterList['b_2'] = float(fin.readline())
        self.modelParameterList['b_3'] = float(fin.readline())
        self.modelParameterList['V1'] = float(fin.readline())
        self.modelParameterList['V3'] = float(fin.readline())
        self.modelParameterList['V5'] = float(fin.readline())
        self.modelParameterList['V7'] = float(fin.readline())
        self.modelParameterList['c_0'] = float(fin.readline())
        self.modelParameterList['c_1'] = float(fin.readline())
        self.modelParameterList['c_2'] = float(fin.readline())
        self.modelParameterList['c_3'] = float(fin.readline())
        self.modelParameterList['c_5'] = float(fin.readline())
        self.modelParameterList['d_0'] = float(fin.readline())
        self.modelParameterList['d_1'] = float(fin.readline())
        self.modelParameterList['d_2'] = float(fin.readline())
        self.modelParameterList['d_3'] = float(fin.readline())

        self.modelParameterList['LinearOffset']    = 0
        self.modelParameterList['QuadraticOffset'] = 0

        if programOptions.fit_offset:

            ##  Skip comment lines
            fin.readline()
            fin.readline()
            fin.readline()

            self.modelParameterList['LinearOffset']    = float(fin.readline())
            self.modelParameterList['QuadraticOffset'] = float(fin.readline())

        fin.close()
        
        print ("\n\tGot parameters:")

        print (self.modelParameterList)
        
    #############################################################################
    ##  A function to clear up any files produced at intermediate steps in
    ##  the fitting algorithm operation. These are input and output files
    ##  from the entanglement energy program
    ##
    def ClearWorkingFiles(self,programOptions):

        parametersFileName = ''

        for nbr_a in programOptions.possibleNa:

            parametersFileName += programOptions.BuildOutDirectoryName()+'parameters_entanglement_energy_model_n_'+str(programOptions.nbr)+'_na_'+str(nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'.dat '
        
        outputFileName = programOptions.BuildOutDirectoryName()+'eigenvalues_entanglement_energy_model*'
        
        os.system('rm '+parametersFileName+' '+outputFileName)

################################################################################
##  Define a class to import and contain DRR431 entanglement Hamiltonian data
##
class DRR431Hamiltonian(FittingModel):

    ############################################################################
    ##  This function initializes a set of random fitting parameters, depending
    ##  on the type of model to be fit. See MinFunction for the definitions of
    ##  the fitting models
    ##
    def __init__(self,programOptions):
    
        FittingModel.__init__(self)

        self.nbrFittingParameters = 4

        self.initialFittingParameters = []
        
        for i in range(0,self.nbrFittingParameters):
            self.initialFittingParameters.append(random.random())
        
        ##  Initialize the full list of model parameters    
            
        self.modelParameterList = {'LinearOffset':0,'Alpha':0,'Beta':0,'Gamma':0}
        
    ############################################################################
    ##	A function to set the diagonalization routines running to evaluate
    ##	the entanglement Hamiltonian spectrum up to a specified sector. 
    ##	Using the Dubail et al Eq 4.31 expression for the 
    ##	entanglement Hamiltonian
    ##
    def GenerateModelSpectrum(self,programOptions,fitSector):

        for nbr_a in programOptions.possibleNa:

            ##	First generate parameter.dat file to contain fitting parameters

            fileName = 'parameters_'+str(programOptions.fitting_model)+'_model_n_'+str(programOptions.nbr)+'_na_'+str(nbr_a)+'.dat'
            
            fout=open(programOptions.BuildOutDirectoryName()+fileName,'w')

            fout.write('##  THIS FILE CONTAINS ENTANGLEMENT HAMILTONIAN FITTING PARAMETER DATA FOR:\n')
            fout.write('##  DUBAIL ET AL. EQ 3.41\n')
            fout.write('##  ALPHA, BETA, GAMMA \n\n')

            fout.write(str(self.modelParameterList['Alpha'])+'\n')
            fout.write(str(self.modelParameterList['Beta'])+'\n')
            fout.write(str(self.modelParameterList['Gamma'])+'\n')

            fout.close()

            ##	Call a C++ program to evaluate the spectrum for these parameters

            os.system(str(programOptions.bin_path)+'//'+str(programOptions.entanglement_hamiltonian_bin)+' --blocks 0 --hamiltonian '+str(programOptions.fitting_model)+' -v 0 -n '+str(programOptions.nbr)+' -a '+str(nbr_a)+' --path '+programOptions.BuildOutDirectoryName()+' --max-sector '+str(fitSector)+' --min-sector '+str(0)+' --power '+str(programOptions.jastrow+1)+' --parameters-file '+fileName)

    ############################################################################
    ##  CalculateOffsets method 
    ##
    def CalculateOffsets(self,includeLinear,programOptions):
    
        offsetList = {}

        for na in programOptions.possibleNa:
            if includeLinear:
                offsetList[na] = self.modelParameterList['LinearOffset']
        
        return offsetList

    ############################################################################
    ##  This function calls external programs to generate the fitting model 
    ##  spectrum, then it determines the current value of the fitting 
    ##  parameter.
    ##
    ##  The particular models that are used for fitting are defined here
    ##
    def MinFunction(self,rawFittingParameters,rsesData,programOptions,fitSector):
        
        ##	eigenvalues, maxSector and baseName are set globally

        ##	Generate entanglement Hamiltonian spectrum 

        message = '\tcurrent iteration: '

        for val in rawFittingParameters:
        
            message += str(val)+'\t'

        print (message)

        self.modelParameterList['LinearOffset'] = rawFittingParameters[0]
        self.modelParameterList['Alpha']        = rawFittingParameters[1]
        self.modelParameterList['Beta']         = rawFittingParameters[2]
        self.modelParameterList['Gamma']        = rawFittingParameters[3]

        self.GenerateModelSpectrum(programOptions,fitSector)

        self.GetSpectrumFromFile(programOptions,fitSector,False,self.CalculateOffsets(True,programOptions))

        ##raw_input("Press Enter to continue...")
        
        ##	Determine least-squares fitting parameter
        
        return self.FittingParameter(rsesData,programOptions)

    ############################################################################
    ##  Function to store final fitting parameter values in a file
    ##
    def StoreFinalFitParameters(self,fittingResult,programOptions,fitSector):

        print ("Final DRR431 Parameters:")
        	
        print (self.modelParameterList)
        
        fout=open(programOptions.BuildOutDirectoryName()+'parameters_DRR431_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_fit_'+str(fitSector)+'.dat','w')
        
        fout.write('##  THIS FILE CONTAINS ENTANGLEMENT HAMILTONIAN FITTING PARAMETER DATA FOR:\n')
        fout.write('##  DUBAIL ET AL. EQ 3.41. \n')
        fout.write('##  OFFSET, ALPHA, BETA, GAMMA \n\n')
        
        fout.write(str(self.modelParameterList['LinearOffset'])+'\n')
        fout.write(str(self.modelParameterList['Alpha'])+'\n')
        fout.write(str(self.modelParameterList['Beta'])+'\n')
        fout.write(str(self.modelParameterList['Gamma'])+'\n')
        
        fout.write('\n## FITTING RESULT:\n\n'+str(fittingResult)+'\n')
        
        fout.close()

    ############################################################################
    ##  Function to obtain the previously stored final fitting parameter values
    ##
    def GetFittingParametersFromFile(self,programOptions,fitSector):
        
        fileName = programOptions.BuildOutDirectoryName()+'parameters_DRR431_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_fit_'+str(fitSector)+'.dat'
        
        print ("\n\tReading in DRR431 model parameters from file "+str(fileName))
        
        try:       
            fin=open(fileName)
            fin.close()
        except: 
            print("ERROR: file "+str(fileName)+" not found")
            quit()
            
        fin=open(fileName) 
        ##  Skip comment lines
        
        fin.readline()
        fin.readline()
        fin.readline()
        fin.readline()
        
        ##  Read in parameter values
        
        self.modelParameterList['LinearOffset'] = float(fin.readline())
        self.modelParameterList['Alpha']  = float(fin.readline())
        self.modelParameterList['Beta']   = float(fin.readline())
        self.modelParameterList['Gamma']  = float(fin.readline())

        fin.close()
        
        print ("\n\tGot parameters:")

        print (self.modelParameterList)
       
    #############################################################################
    ##  A function to clear up any files produced at intermediate steps in
    ##  the fitting algorithm operation. These are input and output files
    ##  from the entanglement energy program
    ##
    def ClearWorkingFiles(self,programOptions):
    
        parametersFileName = programOptions.BuildOutDirectoryName()+'parameters_DRR431_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'.dat'
            
        outputFileName = programOptions.BuildOutDirectoryName()+'eigenvalues_DRR431_model*'

        os.system('rm '+parametersFileName+' '+outputFileName)

################################################################################
##  Define a class to import and contain the extended DRR entanglement 
##  Hamiltonian data
##
class ExtendedDRRHamiltonian(FittingModel):

    ############################################################################
    ##  This function initializes a set of random fitting parameters, depending
    ##  on the type of model to be fit. See MinFunction for the definitions of
    ##  the fitting models
    ##
    def __init__(self,programOptions):
    
        FittingModel.__init__(self)

        if programOptions.nbr_ll == 3:
            self.nbrFittingParameters = 7
        elif programOptions.nbr_ll == 2:
            self.nbrFittingParameters = 6
            
        self.initialFittingParameters = []
        
        for i in range(0,self.nbrFittingParameters):
            self.initialFittingParameters.append(random.random())
        
        ##  Initialize the full list of model parameters    
            
        self.modelParameterList = {'LinearOffset':0,'Alpha0':0,'Beta0':0,'Gamma0':0,'Alpha1':0,'Beta1':0,'Gamma1':0,'Alpha2':0,'Beta2':0,'Gamma2':0,'Cyclotron0':0,'Cyclotron1':0}
        
    ############################################################################
    ##	A function to set the diagonalization routines running to evaluate
    ##	the entanglement Hamiltonian spectrum up to a specified sector. 
    ##	Using the extended DRR model expression for the entanglement Hamiltonian
    ##
    def GenerateModelSpectrum(self,programOptions,fitSector):

        for nbr_a in programOptions.possibleNa:

            ##    First generate parameter.dat file to contain fitting parameters
            
            fileName = 'parameters_'+str(programOptions.fitting_model)+'_model_n_'+str(programOptions.nbr)+'_na_'+str(nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'.dat'
            
            fout=open(programOptions.BuildOutDirectoryName()+fileName,'w')

            fout.write('##  THIS FILE CONTAINS ENTANGLEMENT HAMILTONIAN FITTING PARAMETER DATA FOR:\n')
            fout.write('##  COMPOSITE FERMION MODEL WITH '+str(programOptions.nbr_ll)+' CF LLS\n')
            fout.write('##  ALPHA[0], BETA[0], GAMMA[0], ALPHA[1], BETA[1], GAMMA[1],...,CYCLOTRON ENERGY[0],CYCLOTRON ENERGY[1]  \n\n')

            fout.write(str(self.modelParameterList['Alpha0'])+'\n')
            fout.write(str(self.modelParameterList['Beta0'])+'\n')
            fout.write(str(self.modelParameterList['Gamma0'])+'\n')
            if programOptions.nbr_ll>=2:
                fout.write(str(self.modelParameterList['Alpha1'])+'\n')
                fout.write(str(self.modelParameterList['Beta1'])+'\n')
                fout.write(str(self.modelParameterList['Gamma1'])+'\n')
            if programOptions.nbr_ll>=3:
                fout.write(str(self.modelParameterList['Alpha2'])+'\n')
                fout.write(str(self.modelParameterList['Beta2'])+'\n')
                fout.write(str(self.modelParameterList['Gamma2'])+'\n')
            if programOptions.nbr_ll>=2:    
                fout.write(str(self.modelParameterList['Cyclotron0'])+'\n')
            if programOptions.nbr_ll>=3:
                fout.write(str(self.modelParameterList['Cyclotron1'])+'\n')
            
            fout.close()
            
            ##    Call a C++ program to evaluate the spectrum for these parameters

            if programOptions.statistics == 'bosons':
                
                statisticsArg = 1
            else:
                statisticsArg = -1

            os.system(str(programOptions.bin_path)+'//'+str(programOptions.entanglement_hamiltonian_bin)+' --hamiltonian '+str(programOptions.fitting_model)+' -v 0 -n '+str(programOptions.nbr)+' -a '+str(nbr_a)+' -c '+str(programOptions.nbr_ll)+' --path '+programOptions.BuildOutDirectoryName()+' --blocks '+str(programOptions.highlight_branches)+" --statistics "+str(statisticsArg)+' --max-sector '+str(fitSector)+' --min-sector '+str(0)+' --parameters-file '+fileName)  

    ############################################################################
    ##  CalculateOffsets method 
    ##
    def CalculateOffsets(self,includeLinear,programOptions):
    
        offsetList = {}

        for na in programOptions.possibleNa:
            if includeLinear:
                offsetList[na] = self.modelParameterList['LinearOffset']
        
        return offsetList

    ############################################################################
    ##  This function calls external programs to generate the fitting model 
    ##  spectrum, then it determines the current value of the fitting 
    ##  parameter.
    ##
    ##  The particular models that are used for fitting are defined here
    ##
    def MinFunction(self,rawFittingParameters,rsesData,programOptions,fitSector):

        message = '\tcurrent iteration: '

        for val in rawFittingParameters:
        
            message += str(val)+'\t'

        print (message)

        if programOptions.nbr_ll == 3:

            self.modelParameterList['LinearOffset']     = rawFittingParameters[0]
            self.modelParameterList['Cyclotron0'] = rawFittingParameters[1]
            self.modelParameterList['Cyclotron1'] = 0
            self.modelParameterList['Alpha0']     = rawFittingParameters[2]
            self.modelParameterList['Alpha1']     = rawFittingParameters[3]
            self.modelParameterList['Alpha2']     = rawFittingParameters[4]
            self.modelParameterList['Beta0']      = rawFittingParameters[5]
            self.modelParameterList['Beta1']      = rawFittingParameters[5]
            self.modelParameterList['Beta2']      = rawFittingParameters[5]
            self.modelParameterList['Gamma0']     = rawFittingParameters[6]
            self.modelParameterList['Gamma1']     = rawFittingParameters[6]
            self.modelParameterList['Gamma2']     = rawFittingParameters[6]
	
        elif programOptions.nbr_ll == 2:

            self.modelParameterList['LinearOffset']     = rawFittingParameters[0]
            self.modelParameterList['Cyclotron0'] = rawFittingParameters[1]
            self.modelParameterList['Cyclotron1'] = 0
            self.modelParameterList['Alpha0']     = rawFittingParameters[2]
            self.modelParameterList['Alpha1']     = rawFittingParameters[3]
            self.modelParameterList['Alpha2']     = 0
            self.modelParameterList['Beta0']      = rawFittingParameters[4]
            self.modelParameterList['Beta1']      = rawFittingParameters[4]
            self.modelParameterList['Beta2']      = 0
            self.modelParameterList['Gamma0']     = rawFittingParameters[5]
            self.modelParameterList['Gamma1']     = rawFittingParameters[5]
            self.modelParameterList['Gamma2']     = 0

        self.GenerateModelSpectrum(programOptions,fitSector)

        self.GetSpectrumFromFile(programOptions,fitSector,True,self.CalculateOffsets(True,programOptions))

        ##raw_input("Press Enter to continue...")
        
        ##	Determine least-squares fitting parameter
        
        return self.FittingParameter(rsesData,programOptions)
        
    ############################################################################
    ##  Function to store final fitting parameter values in a file
    ##
    def StoreFinalFitParameters(self,fittingResult,programOptions,fitSector):

        print ("Final extended DRR model Parameters:")

        print (self.modelParameterList)
        
        fout=open(programOptions.BuildOutDirectoryName()+'parameters_extended_DRR_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'_fit_'+str(fitSector)+'.dat','w')
        
        fout.write('##  THIS FILE CONTAINS ENTANGLEMENT HAMILTONIAN FITTING PARAMETER DATA FOR:\n')
        fout.write('##  COMPOSITE FERMION MODEL WITH '+str(programOptions.nbr_ll)+' CF LLS.\n')
        fout.write('##  OFFSET, ALPHA[0], BETA[0], GAMMA[0], ALPHA[1], BETA[1], GAMMA[1],...,CYCLOTRON ENERGY[0],CYCLOTRON ENERGY[1],  \n\n')
        
        fout.write(str(self.modelParameterList['LinearOffset'])+'\n')
        fout.write(str(self.modelParameterList['Alpha0'])+'\n')
        fout.write(str(self.modelParameterList['Beta0'])+'\n')
        fout.write(str(self.modelParameterList['Gamma0'])+'\n')
        if programOptions.nbr_ll>=2:
            fout.write(str(self.modelParameterList['Alpha1'])+'\n')
            fout.write(str(self.modelParameterList['Beta1'])+'\n')
            fout.write(str(self.modelParameterList['Gamma1'])+'\n')
        if programOptions.nbr_ll>=3:
            fout.write(str(self.modelParameterList['Alpha2'])+'\n')
            fout.write(str(self.modelParameterList['Beta2'])+'\n')
            fout.write(str(self.modelParameterList['Gamma2'])+'\n')
        if programOptions.nbr_ll>=2:
            fout.write(str(self.modelParameterList['Cyclotron0'])+'\n')
        if programOptions.nbr_ll>=3:
            fout.write(str(self.modelParameterList['Cyclotron1'])+'\n')

        fout.write('\n## FITTING RESULT:\n\n'+str(fittingResult)+'\n')
        
        fout.close()
    
    ############################################################################
    ##  Function to obtain the previously stored final fitting parameter values
    ##
    def GetFittingParametersFromFile(self,programOptions,fitSector):
 
        fileName = programOptions.BuildOutDirectoryName()+'parameters_extended_DRR_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'_fit_'+str(fitSector)+'.dat'
        
        print ("\n\tReading in extended DRR model parameters from file "+str(fileName))
        
        try:       
            fin=open(fileName)
            fin.close()
        except: 
            print("ERROR: file "+str(fileName)+" not found")
            quit()
        
        fin=open(fileName)
            
        ##  Skip comment lines
        
        fin.readline()
        fin.readline()
        fin.readline()
        fin.readline()
        
        ##  Read in parameter values
        
        self.modelParameterList['LinearOffset']     = float(fin.readline())
        self.modelParameterList['Alpha0']     = float(fin.readline())
        self.modelParameterList['Beta0']      = float(fin.readline())
        self.modelParameterList['Gamma0']     = float(fin.readline())
        if programOptions.nbr_ll>=2:
            self.modelParameterList['Alpha1']     = float(fin.readline())
            self.modelParameterList['Beta1']      = float(fin.readline())
            self.modelParameterList['Gamma1']     = float(fin.readline())
        if programOptions.nbr_ll>=3:
            self.modelParameterList['Alpha2']     = float(fin.readline())
            self.modelParameterList['Beta2']      = float(fin.readline())
            self.modelParameterList['Gamma2']     = float(fin.readline())
        if programOptions.nbr_ll>=2:
            self.modelParameterList['Cyclotron0'] = float(fin.readline())
        if programOptions.nbr_ll>=3:
            self.modelParameterList['Cyclotron1'] = float(fin.readline())

        fin.close()
        
        print ("\n\tGot parameters:")

        print (self.modelParameterList)
        
    #############################################################################
    ##  A function to clear up any files produced at intermediate steps in
    ##  the fitting algorithm operation. These are input and output files
    ##  from the entanglement energy program
    ##
    def ClearWorkingFiles(self,programOptions):
    
        if programOptions.nbr_ll ==1:
        
            parametersFileName = programOptions.BuildOutDirectoryName()+'parameters_extended_DRR_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'.dat'

        else:
            parametersFileName = programOptions.BuildOutDirectoryName()+'parameters_extended_DRR_model_n_'+str(programOptions.nbr)+'_na_'+str(programOptions.nbr_a)+'_ll_'+str(programOptions.nbr_ll)+'.dat'
            
        outputFileName = programOptions.BuildOutDirectoryName()+'eigenvalues_extended_DRR_model*'    

        os.system('rm '+parametersFileName+' '+outputFileName)
