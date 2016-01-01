#! /usr/bin/env python2.7
################################################################################
##                                                                            ##
##						  Author: Simon C. Davenport						  ##
##                                                                            ##
##							Last Modified: 04/11/2014						  ##
##                                                                            ##
##                                                                            ##
##		This file contains a class to handle reading in and processing        ##
##      RSES data calculated with the ewf, reduced density matrix             ##
##      numerically exact methods                                             ##
##                                                                            ##
################################################################################

################################################################################
################################################################################
##  Import modules
import os
import math
import numpy as np
from utilities import plot_utilities

################################################################################
################################################################################
##  Define a base class for importing rses data of various types
##
class RsesData():

    ############################################################################
    ##  Base class constructor
    ##
    def __init__(self):
        self.allEigenvalues = []
        self.allSectors = []
        self.allFitFunctionWeights = []
    
    ############################################################################
    ##  A function to read in the RSES data from a file. The base class
    ##  contains a default implementation
    ##
    def GetSpectrumFromFile(self,programOptions,maxSector):
        
        self.allEigenvalues = []
        self.allSectors = []

        for nbr_a in programOptions.possibleNa:

            if os.path.isfile(programOptions.BuildRsesFileName(0,nbr_a)):

                print ('\n\tReading in RSES data from '+programOptions.BuildRsesFileName(0,nbr_a)+'\n')

                sectors=[]
                eigenvalues=[]
                sectorEigs=[]
                currSector=0
                
                fin=open(programOptions.BuildRsesFileName(0,nbr_a))
                
                for line in fin:
                
                    line = line.split()

                    if line[1]!="nan":
                    
                        sector     = -int(line[0])
                        eigenvalue = float(line[1])
                    
                        if -sector > currSector:

                            ##  Update the sector counter
                        
                            currSector = -sector
                            
                            ##  Sort the eigenvalues for each sector in descending order

                            sectorEigs.sort()
                            
                            ##  Append to the full list of eigenvalues
                            
                            eigenvalues += sectorEigs
                            
                            sectorEigs = []

                        if currSector > maxSector:
                        
                            ##  Done reading now

                            break

                        sectors.append(sector)
                        sectorEigs.append(eigenvalue)

                        ##print(currSector)

                fin.close()
                
                if currSector == maxSector:
                
                    ##  Sort the eigenvalues for each sector in descending order

                    sectorEigs.sort()
                    
                    ##  Append to the full list of eigenvalues
                    
                    eigenvalues += sectorEigs
                  
                ##	check that we have sufficient data for the fitting
                
                if currSector<maxSector:
                
                    print ('\n\tERROR: only found '+str(max(sectors))+' sectors to fit to, but requested '+str(maxSector))
                
                    quit()
                
                self.allEigenvalues.append(eigenvalues)
                self.allSectors.append(sectors)
                
            else:
            
                print ('\n\tERROR: CANNOT FIND FILE '+programOptions.BuildRsesFileName(0,nbr_a))
                
                quit()
        
        print('\t...DONE!\n')
        
        ##  Optionally normalize to zero 
        
        if programOptions.rses_norm_zero:
        
            self.NormalizeToZero()
        
        ##  Optionally normalize from file values
        
        if programOptions.rses_norm_file:
            
            ##  First normalize to zero
            self.NormalizeToZero()
            
            ##  Then get new norms from file
            self.NormalizeFromFile(programOptions)

    ############################################################################
    ##	Normalize the RSES data so that the lowest lying entanglement eigenvalue
    ##  in the lowest sector is set to zero
    ##
    def NormalizeToZero(self):
        
        ##  First obtain the lowest lying eigenvalue in sector 0
        ##  and for each na value
        
        nbrIndex=0
        
        for sectors in self.allSectors:

            minEigenvalue=100000
            i=0

            for sector in sectors:
        
                if sector == 0:
                    
                    if self.allEigenvalues[nbrIndex][i]<minEigenvalue:
                    
                        minEigenvalue = self.allEigenvalues[nbrIndex][i]
                    
                i+=1
            
            ##  Subtract that value from the whole spectrum
        
            i=0
            for eig in self.allEigenvalues[nbrIndex]:
            
                self.allEigenvalues[nbrIndex][i] = eig - minEigenvalue;    
                i+=1
            nbrIndex+=1
            
    ############################################################################
    ##	Normalize the RSES data so that the lowest lying entanglement eigenvalue
    ##  in the lowest sector is set to a value read from an extrnal file
    ##
    def NormalizeFromFile(self,programOptions):
        
        ##  First read in norm values from a file
        
        norms = []

        for nbr_a in programOptions.possibleNa:

            if os.path.isfile(programOptions.BuildRsesNormFileName(programOptions.rses_norm_file_label,nbr_a)):

                print ('\n\tReading in RSES norm data from '+programOptions.BuildRsesNormFileName(programOptions.rses_norm_file_label,nbr_a)+'\n')

                fin=open(programOptions.BuildRsesNormFileName(programOptions.rses_norm_file_label,nbr_a))
                
                ##  Skip the comment line        
                fin.readline()
                
                ##  Get the norm value from the file
                norms.append(float(fin.readline()))
                    
                fin.close()
                
            else:
            
                print ('\n\tERROR: CANNOT FIND FILE '+programOptions.BuildRsesNormFileName(programOptions.rses_norm_file_label,nbr_a))
                quit()
           
        print('\t...DONE!\n')
        
        ##  Rescale the spectrum by the list of norms
        
        nbrIndex=0
        
        for sectors in self.allSectors:

            ##  Subtract that value from the whole spectrum
        
            i=0
            for eig in self.allEigenvalues[nbrIndex]:
            
                self.allEigenvalues[nbrIndex][i] = eig + norms[nbrIndex];    
                i+=1
            nbrIndex+=1
    
    ############################################################################
    ##  A function to write the minimum entanglement energy values of
    ##  sector zero to a new file
    ##
    def EntanglementEnergyMinimaToFile(self,programOptions):

        nbrIndex=0
    
        for nbr_a in programOptions.possibleNa:

            ##  Find the minimum eigenvalue in sector 0

            minEigenvalue=100000
            i=0

            for sector in self.allSectors[nbrIndex]:
        
                if sector == 0:
                    
                    if self.allEigenvalues[nbrIndex][i]<minEigenvalue:
                    
                        minEigenvalue = self.allEigenvalues[nbrIndex][i]
                i+=1

            print ('\n\tWriting RSES norm data to '+programOptions.BuildRsesNormFileName('rses',0,nbr_a)+'\n')

            fout=open(programOptions.BuildRsesNormFileName('rses',nbr_a),'w')
            
            fout.write('##  THIS FILE CONTAINS THE MINIMUM ENTANGLEMENT ENERGY VALUE IN SECTOR 0:\n')
            
            fout.write(str(minEigenvalue)+'\n')
            
            fout.close()
            
            nbrIndex+=1
                
    ############################################################################
    ##	Use the counting from the fitting model RSES to count out the exact 
    ##  number of lowest lying eigenvalues for each sector. This process
    ##	will attempt to remove the upper continuum of states if present.
    ##
    ##  Additionally, generate a list of counting for each sector
    ##  to be used to weight the fitting function appropriately
    ##
    def RemoveContinuum(self,allFitSectors):

        allNewEigs = []
        self.allFitFunctionWeights = []
        
        for k in range(0,len(self.allEigenvalues)):

            eigenvalues = self.allEigenvalues[k]
            sectors     = self.allSectors[k]
            fitSectors  = allFitSectors[k]

            counter = 0
            sectorCounter = 0
            currSector = 0 
            i = 0
            
            newEigs=[]
            counting=[]

            for eig in eigenvalues:

                if sectors[i] == fitSectors[counter]:
                
                    newEigs.append(eig)
                    counter+=1
                    
                    if sectors[i] == currSector:
                    
                        sectorCounter+=1
                        
                    else:
                    
                        counting.append(sectorCounter)
                        currSector-=1
                        sectorCounter=1

                    if counter>=len(fitSectors):

                        currSector-=1
                    
                        break
                    
                i+=1
                
            counting.append(sectorCounter)

            ##  Now determine fitting function weights 

            fitFunctionWeights = []
            
            i=0
            
            for sector in fitSectors:

                fitFunctionWeights.append(math.sqrt(counting[-sector]))
             
                i+=1
             
            ##print (fitFunctionWeights)
            ##quit()
            ##*math.sqrt(weightsList[-sector])
            
            ##quit()
            
            allNewEigs.append(newEigs)
            self.allFitFunctionWeights.append(fitFunctionWeights)     
        
        ##  Overwrite the existing class self.allEigenvalues list
        ##  With the truncated list
        
        self.allSectors = allFitSectors[:]
        self.allEigenvalues = allNewEigs[:]
    
    ############################################################################
    ##  Plot the RSES data stored within the class. Note that matplotlib
    ##  needs to be imported before this function can be called
    ##
    def Plot(self,axes,nbr_a,isInset,programOptions):
        
        numberOffsetY    = 0.5
        markerSize       = 100
        numberOffsetX    = 0
        annotateFontSize = 16
        spacingFactor    = 1
        
        if isInset:
            numberOffsetY = 0.01
            markerSize    = 100
        
        if programOptions.plot_rses:
            markerSize    = 140
            numberOffsetX = 0.5
        
        if programOptions.nbr_ll == 2:
            
            numberOffsetY = 1.75       ##  1.5 for nu=2
        
        if programOptions.plot_degeneracy:
        
            numberOffsetX    = 0.18
            numberOffsetY    = 0.1
            annotateFontSize = 8
            spacingFactor    = 1.3
            
            if isInset:
                numberOffsetY = 0.01
                markerSize    = 100
                spacingFactor = 1.01

        if programOptions.rses_method == 'ewf':
            
            ##fittedDataLabel=r"EWF METHOD"
            fittedDataLabel=r"Numerical Method"
            
        elif programOptions.rses_method == 'exact':
        
            fittedDataLabel=r"EXACT"
        else:
            fittedDataLabel= r"RSES"
        
        ##  Determine which index corresponds to nbr_a
        
        nbrIndex = programOptions.GetPossibleNaIndex(nbr_a)
            
        axes.scatter(self.allSectors[nbrIndex],self.allEigenvalues[nbrIndex],s=markerSize,c='black',marker='_',lw=0.4,label=fittedDataLabel)
    
        ##  Plot degeneracy number label next to each set of levels

        if programOptions.plot_counting or programOptions.plot_degeneracy:
        
            if programOptions.plot_counting:
                branchGap = 1.8     ##  1.0 for nu=2, 1.8 for 2/3
            elif programOptions.plot_degeneracy:
                branchGap = 0.00001     ##  Set small to test for degeneracies
            
            minCount   = 0
            
            if programOptions.nbr_ll == 2:
                maxSector = 5
            else:
                maxSector  = 9
            currLza    = 0
            count      = 0
            sectorCount = 1
            countList  = []
            sectorList = []
            energyList = []
            labelEnergyList = []
            innerEnergyList = []

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
                    
                    #   Determine the energy range for this sector
                    
                    print(sectorCount)
                    
                    print(energyList[-sectorCount:])
                    
                    minRange = np.min(energyList[-sectorCount:])
                    maxRange = np.max(energyList[-sectorCount:])
                    
                    #   Use the range to set equally spaced label 
                    #   positions
                    
                    midRange = 0.5*(minRange+maxRange)
                    
                    ##  rescale
                    
                    if programOptions.nbr_ll == 2:
                        maxRange *= spacingFactor+sectorCount/200.0
                    else:    
                        maxRange *= spacingFactor+sectorCount/300.0
                    minRange = 2*midRange - maxRange
                    
                    for r in np.linspace(minRange,maxRange,sectorCount):
                    
                       labelEnergyList.append(r)
                        
                    innerEnergyList = []
                    currLza -= 1
                    count = 0
                    sectorCount = 1
                    
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

                        innerEnergyList = []
                        count = 0
                        sectorCount += 1
                
                count += 1
                innerEnergyList.append(e)
            
            print("Count of plotted eigenvalues")
            print(countList)
            
            print(labelEnergyList)
            '''
            for i in range(0,len(countList)):
                    
                axes.annotate(
                str(countList[i]),
                xy=(sectorList[i]+numberOffsetX,energyList[i]),
                xytext=(sectorList[i]+3*numberOffsetX,
                labelEnergyList[i]-numberOffsetY),
                arrowprops=dict(arrowstyle='-',color='black',lw=0.001,shrinkA=0.1,connectionstyle='arc,rad=0,angleA=-20'),
                size=annotateFontSize)
            '''
            
    #############################################################################
    ##  Calculate the sum of exp(-eigenvalue) in order to check the overall
    ##  rses normalization
    ##
    def CalculateNormalizations(self,normSector,programOptions):

        print("\n\tCHECKING SUM OF EXP(-E) AND LOWEST ENTANGLEMENT ENERGY STATE\n")

        self.GetSpectrumFromFile(programOptions,normSector)

        currTotal = 0.0
        nbrIndex = 0
        norms = []

        for nbr_a in programOptions.possibleNa:

            counter = 0
            sectorCounter = 0
            currSector = 0
            
            sectorTotals = []
           
            i = 0

            for eig in self.allEigenvalues[nbrIndex]:

                if i==0:
                    
                    norm = self.allEigenvalues[nbrIndex][i]
                
                    norms.append(norm)

                if self.allSectors[nbrIndex][i] == currSector:
                
                    currTotal += math.exp(-eig)

                else:
                
                    sectorTotals.append(currTotal)
                    
                    currTotal += math.exp(-eig)

                    currSector -= 1
                    
                i += 1
            
            nbrIndex += 1

            print("\n\tSector totals of exp(-E):",sectorTotals)
       
        print("\n\tOverall total of exp(-E):",currTotal)
                
        print ("\n\tLowest entanglement energy states:",norms)
        
        print ("\n")
        
        return norms

    #############################################################################
    ##  Count out the number of RSES states below a given entanglement energy
    ##
    def GetCounting(self,maxEnergy):
    
        counting = []
        nbrIndex=0
        
        print(self.allSectors)
        
        for sectors in self.allSectors:

            sectorCounting = 0
            i=0

            for sector in sectors:
                if self.allEigenvalues[nbrIndex][i]<maxEnergy:            
                    sectorCounting+=1           
                i+=1
            counting.append(sectorCounting)
            
        return counting      
    
################################################################################
##  Define a class to import and contain data from the RSES calculated with
##  the entanglement wave function (ewf) method
##
class EwfRses(RsesData):
    pass

################################################################################
##  Define a class to import and contain data from the RSES calculated with
##  the pseudo reduced density matrix (rdm) method
##
class RdmRses(RsesData):

    ############################################################################
    ##  A function to read in the RSES data from a file. Pass a ProgramOptions
    ##  class as the programOptions argument
    ##
    def GetSpectrumFromFile(self,programOptions,maxSector):
    
        self.allEigenvalues = []
        self.allSectors = []

        for nbr_a in programOptions.possibleNa:

            #	find first sector file

            condition=1
            sector=0

            while condition==1:
                if os.path.exists(programOptions.BuildRsesFileName(sector,nbr_a)):
                    condition=0
                    break
                else:
                    condition=1
                    sector+=1
                    
                    ##print (rsesName+str(sector)+'.dat')
                    
                    if sector>500:
                    
                        print ('CANNOT FIND FOLDER '+programOptions.rdm_rses_path)
                    
                        quit()
                    
            firstSector=sector+programOptions.skipSector 	

            eigenvalues=[]
            sectors=[]

            while condition!=1:

                fin=open(rsesName+str(sector)+'.dat')

                for line in fin:

                    eigenvalue = -float(line)

                    eigenvalues.append(eigenvalue)
                    sectors.append(-sector+firstSector)

                fin.close()

                print ('read sector '+str(sector))

                sector+=1

                #   check for next file

                if os.path.isfile(rsesName+str(sector)+'.dat'):
                    condition = 0
                else:
                    condition = 1  
                    
                if sector-firstSector>maxSector:	
                    condition = 1

            self.allEigenvalues.append(eigenvalues)
            self.allSectors.append(allSectors)
        
        print('\t...DONE!\n')
        
################################################################################
##  Define a class to import and contain data from the RSES calculated with
##  a numerical exact method
##
class ExactRses(RsesData):
    pass
