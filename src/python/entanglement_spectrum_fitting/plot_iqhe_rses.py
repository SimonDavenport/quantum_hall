#! /usr/bin/env python2.7
################################################################################
##                                                                            ##
##						  Author: Simon C. Davenport						  ##
##                                                                            ##
##							Last Modified: 01/04/2015						  ##
##                                                                            ##
##                                                                            ##
##		Calculate and plot the RSES for integer quanutm Hall states           ##
##      and its Tayor expansion                                               ##
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

##  Import additional standard python libraries
import math
import mpmath
import numpy as np

import os
import argparse

from scipy.integrate import quad
from scipy import interpolate

def FactorialRatio(input1,input2):

    if input1==input2:
        return 1.0

    if input1<0:
        return 0

    cmp=1

    if input1==0:
	    output=input2
	    while input2-cmp>1:
	        output*=(input2-cmp)
	        cmp+=1
	    output=1.0/output;

    if input1>input2:
	    output=input1;
	    while input1-cmp>input2:
	        output*=(input1-cmp)
	        cmp+=1
    elif input2>input1:
	    output=input2;
	    while input2-cmp>input1:
	        output*=(input2-cmp)
	        cmp+=1
	    output=1.0/output;
	
    return output

def MonopoleHarmonic(q,l,m,x):

    total = 0;
    maximum = int(q+l-m);
    
    for s in range(0,maximum+1):

        ##  Calculate -1^(q+l-m-s)
    
        if int(q+l-m-s) & 1:
            sign = -1
        else:
            sign = 1

        total +=  sign*mpmath.binomial(l,s)*mpmath.binomial((int)(2*q+l),(int)(q+l-m-s))*((1-x)**(l-s))*((1+x)**s)

    total *= ((1.0-x)**((q-m)/2.0)) * ((1.0+x)**((q+m)/2.0))
    
    return total

def MonopoleHarmonicNorm(q,l,m):

    norm = (2.0**(-q-l))*math.sqrt((2*(q+l)+1)/(4*math.pi));
    
    if (q+l-m)>(q+l+m):
        norm*=mpmath.sqrt(FactorialRatio((q+l-m),(2.0*q+l))*FactorialRatio((q+l+m),l))
    else:
        norm*=mpmath.sqrt(FactorialRatio((q+l-m),l)*FactorialRatio((q+l+m),(2.0*q+l)))

    return norm

def Integrand(x,q,m,ll1,ll2):

    return MonopoleHarmonic(q,ll1,m,x)*MonopoleHarmonic(q,ll2,m,x)

################################################################################
################################################################################
######      MAIN        ########################################################

if __name__=='__main__':

    formatting = 1

    if formatting ==1:
        N  = 100             ##  Number of particles
        LL = 1               ##  Number of Landau levels
        Q = (N-LL**2)/(2*LL) ## Monopole strength
        nStates = 2*(Q+LL)-1
    elif formatting==2:
        N  = 60              ##  Number of particles
        LL = 2               ##  Number of Landau levels
        Q = (N-LL**2)/(2*LL) ## Monopole strength
        nStates = 2*(Q+LL)-3
    
    ##  Construct and diagonalize the overlap matrix M^m_{sigma,sigma'}

    epsilon = np.zeros(shape=(nStates,LL))
    '''
    fin=open("current_data.dat")
    
    for m in range(0,int(nStates)):
        for i in range(0,LL):
            epsilon[m,i] = float(fin.readline())
        
    fin.close()
    '''
    
    if formatting ==1:
        zAngularMomentum = -Q-LL+1
    elif formatting==2:
        zAngularMomentum = -Q-LL+2
    
    w = mpmath.matrix(LL,1)
    
    w[0] = 0
    #w[1] = 0
    
    for m in range(0,int(nStates)):
    
        print(str(math.ceil(float(100*m)/nStates))+'%\r')
        
        if LL==1 and w[0]>0.6:
            for i in range(0,LL):
            
                epsilon[m,i] = -epsilon[nStates-m-1,i]
                
                #print(epsilon[m,i]) 
                #print(mpmath.log(w/(1-w))) 
        else:
        
            M = mpmath.matrix(LL,LL)

            print("zAngularMomentum="+str(zAngularMomentum))

            for i in range(0,LL):
            
                iNorm = MonopoleHarmonicNorm(Q,i,zAngularMomentum)
                ##print(iNorm)
            
                for j in range(0,LL):
                
                    jNorm = MonopoleHarmonicNorm(Q,j,zAngularMomentum)
                    
                    #print(jNorm)
                    
                    ##args = (Q,zAngularMomentum,i,j)

                    ##print(jNorm)
                    ##print(integrate.quad(Integrand,0.0,1.0,args)[0])
                    
                    M[i,j] = 2.0*math.pi*iNorm*jNorm*mpmath.quad(lambda x: Integrand(x,Q,zAngularMomentum,i,j),[0.0,1.0])
                    
                    ##M[i,j] = 2.0*math.pi*iNorm*jNorm*quad(lambda x: Integrand(x,Q,zAngularMomentum,i,j),0.0,1.0)
                    
                    ##print(2.0*math.pi*iNorm*jNorm*integrate.quad(Integrand,0.0,1.0,args)[0])
            print(M)
            
            if LL==1:
                w[0] = M[0,0]
            elif LL==2:
                desc = mpmath.sqrt(M[0,0]**2 + 4*M[0,1]*M[1,0]-2*M[0,0]*M[1,1]+M[1,1]**2)
                
                w[0] = 0.5*(M[0,0]+M[1,1]-desc)
                w[1] = 0.5*(M[0,0]+M[1,1]+desc)
            
            ##w = M[0,0]
            
            ##print(w)

            for i in range(0,LL):

                if mpmath.isinf(mpmath.log(w[i]/(1-w[i]))):
                    epsilon[m,i] = mpmath.nan
                else:
                    epsilon[m,i] = mpmath.log(w[i]/(1-w[i]))
                
                ##print(epsilon[m,i]) 
                
                ##mpmath.log(w[i]/(1-w[i]))
                    
            zAngularMomentum += 1.0
        
    fout=open("current_data.dat",'w')
    
    for i in range(0,LL):   
        for val in epsilon[:,i]:
            print(val)
            fout.write(str(val)+"\n") 
    fout.close()
    #'''
    font_size = 36
    
    plotX = np.arange(-1,1+2.0/(nStates-1),2.0/(nStates-1))
    
    if formatting==2:
        ##  Fix up some numerical issues with the tail of the data on the
        ##  upper branch
        epsilon[-1,1] = -epsilon[0,0]
        epsilon[-2,1] = -epsilon[1,0]
        ##epsilon[-3,1] = -epsilon[2,0]

    print(epsilon[:,0])
    if LL>1:
        print(epsilon[:,1])
    print(plotX)

    ##print(len(epsilon[:,0]))
    ##print(len(plotX))
    
    main_axes = plot_utilities.plt.axes([0.09,0.14,0.75,0.75])
    
    if formatting==1:
        main_axes.plot(plotX,epsilon[:,0],c=plot_utilities.colour_list[1],label=r'$\epsilon_{Q,m,0}$',lw=3)
    elif formatting==2:
        main_axes.plot(plotX,epsilon[:,0],c=plot_utilities.colour_list[1],label=r'$\epsilon_{Q,m,\sigma}$',lw=3)
    
        main_axes.plot(plotX,epsilon[:,1],c=plot_utilities.colour_list[1],lw=3)

    ##  Plot various levels of approximation
    
    plotY = []
    
    for x in plotX:
        
        if formatting==1:
            plotY.append(15.7848*x)
        elif formatting==2:
            plotY.append(-2.22159+11.21*x)

    if formatting==1:
        main_axes.plot(plotX,plotY,'-.',c=plot_utilities.colour_list[2],label=r'$b m$',lw=3)
    elif formatting==2:
        main_axes.plot(plotX,plotY,'-.',c=plot_utilities.colour_list[2],label=r'$b m +[\sigma-1/2] e$',lw=3)
    
    if formatting==2:
        plotY = []
        
        for x in plotX:
                plotY.append(2.22159+11.21*x)

        main_axes.plot(plotX,plotY,'-.',c=plot_utilities.colour_list[2],lw=3)

    plotY = []
    
    for x in plotX:
        if formatting==1:
            plotY.append(15.7848*x+71.6213*(x**3))
        elif formatting==2:
            plotY.append(-2.22159+11.21*x-7.34286*(x**2))

    if formatting==1:
        main_axes.plot(plotX,plotY,'--',c=plot_utilities.colour_list[3],label=r'$b m + d m^3$',lw=3)
    elif formatting==2:
        main_axes.plot(plotX,plotY,'--',c=plot_utilities.colour_list[3],label=r'$b m  + [\sigma-1/2] (e + g m^2)$',lw=3)
    
    if formatting==2:
        plotY = []
        
        for x in plotX:
            plotY.append(2.22159+11.21*x+7.34286*(x**2))

        main_axes.plot(plotX,plotY,'--',c=plot_utilities.colour_list[3],lw=3)
    
    plotY = []
    
    for x in plotX:
        
        if formatting==1:
            plotY.append(15.7848*x+71.6213*(x**3)+10.6528*(x**5))
        elif formatting==2:
            plotY.append(-2.22159+11.21 *x-7.34286*(x**2)+7.11101*(x**3))
    
    if formatting==1:
        main_axes.plot(plotX,plotY,':',c=plot_utilities.colour_list[5],label=r'$b m + d m^3 + f m^5$',lw=3)
    elif formatting==2:
        main_axes.plot(plotX,plotY,':',c=plot_utilities.colour_list[5],label=r'$b m + d m^3 + [\sigma-1/2] (e + g m^2) $',lw=3)
    
    if formatting==2: 
        plotY = []
        
        for x in plotX:
            plotY.append(2.22159+11.21 *x+7.34286*(x**2)+7.11101*(x**3))

        main_axes.plot(plotX,plotY,':',c=plot_utilities.colour_list[5],lw=3)
    
    if formatting==1: 
        main_axes.legend(loc='best',numpoints=1,ncol=1,fontsize=font_size-12)
    elif formatting==2:
        main_axes.legend(loc='best',numpoints=1,ncol=1,fontsize=font_size-16)
        
    main_axes.set_xlim([-1,1])
    main_axes.set_xlabel(r'$m/Q$',fontsize=font_size)
    main_axes.set_xticks(np.arange(-1,1.2,0.5))
    main_axes.set_xticklabels([r'$-\textbf{1.0}$',r'$-\textbf{0.5}$',r'$\textbf{0.0}$',r'$\textbf{0.5}$',r'$\textbf{1.0}$'],fontsize=font_size-10)

    if formatting==1:
        main_axes.set_ylim([-100,115])
        ##main_axes.set_ylabel(r'$\epsilon_{m,0}$',fontsize=font_size)
        main_axes.set_yticks(np.arange(-100,125,25))
        main_axes.set_yticklabels([r'$-\textbf{100}$',r'$-\textbf{75}$',r'$-\textbf{50}$',r'$-\textbf{25}$',r'$\textbf{0}$',r'$\textbf{25}$',r'$\textbf{50}$',r'$\textbf{75}$',r'$\textbf{100}$'],fontsize=font_size-10)

    elif formatting==2:
        main_axes.set_ylim([-30,50])
        main_axes.set_yticks(np.arange(-30,40,10))
        main_axes.set_yticklabels([r'$-\textbf{30}$',r'$-\textbf{20}$',r'$-\textbf{10}$',r'$\textbf{0}$',r'$\textbf{10}$',r'$\textbf{20}$',r'$\textbf{30}$'],fontsize=font_size-10)
        
        main_axes.annotate(r"$\sigma=1$",xy=(-0.3,6),fontsize=font_size-10)
        main_axes.annotate(r"$\sigma=0$",xy=(0.2,-6),fontsize=font_size-10)

    main_axes.axhline(y=0,linestyle='--',color='black')
    main_axes.axvline(x=0,linestyle='--',color='black')

    main_axes.tick_params(direction='in',pad=10,axis='y')
    main_axes.tick_params(direction='in',pad=10,axis='x')

    if formatting==1:
        plot_utilities.plt.savefig("nu_1_entanglement_energy_function.pdf",bbox_inches='tight')
    elif formatting==2:
        plot_utilities.plt.savefig("nu_2_entanglement_energy_function.pdf",bbox_inches='tight')

    '''

    mpmath.dps = 100

    nbrList = [100,200,300,400,500]

    for N in nbrList:
    
        #N  = 30     ##  Number of particles
        LL = 1      ##  Number of Landau levels
        Q = (N-LL**2)/(2*LL) ## Monopole strength
        nStates = 1+2*Q+LL-1
        
        ##  Construct and diagonalize the overlap matrix M^m_{sigma,sigma'}

        epsilon = np.zeros(shape=(nStates,LL))

        zAngularMomentum = -Q-LL+1
        
        w = 0
        
        for m in range(0,nStates):
        
            print(str(math.ceil(float(100*m)/nStates))+'%\r')
            
            if w>0.7:
                for i in range(0,LL):
                
                    epsilon[m,i] = -epsilon[nStates-m-1,i]
                    
                    print(epsilon[m,i]) 
                    print(mpmath.log(w/(1-w))) 
            else:
            
                M = mpmath.matrix(LL,LL)
            
                for i in range(0,LL):
                
                    iNorm = MonopoleHarmonicNorm(Q,i,zAngularMomentum)
                    ##print(iNorm)
                
                    for j in range(0,LL):
                    
                        jNorm = MonopoleHarmonicNorm(Q,j,zAngularMomentum)
                        
                        #print(jNorm)
                        
                        ##args = (Q,zAngularMomentum,i,j)

                        ##print(jNorm)
                        ##print(integrate.quad(Integrand,0.0,1.0,args)[0])
                        
                        M[i,j] = 2.0*math.pi*iNorm*jNorm*mpmath.quad(lambda x: Integrand(x,Q,zAngularMomentum,i,j),[0.0,1.0])
                        
                        ##print(2.0*math.pi*iNorm*jNorm*integrate.quad(Integrand,0.0,1.0,args)[0])

                ##print(M)
                
                ##w,v = mpmath.eig(M)
                
                w = M[0,0]
                
                ##print(w)

                for i in range(0,LL):
                
                    epsilon[m,i] = mpmath.log(w/(1-w))
                    
                    ##print(epsilon[m,i]) 
                    
                    ##mpmath.log(w[i]/(1-w[i]))
                                
                    
            zAngularMomentum += 1.0
        
        ##print(epsilon[:,0])
        ##print(np.arange(-0.5,0.5+1.0/(nStates-1),1.0/(nStates-1)))
        
        plot_utilities.plt.plot(np.arange(-0.5,0.5+1.0/(nStates-1),1.0/(nStates-1)),epsilon[:,0])
    
    plot_utilities.plt.show()
    
    '''
    
