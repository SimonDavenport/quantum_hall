#! /usr/bin/env python2.7
################################################################################
##
##						  Author: Simon C. Davenport
##
##							Last Modified: 04/11/2014
##
##
##		This file defines a python module containing some basic formatting 
##      and colouring options to be used with pthon's Matplotlib library.
##
##      Import the library with:
##      import sys
##      sys.path.append(r'/scratch/scd51/Dropbox/physics_programs/python')
##
##      from utilities import plot_setup
##
##      Then refer to the library variables with plot_setup.NAME
##
################################################################################

################################################################################
################################################################################
##  Import modules

##  Import numpu for numpy linspace
import numpy as np

## Import the os library to allow us to check if displaying the plot is possible
import os

## Import matplotlib library
import matplotlib

##  Check to see if the DISPLAY environment variable is set. If it is
##  then plt.show() will be possible
if "DISPLAY" not in os.environ.keys():
    
    ## Turn off X-forwarding if not available [plt.show() disabled]
    ## Must be before importing matplotlib.pyplot or pylab!
    
    matplotlib.use('Agg')
   
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MultipleLocator

################################################################################
################################################################################
##	Specify figure fonts

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 19}

##  Set latex fonts	
matplotlib.rc('text', usetex=True)	

##  Define plot file extension
file_ext='.pdf'

################################################################################
################################################################################
##  Colour-blindness friendly pallette

##  15 colours are defined here
colour_list = ['#bd2309', '#bbb12d', '#1480fa', '#14fa2f','#faf214',
               '#ea2ec4','#2edfea',  '#ea2e40', '#cdcdcd','#577a4d',
               '#2e46c0', '#f59422', '#219774', '#8086d9', '#000000']             

################################################################################
################################################################################
##  Specify a list of plot markers
     
##  15 markers are defined here
marker_list = ['.','^','s','o','>','<','*','d','h','1','2','3','4','p','x']

################################################################################
################################################################################
##  Define a discrete colour map for use in plot colour bars
##
def DiscreteColourMap(valueList):

    ##  Find the unique set of values in the list and map these to a 
    ##  discrete colour map

    colourList = []

    seen = {}
    unique = []
    
    for item in valueList:
        marker = item

        if marker in seen: continue
        seen[marker] = 1
        unique.append(item)
    
    if len(unique)>15:
        print("TOO MANY COLOURS REQUIRED!")
        quit()
        
    ##dictionary = dict(zip(unique,COLOUR_LIST[:len(unique)]))

    ##return [dictionary[k] for k in valueList]

    unique.sort()

    dictionary = dict(zip(unique,range(0,len(unique))))

    ##  Define a colour map

    colList = colour_list[:len(unique)]

    myColourMap = matplotlib.colors.ListedColormap(colList)

    # define the bins and normalize
    bounds = np.linspace(0,len(colList),len(colList)+1)
    myNorm = matplotlib.colors.BoundaryNorm(bounds, myColourMap.N)
    
    ##print(COLOUR_LIST[:len(unique)])
    ##print([dictionary[k] for k in valueList])
    
    return [dictionary[k] for k in valueList],myColourMap,bounds,myNorm

