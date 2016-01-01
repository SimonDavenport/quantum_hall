#! /usr/bin/env python2.7
################################################################################
##
##						  Author: Simon C. Davenport
##
##							Last Modified: 04/11/2014
##
##
##		This file contains some basic curve fitting functions to fit data
##      to exponential, or polynomial functions and determine an 
##      extrapolation to some specified point.
##
##      Import the library with:
##      import sys
##      sys.path.append(r'/scratch/scd51/Dropbox/physics_programs/python')
##
##      from utilities import curve_fitting
##
##      Then refer to the library variables with curve_fitting.NAME
##
################################################################################

################################################################################
################################################################################
##  Import modules

from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import numpy as np
import math

################################################################################
##    Extrapolation function in a+b/x
##
##    Use as fitData,errors,extrap = InverseLinearFit(xData,yData,extrapPoint)
##
def InverseLinearFit(xData,yData,extrapPoint):    
    
    xData=np.transpose(xData)    
    yData=np.transpose(yData)

    fitfunc = lambda x,a,b: a+b/x
    
    params,covar=curve_fit(fitfunc,xData,yData,None)
    
    output='\n\tFitted data:\n'
    
    for i in range(0,len(xData)):
    
        output+='\n\t'+str(xData[i])+'\t'+str(yData[i])
    
    output+='\n\n\tInverse linear fit parameter values:\n'
    output+='\n\tIntercept\t'+str(params[0])+' +/- '+str(math.sqrt(covar[0][0]))
    output+='\n\tGradient\t'+str(params[1])+' +/- '+str(math.sqrt(covar[1][1]))

    print (output)

    fitdata=fitfunc(xData,params[0],params[1])

    return fitdata,math.sqrt(covar[0][0]),fitfunc(extrapPoint,params[0],params[1])

################################################################################
##    Linear extrapolation functionin a+bx
##
##    Use as fitData,errors,extrap = LinearFit(xData,yData,extrapPoint)
##
def LinearFit(xData,yData,extrapPoint):    
    
    xData=np.transpose(xData)
    yData=np.transpose(yData)

    fitfunc = lambda x,a,b: a+b*x
    
    params,covar=curve_fit(fitfunc,xData,yData,None)
    
    output='\n\tFitted data:\n'
    
    for i in range(0,len(xData)):
    
        output+='\n\t'+str(xData[i])+'\t'+str(yData[i])
    
    output+='\n\n\tLinear fit parameter values:\n'
    output+='\n\tIntercept\t'+str(params[0])+' +/- '+str(math.sqrt(covar[0][0]))
    output+='\n\tGradient\t'+str(params[1])+' +/- '+str(math.sqrt(covar[1][1]))

    print (output)

    fitdata=fitfunc(xData,params[0],params[1])

    return fitdata,math.sqrt(covar[0][0]),fitfunc(extrapPoint,params[0],params[1])
    
################################################################################
##    Extrapolation function in a+b/x+c/x^2
##
##    Use as fitData,errors,extrap = InverseQuadraticFit(xData,yData,extrapPoint)
##
def InverseQuadraticFit(xData,yData,extrapPoint):    
    
    xData=np.transpose(xData)    
    yData=np.transpose(yData)

    fitfunc = lambda x,a,b,c: a+b/x+c/x**(2)
    
    params,covar=curve_fit(fitfunc,xData,yData,None)
    
    fitdata=fitfunc(xData,params[0],params[1],params[2])
    
    output='\n\tFitted data:\n'
    
    for i in range(0,len(xData)):
    
        output+='\n\t'+str(xData[i])+'\t'+str(yData[i])+'\t'+str(fitdata[i])
    
    output+='\n\t'+str(extrapPoint)+'\t\t\t'+str(fitfunc(extrapPoint,params[0],params[1],params[2]))
    
    output+='\n\n\tInverse quadratic fit parameter values:\n'
    output+='\n\tIntercept\t'+str(params[0])+' +/- '+str(math.sqrt(covar[0][0]))
    output+='\n\tGradient\t'+str(params[1])+' +/- '+str(math.sqrt(covar[1][1]))
    output+='\n\tCurvature\t'+str(params[2])+' +/- '+str(math.sqrt(covar[2][2]))

    print (output)

    return fitdata,math.sqrt(covar[0][0]),fitfunc(extrapPoint,params[0],params[1],params[2])

################################################################################
##    Exponetial extrapolation function a exp(-x) + b
##
##    Use as extrapXdata,fitData,result,error = ExponentialFit(xData,yData)
##
def ExponentialFit(xData,yData,extrapPoint):
    
    xData=np.transpose(xData)    
    yData=np.transpose(yData)

    fitfunc = lambda x,a,b: (a)*np.exp(-x)+b
    
    params,covar=curve_fit(fitfunc,xData,yData,None)
    
    fitdata=fitfunc(xData,params[0],params[1])
    
    output='\n\tFitted data:\n'
    
    for i in range(0,len(xData)):
    
        output+='\n\t'+str(xData[i])+'\t'+str(yData[i])
    
    output+='\n\n\tExp fit parameter values:\n'
    output+='\n\ta Exp(-x)\t'+str(params[0])+' +/- '+str(math.sqrt(covar[0][0]))
    output+='\n\tConstant\t'+str(params[1])+' +/- '+str(math.sqrt(covar[1][1]))

    print (output)

    return fitdata,math.sqrt(covar[1][1]),fitfunc(extrapPoint,params[0],params[1])
    
