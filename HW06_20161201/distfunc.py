'''
Created on Nov 27, 2016
@author: Mitch

Module containing probability density functions to be called in rejection sampling
'''
import math

def quartercircle(x):
    # Equation defining quarter of a circle
    y = math.sqrt(1-x**2)
    return y

def deriv_atan(x):
    # Equation defining the derivative of arctan
    y = 1/(1+x**2)
    return y