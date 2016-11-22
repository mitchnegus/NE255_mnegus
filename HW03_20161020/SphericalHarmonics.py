'''
Created on Oct 15, 2016

@author: Mitch
'''

import math
import cmath
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import size


def func_derivs(n,x,l):
# Function returning the nth derivative (up to n=4) of the function (x^2-1)^l evaluated at x
    if n == 0:
        Y = (x**2-1)**l
    elif n == 1:
        Y = 2*(l)*(x)*(x**2-1)**(l-1)
    elif n == 2:
        Y = 2*(l)*(2*(l-1)*(x**2)*(x**2-1)**(l-2)+(x**2-1)**(l-1))
    elif n == 3:
        Y = 2*(l)*(2*(l-1)*(2*(l-2)*(x**3)*(x**2-1)**(l-3)+2*(x)*(x**2-1)**(l-2))+2*(l-1)*x*(x**2-1)**(l-2))
    elif n == 4:
        Y = 2*(l)*(2*(l-1)*(2*(l-2)*(2*(l-3)*(x**4)*(x**2-1)**(l-4)+3*(x**2)*(x**2-1)**(l-3))+2*(2*(l-2)*(x**2)*(x**2-1)**(l-3)+(x**2-1)**(l-2)))   +2*(l-1)*(2*(l-2)*(x**2)*(x**2-1)**(l-3)+(x**2-1)**(l-2)))
                
    return Y

def ALP_mPos(l,m,x):
# A function returning the associated Legendre Polynomials for positive m values
    dydx = func_derivs(l+m,x,l)
    P = ((-1)**m)/(2**l*math.factorial(l))*(1-x**2)**(m/2)*dydx
    return P

def assocLegPoly(l,m,x):
# Function returning associated Legendre Polynomials (differentiating between (+) and (-) m values)
    if m >= 0:
        P = ALP_mPos(l,m,x)
    else:
        m = abs(m)
        P = (-1)**m*(math.factorial(l-m)/math.factorial(l+m))*ALP_mPos(l,m,x)
    return P

def SphHarmonics(l,m,x,phi):
# Function returning SH using associated Legndre Polynomials
    P = assocLegPoly(l,m,x)
    Y = (-1)**m*math.sqrt((2*l+1)/(4*math.pi)*math.factorial(l-m)/math.factorial(l+m))*P*cmath.exp(1j*m*phi)
    return Y

def mod_squared(Z):
# Function calculating the square modulus of complex Z
    modsquaredZ = Z.real**2 + Z.imag**2
    return modsquaredZ



fig = plt.figure(figsize=(18,12))
fig.suptitle('Spherical Harmonics',fontsize=40)
l = [0,1,2]

# Discretization of theta & phi
inc = 101
Theta = np.linspace(0.0001,math.pi-0.0001,inc)  # Offset included to avoid derivative where 0 is raised to negative power
Phi = np.linspace(0,2*math.pi,inc)

for l_i in l:
    for m in range(l_i+1):
        # Calculation of Y_lm(theta,phi)
        #     ***Note: mod^2 is independent of phi, so it is not looped over (value of zero used throughout)
        modsquaredYlm = np.zeros(inc)
        Ylmtable = []
        i=0
        for theta in Theta:
            Ylm = SphHarmonics(l_i, m, math.cos(theta), 0)
            Ylmtable.append(Ylm)
            modsquaredYlm[i] = mod_squared(Ylm)
            i += 1
        
        # Convert Theta,Phi,R=Y_lm to X,Y,Z coordinates 
        X = np.outer(modsquaredYlm*np.sin(Theta),np.cos(Phi))
        Y = np.outer(modsquaredYlm*np.sin(Theta),np.sin(Phi))
        Z = np.outer(modsquaredYlm*np.cos(Theta),np.ones(size(Phi)))

        subplot = 330+(1+3*l_i+m)
        ax = fig.add_subplot(subplot, projection='3d')
        ax.plot_surface(X,Y,Z,rstride=2,cstride=2,cmap=cm.RdYlGn,linewidth=0)
        
        ax.text2D(-0.09,0.06,r'$|Y^{%i}_{%i}|^{2}$' % (m,l_i),{'fontsize':28})
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_xticks([],[])
        ax.set_yticks([],[])
        ax.set_zticks([],[])
        

fig.savefig('<file_location>')
plt.show()


