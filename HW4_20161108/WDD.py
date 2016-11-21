import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colorbar as cb
import matplotlib as mpl
import math


def WDD_calc1D(q_e,Sigma_t,Sigma_s,alpha,mu,h,incoming,scattering):
    # Function to calculate 1D flux using weighted diamond difference, given incoming flux
    #    Cell-center flux        Psi_x = (s_x + (2/(1+-alpha_x))(|mu|/h)Psi_(x-+1/2))/(sigma_x + (2/(1+-alpha_x))(|mu|/h))
    #    Outgoing flux           Psi_(x+-1/2) = (2/(1+-alpha_x))Psi_x - ((1-+alpha_x)/(1+-alpha_x))Psi_(x-+1/2)
    
    alpha = alpha*mu/abs(mu)            # sign before alpha depends on mu>0 or mu<0
    Psi_i = incoming
    A = 2/(1+alpha)
    
    B = A*(abs(mu)/h)
    Psi_c = (q_e + B*Psi_i + scattering)/(Sigma_t + B)                   # Cell-center flux calculation
    
    Psi_o =  A*Psi_c - (1 - alpha)/(1 + alpha)*Psi_i    # Outgoing flux calculation
    
    center = Psi_c
    outgoing = Psi_o
    
    return center, outgoing


def Scatter(Sigma_s,ccflux,n):
    # Function to calculate scattering into angle a of the nth cell from all angles
    A = len(ccflux[0,:])
    scattersum = 0
    for ccf_a in ccflux.T:
        scattersum += (1/A)*Sigma_s*ccf_a[n]    
    return scattersum


def Sweep(q_e,Sigma_t,Sigma_s,alpha,mu,a,h,ncells,ccflux,influx):
    # Function to sweep over mesh for a given angle mu
    if mu>0:        # sweep left to right (count up cellorder)
        cellorder = range(ncells)
    elif mu<0:      # sweep right to left (count down cellorder)
        cellorder = reversed(range(ncells))
    for n in cellorder:
        scattering = Scatter(Sigma_s,ccflux,n)
        centerflux, outflux = WDD_calc1D(q_e,Sigma_t,Sigma_s,alpha,mu,h,influx,scattering)
        ccflux[n,a] =  centerflux
        influx = outflux
        
    return ccflux,outflux
 

if __name__ == "__main__":
    # Define parameters
    q_e = 1.0
    Sigma_t = 1
    Sigma_s = 0.9
    Alpha = [0]
    mu = [0.2,0.7]
    H = [0.08,0.1,0.125,0.2,0.4]
    lowbound = 0.0
    upbound = 2.0
    startflux = 2.0    
    
    #An isotropic source means:
    q_e = q_e/(2*len(mu))
    
    for alpha in Alpha:
        print('\nalpha = ',alpha)
        fig = plt.figure(figsize=(12,10))
        subplotnum = 100*len(H)+10
     
        for h in H:
            print('\nMesh Spacing: ',h)
            
            # Generate an array to store calculated flux values (1 column each for mu_a and -mu_a
            ncells = int((upbound-lowbound)/h)
            ccflux = np.zeros((ncells,2*len(mu)))
            
            reflected = [0 for a in range(len(mu))]
            
            # Iterate until the scalar flux converges
            scalarflux = np.zeros(ncells)
            scalarfluxnorm2 = -1
            f = 0
            while (scalarfluxnorm2 > 0.000001 or scalarfluxnorm2 == -1 ) and f<1000:
                # Sweep (Left->Right), 1 sweep per angle
                for a in range(len(mu)):
                    ccflux,outflux = Sweep(q_e,Sigma_t,Sigma_s,alpha,mu[a],a,h,ncells,ccflux,startflux)
                    reflected[a] = outflux  # save the reflected flux to be used as the inputs in the other direction
                # Sweep (Right->Left), 1 sweep per angle
                for a in range(len(mu)):
                    ccflux,outflux = Sweep(q_e,Sigma_t,Sigma_s,alpha,-mu[a],(a+len(mu)),h,ncells,ccflux,reflected[a])
                
                
                scalarflux_prev = scalarflux
                scalarflux = np.sum(ccflux,axis=1)
                
                scalarfluxnorm2 = math.sqrt(sum((scalarflux - scalarflux_prev)**2))
                
                
                f+=1
            
            print('Scalar Flux:\n',scalarflux)
            print('Norm2: ',scalarfluxnorm2)
                
            # Plotting parameters
            subplotnum += 1
            ax = fig.add_subplot(subplotnum)
            ax.text(0.5,0.8,'Center-cell Scalar Flux, h=%s' %(str(h)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            ax.text(0.9,0.1,r'$\alpha$=%s' %(str(alpha)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            ax.set_ylabel('$\psi_{i}$',rotation=0)
            plt.plot(scalarflux)
            plt.plot([0 for i in range(ncells)])
            
        fig.savefig('C:/Users/Mitch/Documents/Cal/1 - 2016 Fall/NUCENG 255 - Numerical Simulation in Radiation Transport/Homework/HW_Drafts/HW_Figures/HW04_Prob3d/WDD_alpha%s.jpg' %(str(alpha)))
    