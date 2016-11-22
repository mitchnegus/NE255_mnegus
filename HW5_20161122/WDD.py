import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colorbar as cb
import matplotlib as mpl
import math


def WDD_calc1D(q_e,Sigma_t,Sigma_s,alpha,group,mu,h,incoming,scattering):
    # Function to calculate 1D flux using weighted diamond difference, given incoming flux
    #    Cell-center flux        Psi_x = (s_x + (2/(1+-alpha_x))(|mu|/h)Psi_(x-+1/2))/(sigma_x + (2/(1+-alpha_x))(|mu|/h))
    #    Outgoing flux           Psi_(x+-1/2) = (2/(1+-alpha_x))Psi_x - ((1-+alpha_x)/(1+-alpha_x))Psi_(x-+1/2)
    
    alpha = alpha*mu/abs(mu)            # sign before alpha depends on mu>0 or mu<0
    Psi_i = incoming
    A = 2/(1+alpha)
    
    B = A*(abs(mu)/h)
    Psi_c = (q_e[group] + B*Psi_i + scattering)/(Sigma_t[group] + B)                   # Cell-center flux calculation
    
    Psi_o =  A*Psi_c - (1 - alpha)/(1 + alpha)*Psi_i    # Outgoing flux calculation
    
    center = Psi_c
    outgoing = Psi_o
    
    return center, outgoing


def Scatter(group,Sigma_s,mg_ccflux,n):
    # Function to calculate scattering into angle a of the nth cell from all angles
    A = len(mg_ccflux[group,0,:])
    ngroups = len(mg_ccflux)
    
    scattersum = 0
    for primegroup in range(ngroups):
        for ccf_a in mg_ccflux[primegroup,n,:]:
            scattersum += (2/A)*Sigma_s[group][primegroup]*ccf_a
    return scattersum


def Sweep(q_e,Sigma_t,Sigma_s,alpha,group,mu,a,h,ncells,mg_ccflux,influx):
    
    calc_ccflux = np.zeros_like(mg_ccflux[group,:,a]) # an array to store the newly calculated cell-center fluxes
    # Function to sweep over mesh for a given angle mu
    if mu>0:        # sweep left to right (count up cellorder)
        cellorder = range(ncells)
    elif mu<0:      # sweep right to left (count down cellorder)
        cellorder = reversed(range(ncells))
    for n in cellorder:
        scattering = Scatter(group,Sigma_s,mg_ccflux,n)
        centerflux, outflux = WDD_calc1D(q_e,Sigma_t,Sigma_s,alpha,group,mu,h,influx,scattering)
        print(centerflux)
        calc_ccflux[n] =  centerflux
        influx = outflux
        
    return calc_ccflux,outflux
 

if __name__ == "__main__":
    # Define parameters
    q_e = np.array([1.5, 0.0, 0.2]) # external source (by energy group)
    Sigma_t = np.array([0.5,0.8,1.0]) # total absorption cross section (by energy group)
    Sigma_s = np.array([[0.1,0.0,0.0],[0.3,0.1,0.1],[0.1,0.3,0.3]]) # scattering cross section (from energy g[row] to g'[column])
    Alpha = [0] # weighted diamond difference "weight"
    mu = [0.2,0.7] # angles
    H = [0.1] # mesh spacing
    lowbound = 0.0 # lower (left) boundary location
    upbound = 2.0 # upper (right) boundary location
    startflux = np.array([0.5,0,0]) # incoming flux from lower (left) boundary (by energy group)
    
    # Check that all arrays have equal numbers of groups
    if len(q_e) != len(Sigma_t) or len(q_e) != len(Sigma_s[0]) or len(q_e) != len(startflux):
        print('ERROR: All multigroup arrays must contain the same number of groups.')
        
    #An isotropic source means:
    q_e = q_e/(2*len(mu))
    
    # Run for various values of alpha (weighting parameter)
    for alpha in Alpha:
        print('\nalpha = ',alpha)
        
        # General plotting parameters
        fig = plt.figure(figsize=(12,10))
        subplotnum = 100*len(H)+10
        
        
        # Run for various values of h (mesh spacing)
        for h in H:
            print('\nMesh Spacing: ',h)
            ncells = int((upbound-lowbound)/h) # calculate the number of cells used in the mesh
            ngroups = len(q_e)
            
            # More specific plotting parameters
            subplotnum += 1
            ax = fig.add_subplot(subplotnum)
            # Title
            ax.text(0.5,0.8,'Center-cell Scalar Flux, h=%s' %(str(h)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            # Label
            ax.text(0.9,0.1,r'$\alpha$=%s' %(str(alpha)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
            ax.set_ylabel('$\psi_{i}$',rotation=0)
            plt.plot([0 for i in range(ncells)]) # set "zero" on plot
            
            mg_ccflux = np.zeros((ngroups,ncells,2*len(mu))) # an array to store the cell-center flux for each angle, mesh cell, and energy group
            mg_scalarflux = np.zeros((ngroups,ncells))
            
            # Perform multigroup iterations
            for group in range(ngroups):
                
                g_ccflux = np.zeros_like(mg_ccflux[group]) # an array to store the group's newly calculated cell-center fluxes
                
                # Prepare reflecting boundary conditions (construct an array to store values of flux at boundary)
                reflected = [0 for a in range(len(mu))] 
                
                # Iterate until the (group) scalar flux converges
                g_scalarflux = np.zeros(ncells)
                g_scalarfluxnorm2 = -1 # flag start of run
                f = 0 # limit for while loop
                while (g_scalarfluxnorm2 > 0.0001 or g_scalarfluxnorm2 == -1 ) and f<1000:
                    
                    calc_ccflux = np.zeros_like(g_ccflux) # an array to store the newly calculated cell-center fluxes
                    # Sweep (Left->Right), 1 sweep per angle
                    for a in range(len(mu)):
                        centerflux,outflux = Sweep(q_e,Sigma_t,Sigma_s,alpha,group,mu[a],a,h,ncells,mg_ccflux,startflux[group])
                        reflected[a] = outflux  # save the reflected flux to be used as the inputs in the other direction
                        calc_ccflux[:,a] += centerflux
                    # Sweep (Right->Left), 1 sweep per angle
                    for a in range(len(mu)):
                        a_index = a + len(mu)
                        centerflux,outflux = Sweep(q_e,Sigma_t,Sigma_s,alpha,group,-mu[a],a_index,h,ncells,mg_ccflux,reflected[a])
                        calc_ccflux[:,a_index] += centerflux
                    
                    g_scalarflux_prev = g_scalarflux
                    g_scalarflux = np.sum(calc_ccflux,axis=1)
                    
                    g_scalarfluxnorm2 = math.sqrt(sum((g_scalarflux - g_scalarflux_prev)**2))
                    
                    g_ccflux = calc_ccflux
                    f+=1
                    
                mg_ccflux[group] = g_ccflux
                mg_scalarflux[group] = g_scalarflux
                
                print('Group %s Scalar Flux:\n' %str(group),g_scalarflux)
                print('Group %s Norm2: ' %str(group),g_scalarfluxnorm2)
                
                plt.plot(mg_scalarflux[group])
                
        plt.show()
        #fig.savefig('C:/Users/Mitch/Documents/Cal/1 - 2016 Fall/NUCENG 255 - Numerical Simulation in Radiation Transport/Homework/HW_Drafts/HW_Figures/HW04_Prob3d/WDD_alpha%s.jpg' %(str(alpha)))
        