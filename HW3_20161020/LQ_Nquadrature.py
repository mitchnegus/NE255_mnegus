'''
Created on Oct 19, 2016

@author: Mitch
'''
import math

#ew_file = open(eqweight_filename)
def opensplit(filename):
    infile = open(filename,'r')
    instring = infile.read().strip()
    inlist = instring.split('\n')
    return inlist

def select_SN(inlist,N):
    for i in range(len(inlist)):
        if inlist[i][0] == 'S' and int(inlist[i][1:]) == N:
            SNsect = inlist[i+1:i+int(N/2)+1]
            break
    return SNsect
    
def read_weightvals(filename,N):
    wv_list = opensplit(filename)
    n_mu_w = select_SN(wv_list,N)
    for i in range(len(n_mu_w)):
        n_mu_w[i] = n_mu_w[i].strip("\t")
        n_mu_w[i] = n_mu_w[i].split("\t")
    return n_mu_w

def read_equalweights(filename,N):
    ew_list = opensplit(filename)
    weightshape = select_SN(ew_list,N)
    weightcount = {}
    for i in range(len(weightshape)):
        weightshape[i] = weightshape[i].split()
        for j in range(len(weightshape[i])):
            weightshape[i][j] = int(weightshape[i][j])
            if weightshape[i][j] not in weightcount:
                weightcount[weightshape[i][j]]=1
            else:
                weightcount[weightshape[i][j]]+=1
    return(weightcount)


def LQ_N_integrate(N,n_w,weightcounts):
        
    wsum = 0
    for n in n_w:
        factor = weightcounts[n]
        weight = n_w[n]
        wsum += factor*weight
    # Define the normalization factor
    A = math.pi/2
    integral = 8*A*wsum
    
    return integral

'''
def LQ_N_integrate_wFunc(N,n_w,weightcounts,func):
    # IN PROGRESS
    wsum = 0
    for n in n_w:
        factor = weightcounts[n]
        weight = n_w[n]
        wsum += factor*weight
    # Define the normalization factor
    A = math.pi/2
    integral = 8*A*wsum
    
    return integral
'''

weightval_filename = 'L&M_LQ_N_QuadSets.txt'
eqweight_filename = 'L&M_LQ_N_EqualWeights.txt'

for N in [4,6,8,12,16]:
    print('\nS',N)
    n_mu_w = read_weightvals(weightval_filename,N)
    weightcounts = read_equalweights(eqweight_filename,N)
    
    n_w = {}
    for i in n_mu_w:
        if len(i) == 3:
            n_w[int(i[0])]=float(i[2])
    integral = LQ_N_integrate(N,n_w,weightcounts)
    print('Integral = '+str(integral/math.pi)+'*pi')
