'''
Created on Nov 27, 2016
@author: Mitch

Program to perform rejection sampling in order to calculate pi
'''
import random
import math
import distfunc as df

def scoresample(x,y,pdf):
    # Function to score all y < PDF(x)
    if y < pdf(x):
        return 1
    else:
        return 0

    
if __name__ == "__main__":
    
    pival = 3.14159
    sampleset = [10,100,1000,10000]
    PDF = df.deriv_atan
    
    for samples in sampleset:
        score = 0
        ysum = 0
        ysumsq = 0
        for s in range(1,samples+1):
            x = random.random()
            y = random.random()
            # Accumulate total score
            thisscore = scoresample(x,y,PDF)
            score += thisscore
            if thisscore:
                ysum += y
                ysumsq += y**2
            # Score ratio = accepted/(accepted+rejected); score should approach the integral of the PDF for many samples
            scoreratio = score/s
            
        relerr = math.sqrt(abs((1/samples**2)*ysumsq/(pival)**2-(1/samples)))
        pi_calc = 4*scoreratio
        relerr = 4*relerr
        print(pi_calc,relerr)    
        
        