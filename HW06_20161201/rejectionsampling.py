'''
Created on Nov 27, 2016
@author: Mitch

Program to perform rejection sampling in order to calculate pi
'''
import random
import distfunc as df

def scoresample(x,y,pdf):
    # Function to score all y < PDF(x)
    if y < pdf(x):
        return 1
    else:
        return 0

    
if __name__ == "__main__":
    
    pi_approx = 3.14159
    sampleset = [10,100,1000,10000]
    PDF = df.deriv_atan
    
    for samples in sampleset:
        score = 0
        for s in range(1,samples+1):
            x = random.random()
            y = random.random()
            # Accumulate total score
            score += scoresample(x,y,PDF)
            # Score ratio = accepted/(accepted+rejected); score should approach the integral of the PDF for many samples
            scoreratio = score/s
            
            picalc = 4*scoreratio      # multiply by 4 since the integrals of the given pdfs each equal pi/4
            relerr = abs(picalc-pi_approx)/s
        print(picalc,relerr)