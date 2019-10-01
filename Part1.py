import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate

    # Constants
c1 = 0.0380998 # nm^2 eV h-^2/2m
c2 = 1.43996 # nm eV e^2/4pieps0
r0 = 0.0529177 # nm
h  = 6.62606896e-34 # J s
c  = 299792458 # m/s

#calculate the coulomb force
def coulomb_force(r,alpha):
    #calculate the coulomb force
    coulomb_force=-c2*(1/r**2)*(r/r0)**alpha
    return coulomb_force

#works out electric potential
def electric_potential(alpha,r):
    #uses scipy quads to integrate the force
    potential=integrate.quad(coulomb_force,r,np.inf,alpha)
    return potential[0]

#plots the results
def plot(xmin,xmax,step,alpha):
    #builds the range of r
    rrange=np.arange(xmin,xmax+step,step)
    vrange=[]
    for r in rrange:
        #builds the range of v
        v=electric_potential(alpha,r)
        vrange.append(v)
    plt.plot(rrange,vrange)
    plt.xlabel("r in nm")
    plt.ylabel("V in eV")
    plt.title("Variation of Electric Potential with distance from the origin")
    plt.show()

def analytic_electric_potential(r,alpha):
    #uses the analytic method to find v for a given r
    v_anal=c2*1/((alpha-1)*r0**alpha)*r**(alpha-1)
    return v_anal

def compare_anal_to_numeric(rmin,rmax,alpha,step):
    rrange=np.arange(rmin,rmax+step,step)
    diff_range=[]
    for r in rrange:
        v_numeric=electric_potential(alpha,r)
        v_anal=analytic_electric_potential(r,alpha)
        diff=abs(v_anal-v_numeric)
        diff_range.append(diff)
    return max(diff_range)
        
    
rmin=0.01
rmax=1
alpha=0.01
step=0.01

plot(rmin,rmax,step,alpha)
print(compare_anal_to_numeric(rmin,rmax,alpha,step))



    


    