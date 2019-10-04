import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate
from numpy import linalg as LA
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from tkinter import *
from tkinter import ttk
from tkinter import messagebox

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

 #       
def find_H(N,alpha,step):
    delta=1/dict_of_vals['step']**2*diags([1, -2, 1], [-1, 0, 1], shape=(N,N))
    V_ii=[]
    for n in range(1,dict_of_vals['N']+1):
        V_ii.append(electric_potential(alpha,n*step))
    V=diags(V_ii)
    H=-c1*delta+V
    return H
   
def two_lowest_eigens(matrix):
    vals=np.real(eigs(matrix,k=2,which='SR')[0])
    return vals

def submit(root,entries):
    global dict_of_vals
    vals={}
    errpoints=[]
    error=False
    for items in entries.items():
        try:
            if items[0]=="N":
                vals[items[0]]=int(items[1].get())
            else:
                vals[items[0]]=float(items[1].get())
        except:
            error=True
            errpoints.append(items[0])
        finally:
            items[1].delete(0,'end')
    if error:
        error_message=""
        for e in range(0,len(errpoints)):
            if e==0:
                error_message+=(" "+errpoints[e])
            else:
                error_message+=(", "+errpoints[e])
        messagebox.showinfo("Entry Error","Please check entry for"+error_message)
    dict_of_vals=vals


def alpnumbutton(frame):
    diff=compare_anal_to_numeric(dict_of_vals['rmin'],dict_of_vals['rmax'],dict_of_vals['alpha'],dict_of_vals['step'])
    string="The Maximum difference between the analytical and numerical values is:\n "+str(diff)
    ttk.Label(frame,text=string,justify="center").pack()

def eigenbutton(frame):
    H=find_H(dict_of_vals['N'],dict_of_vals['alpha'],dict_of_vals['step'])
    vals=two_lowest_eigens(H)
    string="The two lowest eigenvalues are:\n"+str(vals[0])+"ev and "+str(vals[1])+"ev"
    ttk.Label(frame,text=string,justify="center").pack()


def main():

    root=Tk()
    notebook=ttk.Notebook()
    notebook.pack()
    vars_frame=ttk.Frame(notebook)
    notebook.add(vars_frame,text="Enter Variables")
    ttk.Label(vars_frame,text="Minimum r: ").grid(row=0,column=0,padx=5,sticky="sw")
    ttk.Label(vars_frame,text="Maximum r: ").grid(row=1,column=0,padx=5,sticky="sw")
    ttk.Label(vars_frame,text="α: ").grid(row=2,column=0,padx=5,sticky="sw")
    ttk.Label(vars_frame,text="step: ").grid(row=3,column=0,padx=5,sticky="sw")
    ttk.Label(vars_frame,text="N: ").grid(row=4,column=0,padx=5,sticky="sw")
    entries={"rmin":ttk.Entry(vars_frame),"rmax":ttk.Entry(vars_frame),"alpha":ttk.Entry(vars_frame),"step":ttk.Entry(vars_frame),"N":ttk.Entry(vars_frame)}
    ttk.Button(vars_frame,text="submit",command= lambda: submit(root,entries)).grid(row=5,column=0)
    i=0
    for items in entries.items():
        items[1].grid(row=i,column=1)
        i+=1

    graph_frame=ttk.Frame(notebook)
    ttk.Button(graph_frame,text="draw graph",command= lambda: plot(dict_of_vals['rmin'],dict_of_vals['rmax'],dict_of_vals['step'],dict_of_vals['alpha'])).pack()
    notebook.add(graph_frame,text="plot graph")
    compare_frame=ttk.Frame(notebook)
    ttk.Button(compare_frame,text="compare analytical value to numerical value",command= lambda: alpnumbutton(compare_frame)).pack()
    notebook.add(compare_frame,text="compare analytical to numerical")
    eigen_frame=ttk.Frame(notebook)
    ttk.Button(eigen_frame,text="Find Two Lowest Eigenvalues (first two energy levels)",command= lambda: eigenbutton(eigen_frame)).pack()
    notebook.add(eigen_frame,text="Eigenvalues")
    root.mainloop()


main()
    


    