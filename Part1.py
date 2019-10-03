import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as integrate
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
        
    
# rmin=0.01
# rmax=1
# alpha=0.01
# step=0.01

# plot(rmin,rmax,step,alpha)
# print(compare_anal_to_numeric(rmin,rmax,alpha,step))


def submit(root,entries):
    global dict_of_vals
    vals={}
    errpoints=[]
    error=False
    for items in entries.items():
        try:
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
    entries={"rmin":ttk.Entry(vars_frame),"rmax":ttk.Entry(vars_frame),"alpha":ttk.Entry(vars_frame),"step":ttk.Entry(vars_frame)}
    ttk.Button(vars_frame,text="submit",command= lambda: submit(root,entries)).grid(row=4,column=0)
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
    root.mainloop() 

main()
    


    