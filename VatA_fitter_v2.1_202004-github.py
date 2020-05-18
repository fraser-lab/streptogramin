"""
Written by: Jenna Pellegrino
Last updated: 20200428
Language: python 2
Execution: right now, you must edit the script with the file name(s) in order to execute
Purpose: 
    0) to import triplicate kinetics data in umol/min/mg
    1) to scatter plot all data points with error bars from calculated mean and stdev
    2) to pass all data points to Michaelis Menton eqn for curve fitting and to extract kinetics constants
    3) to export all those kinetics constants to a csv
    4) to use those kinetics constants in plotting best-fitting curve to the data
    5) to save plot as a figure
Limitations:
- hardcoded input
- hardcoded output figure name
- input data csv must be in the following format:
    {nothing} | Rep1   | Rep2   | Rep3     # these are mostly for you. first col must be x. all others are y replicates.
    [conc] uM | {name} | {name} | {name}   # {name} is the name of your compound
    {data}    | {data} | {data} | {data}   # the values start in this row
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy
from numpy import genfromtxt
from scipy.optimize import curve_fit
from scipy import stats


class Library(dict):
  def __init__ (self):
    self = dict()
  def add(self, name, marker, color, label):
    self[name]=Compound(name, marker, color, label)
    #self[name].add(name, marker, color, label)


class Compound():
  def __init__ (self, name="", marker="", color="", label=""):
    self.name = str(name)
    self.marker = str(marker)
    self.color = str(color)
    self.label = str(label)
  def add(self, name, marker, color, label):
    self.name = str(name)
    self.marker = str(marker)
    self.color = str(color)
    self.label = str(label)
  def __str__(self):
    return str([self.name, self.marker, self.color, self.label])

library = Library()
library.add("flopristin",'^','xkcd:amber',r'flopristin ($\bf{4}$)') #this matches color in mouse data
library.add("F1037",'s','xkcd:windows blue',r'$\bf{47}$')  #this matches color in mouse data


"""
REGULAR OL MICHAELIS MENTON EQUATION, x + 2 params
x = substrate concentration
Vmax = max enzyme velocity
Km = Km
"""
def michaelismenton(x, Vmax, Km):
    y = Vmax * (x) / (Km + x)
    return y

def run(inFile):
    df = pd.read_csv(file,header=1)
    
    """
    CONVERT 1 X-COL AND N Y-COLS TO 2 1D ARRAYS
    - one array will be all x values; other array will be all y values
    - this is necessary for feeding to the curve_fit function; it only takes one x[] and one y[]
    """
    x_all = []
    y_all = []
    for rep_col in range(1,len(df.columns)):
        #print rep
        x_all += df.iloc[:,0].tolist()
        y_all += df.iloc[:,rep_col].tolist()
    
    """
    GRAB COMPOUND NAME FOR FILE AND HEADER NAMING PURPOSES
    """
    compound_name = df.columns[1]
    
    """
    CALCULATE AVERAGE AND STANDARD DEVIATION
    - will use these to scatter plot and to add the error bars
    - will NOT use these for curve_fit
    """
    df[compound_name+"_avg"]=df.iloc[:,1:].mean(axis=1)
    df[compound_name+"_std"]=df.iloc[:,1:].std(axis=1)
    

    """
    FIT ALL X AND Y POINTS TO MICHAELIS-MENTON CURVE
    """
    try:
        popt, pcov = curve_fit(michaelismenton, x_all, y_all)
        ##popt = optimal values of the variables to sigfit: Vmax, Km (in that order)
        ##pcov = covariance of popt; the diagonals provide the variance of the parameter estimate
        stdev = np.sqrt(np.diagonal(pcov))
        ##stdev computes the one standard deviation errors on the sigfit variables: Vmax, Km (in that order)

    #except TypeError:
    #   print "TypeError for {name}: only integer scalar arrays can be converted to a scalar index".format(name=name)

    except ValueError:
        print "ValueError for {name}: array must not contain infs or NaNs".format(name=compound_name)
        print "Editing array to remove NaN values and passing back to michaelismenton for fitting."
        #you cannot pass curve_fit data with NaN; remove those nonexisting points from the dataset. you can do that with a bool.
        #https://stackoverflow.com/questions/33876226/scipy-curve-fit-fails-on-easy-linear-fit

        valid = ~(np.isnan(x_all) | np.isnan(y_all))
        valid #bool list that asks if value at index in x_ <OR> y_ is NaN. it will say "True" if either x_ or y_ is NaN at that index, hence |

        #https://stackoverflow.com/questions/14537223/remove-items-from-a-list-using-a-boolean-array
        #[d for (i) in xy["y_F1037"] if valid[i]==True]
        heyo = [y_all[d] for d in range(0,len(y_all)) if valid[d]==True]
        hexo = [x_all[d] for d in range(0,len(x_all)) if valid[d]==True]

        #pass these edited x and y values to michaelismenton
        popt, pcov = curve_fit(michaelismenton, hexo, heyo)
        stdev = np.sqrt(np.diagonal(pcov))

    """
    Calculate Kcat
    Add Vmax, Km, and Kcat to a column specified by the compound present in assay
    Add stdev for Vmax, Km, and Kcat
    """
    MW_VatA = 24287.01 ##MW of TEV cleaved VatA; 26267.12 is MW of VatA-6HIS-TEV (Daltons = g/mol)
    Kcat = popt[0]/60*1000*MW_VatA/10**6
    Kcat_stdev = stdev[0]/60*1000*MW_VatA/10**6
    out_df["{name}".format(name=compound_name)]=[popt[0],stdev[0],popt[1],stdev[1],Kcat,Kcat_stdev]
    print popt
    print stdev


    """
    Setting up min and max x in np.linspace.
    If you have to change anything, only change args of linspace.
    """
    fit_x = np.linspace(0, 0.5, 10000) #min x, max x, how many points in between 0 and 10
    fit_y = michaelismenton(fit_x, *popt)


    """
    Setting up the plot area (min and max y)
    Plotting your data and the fit curve
    - scatter average y data
    - error bar using the stdev
    - plot the curve using the calculated Michaelis Menton Vmax and Km params
    """
    plt.scatter(df["[conc] mM"],df[compound_name+"_avg"],
                color=library[compound_name].color,
                marker=library[compound_name].marker,
                label=library[compound_name].label,
               )
    plt.errorbar(df["[conc] mM"],df[compound_name+"_avg"],yerr=df[compound_name+"_std"],
                ls="none", #linestyle; 'none' so it won't connect the dots
                color=library[compound_name].color,
                marker='_', #adds marker for mean and makes it a horizontal line
                label=None
                )
    plt.plot(fit_x, fit_y, 
             '-', 
             color=library[compound_name].color
            )


    plt.xlabel("concentration (mM)", fontsize=13)
    plt.ylabel("velocity ($\mu$mol/min/mg)", fontsize=13)
    plt.legend(loc='lower right', fontsize=10, labelspacing=0.2)


if __name__ == "__main__":
    
    """
    SET UP PLOTTING AREA
    - all input files will share one plot space, so initiate this before opening any files
    """
    fig, ax = plt.subplots(figsize=(5,5)) #you can put figsize=(w,h) in the (); note that it stretches the axes
    handles, labels = ax.get_legend_handles_labels()
    
    """
    SET UP OUTPUT DF TO INCLUDE ALL MICHAELIS MENTON CONSTANTS AND VALUES
    """
    out_df = pd.DataFrame() ##for new exported csv containing calculated values
    out_df["variables"]=['Vmax','Vmax stdev','Km','Km stdev','Kcat', 'Kcat stdev']
    out_df["units"]=['umol/min/mg','umol/min/mg','mM','mM','1/s','1/s']

    """
    LIST THE INPUT FILES (CSV FORMAT)
    """
    in_Files = [
            'FileNameHere_1.csv',
            'FileNameHere_2.csv'
           ]
    
    """
    PASS EACH INPUT FILE ONE BY ONE TO THE RUN FUNCTION
    """
    for file in in_Files:
        run(file)
    
    plt.savefig("FigureNameHere.png", bbox_tight="inches",dpi=300)
    plt.show()
    



