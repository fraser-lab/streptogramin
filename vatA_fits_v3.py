#=================================================================================
# Programmer: Jenna Pellegrino
#
# edited from sigmoid_v6_labelsOFF.py
#
# for the sigfit equation:
# https://www.graphpad.com/guides/prism/7/curve-fitting/index.htm?reg_allosteric_enzyme.htm#
#
# for the linear regression:
# https://plot.ly/matplotlib/linear-fits/
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
#
# Date Made:	17 May 2019
# Updated:      11 June 2019
# Changes:      added more explanation; exports a new csv with calculated values and stdev for
#               Vmax, h, Km ("Khalf"), and Kcat values and stdev
# Language:     python
# Purpose:      to take replicate vatA assay data and plot it to the Allosteric sigmoidal curve;
#               to calculate Vmax, h, and Km values (along with their stdev)
#
# Requirements: csv, think you can only run things that are set up exactly like 'Prism_vatA_fitme.csv'
# Notes: not set up to run on just anything. will need to specify columns to use for y_unique and y_unique_std
# Execution:    python vatA_fits_v2.py Prism_vatA_fitme.csv
#=============================================================================================

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

print sys.argv[1]

"""
Allosteric sigmoidal, x + 3 params:
x = substrate concentration
Vmax = max enzyme velocity
h = hill slope
Km = Km
"""
def sigfit(x, Vmax, h, Km):
    y = Vmax * (x**h) / ((Km**h) + (x**h))
    return y

"""
Opening and reading the csv, transposing columns, and doing log10 of y values
"""
def run(inFile):
    df = pd.read_csv(inFile, header=1)

    """
    Take average of replicate data. Add new column.
    """
    df["VM2"] = df[["VM2_1","VM2_2"]].mean(axis=1)
    df["flopristin"] = df[["flopristin_1","flopristin_2"]].mean(axis=1)
    df["0224"] = df[["0224_1","0224_2"]].mean(axis=1)
    df["F0224"] = df[["F0224_1","F0224_2"]].mean(axis=1)
    df["F1037"] = df[["F1037_1","F1037_2"]].mean(axis=1)

    """
    Get standard deviation (stdev) of the replicates. Add new column.
    """
    df["VM2_std"] = df[['VM2_1','VM2_2']].std(axis=1)
    df["flopristin_std"] = df[["flopristin_1","flopristin_2"]].std(axis=1)
    df["0224_std"] = df[["0224_1","0224_2"]].std(axis=1)
    df["F0224_std"] = df[["F0224_1","F0224_2"]].std(axis=1)
    df["F1037_std"] = df[["F1037_1","F1037_2"]].std(axis=1)
    # print df

    """
    Grab headers from averaged replicate data (y_unique); one header per compound.
    Grab headers from stdev of replicate data (y_unique_std); one header per compound.
    Iterate through the columns of y data by calling their headers (y_unique).
    Feed each y column of data with x data to sigfit function.
    Using counters for the marker shape and color.
    """
    column_headers = list(df)
    y_unique = column_headers[11:-5]
    y_unique_std = column_headers[-5:]
    #print y_unique
    #print y_unique_std
    #print df

    x_data = df["[drug] mM"]
    marker = 0
    markers = ['o','^','s','P','D']
    color = 0
    colors = ['teal', 'orange', 'r', 'm', 'limegreen']
    label = 0
    labels = ['VM2', r'flopristin ($\bf{4}$)', '40q', '46', '47'] #$\ bold font {4} $
    #plt.figure(figsize=(5, 9)) ##I added this later to resize the image, since we wanted a tall rectangle
    plt.figure(figsize=(5, 5)) ##I added this later to resize the image, since we wanted a tall rectangle

    out_df = pd.DataFrame() ##for new exported csv containing calculated values
    out_df["variables"]=['Vmax','Vmax stdev','h','h stdev','Km','Km stdev','Kcat', 'Kcat stdev'] #had to have empty cell so col length was the same across the board
    for i in range(len(y_unique)):
        y_data = df[y_unique[i]]
        y_data_std = df[y_unique_std[i]]
        # print y_data
        # print y_data_std

        popt, pcov = curve_fit(sigfit, x_data, y_data, p0=[15, 1, 0.02], bounds=([5.,0.5,0.0001],[50.,2.,0.1])) ##guesses for Vmax, h, Km
        ##popt = optimal values of the variables to sigfit: Vmax, h, Km (in that order)
        ##pcov = covariance of popt; the diagonals provide the variance of the parameter estimate
        stdev = np.sqrt(np.diagonal(pcov))
        ##stdev computes the one standard deviation errors on the sigfit variables: Vmax, h, Km (in that order)

        """
        Calculate Kcat
        Add Vmax, h, Km, and Kcat to a column specified by the compound present in assay
        Add stdev for Vmax, h, Km, and Kcat
        """
        MW_vatA = 26267.12 ##molecular weight of vatA-6HIS-TEV (Daltons = g/mol)
        Kcat = popt[0]/60*1000*MW_vatA/10**6
        Kcat_stdev = stdev[0]/60*1000*MW_vatA/10**6
        out_df["{compound}".format(compound=y_unique[i])]=[popt[0],stdev[0],popt[1],stdev[1],popt[2],stdev[2],Kcat,Kcat_stdev]
        print popt
        print stdev

        """
        Passing things to the sigfit function.
        Setting up min and max x in np.linspace.
        If you have to change anything, only change args of linspace.
        """
        x = np.linspace(0, 0.35, 10000) #min x, max x, how many points in between 0 and 10
        y = sigfit(x, *popt)

        """
        Setting up the plot area (min and max y)
        Plotting your data and the fit curve
        """
        plt.plot(x_data, y_data, markers[marker], color=colors[color], label=labels[label])
        plt.errorbar(x_data, y_data, y_data_std, color=colors[color], linestyle='None', label=None)
        plt.plot(x, y, '-', color=colors[color])
        plt.xlabel("concentration [mM]", fontsize=13)
        plt.ylabel("velocity ($\mu$mol/min/mg)", fontsize=13)
        plt.legend(loc='best', fontsize=10)

        marker+=1
        color+=1
        label+=1

    out_df.to_csv("Prism_vatA_fitme_CALCULATED_KINETICS_VARIABLES_20190621.csv",index=False)
    #plt.savefig("vatA_fits.png", dpi=300)
    plt.savefig("vatA_fits_20190621.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    run(sys.argv[1])
