#=================================================================================
# Programmer: Jenna Pellegrino
#
# inspired by this code:
# https://gist.github.com/andrewgiessel/5684769
#
# for the sigfit_hill equation:
# https://www.graphpad.com/guides/prism/7/curve-fitting/index.htm?reg_dr_inhibit_variable_2.htm
#
# Date Made:	16 October 2018
# Last Updated: 19 February 2019
# Language:     python 2
# Purpose:      Takes triplicate data. Plots all data points.
#               Fits data to a sigmoid. Plots the fit sigmoid with vertical and horizontal
#               dropdown lines to help you see where the calculated IC50 is located.
# Requirements: csv, header row, X vals in col 0 and Y vals in col 1-3
#               [drug] uM rep1 rep2 rep3 <<< NO SPACES for rep#!!!!!
#               # blank blank blank <<< need a # where I have #, even though it's the blank
#               0 # # # <<< first line with concentrations
#               # # # #
# Notes:        (1) the -l and -o are optional. if you're going to use one of
#               them, use -l. -o will duplicate -l. without both, -l and -o will
#               be input file name without the extension.
#               (2) if you do use -l and/or -o, the number of args must be the
#               same across -i -l and -o
# Execution:    python sigmoid_v6.py -i 0224.csv F0224.csv
#               python sigmoid_v6.py -i 0224.csv F0224.csv -l 0224 F0224
#               python sigmoid_v6.py -i 0224.csv F0224.csv -o 0224 F0224 -l 0224_fig F0224_fig
#=============================================================================================

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from numpy import genfromtxt
from scipy.optimize import curve_fit

import argparse
parser = argparse.ArgumentParser(description="no desc")
parser.add_argument('-i', nargs="+", help='input file in csv') ##nargs + takes 1 or more
parser.add_argument('-o', nargs="*", type=str, help='output figure png') ##nargs * take 0 or more
parser.add_argument('-l', nargs="*", type=str, help='label of drug for figure (eg F0224)') ##nargs * take 0 or more
args = parser.parse_args()

"""
Dose-response inhibition curve w hill slope, x + 4 params:
x = concentration
t = top of slope; essentially y when x == 0
b = bottom of slope; essentially y when x == infitinity
i = IC50; as the last of the variables, it will be popt[-1]
"""

def sigfit(x,t,b,i):
    y = b + (t-b)/(1 + (x/i))
    return y


"""
Opening and reading the csv, transposing columns, and doing log10 of y values
"""
def run(inFile, outFile, label):
    df = pd.read_csv(inFile, header=0)

    """
    Subtract out the blanks; create new column
    """
    df["rep1-bkg"] = df["rep1"][1:] - df["rep1"][0]
    df["rep2-bkg"] = df["rep2"][1:] - df["rep2"][0]
    df["rep3-bkg"] = df["rep3"][1:] - df["rep3"][0]
    #print df

    """
    Extract all x and y data. Feed to sigfit_hill function.
    """
    all_x_data = np.concatenate([df["[drug] uM"][1:], df["[drug] uM"][1:], df["[drug] uM"][1:]])
    all_y_data = np.concatenate([df["rep1-bkg"][1:], df["rep2-bkg"][1:], df["rep3-bkg"][1:]])

    popt, pcov = curve_fit(sigfit, all_x_data, all_y_data)
    IC50 = popt[-1]
    ##popt = optimal values of the variables to sigfit_hill: t, b, h, i (in that order)
    ##pcov = covariance of popt; the diagonals provide the variance of the parameter estimate
    stdev = np.sqrt(np.diagonal(pcov))
    print "IC50 = ", IC50, " +- ", stdev[-1]
    print stdev

    """
    Passing things to the sigfit_hill function
    Setting up min and max x in np.linspace
    If you have to change anything, only change args of linspace
    """
    x = np.linspace(0, 36, 10000) #min x, max x, how many points in between 0 and 10
    y = sigfit(x, *popt)

    """
    Setting up the plot area (min and max y)
    Plotting your data and the fit curve
    """
    plt.plot(all_x_data, all_y_data, 'o')
    plt.plot(x, y, '-')
    plt.xscale("log")
    plt.xlabel("log10 of concentration ("+u"\u03bc"+"M)")
    plt.ylabel("fluorescence at 535 nm")


    """
    Adding dashed vertical and horizontal lines to show where the IC50 (popt[-1]) is
    """
    vertical_x = [popt[-1],popt[-1]]
    vertical_y = [0, sigfit(popt[-1],*popt)]
    plt.plot(vertical_x, vertical_y, "k--") #rounding IC50 to 2 decimal places

    horizontal_x = [0, popt[-1]]
    horizontal_y = [sigfit(popt[-1],*popt), sigfit(popt[-1],*popt)]
    plt.plot(horizontal_x, horizontal_y, "k--")

    plt.legend(loc='best')

    rounded_IC50 = str(round(float(IC50),2))
    rounded_stdev = str(round(stdev[-1], 4))
    informative_title = str(outFile+": IC50 = "+rounded_IC50+" +/- "+rounded_stdev)
    plt.title(informative_title)

    plt.savefig(outFile+"_IC50="+str(round(float(IC50),3))+"uM_stdev="+str(rounded_stdev)+"_allData_labelsOFF.png", dpi=300)
    plt.show()

"""
Grouping the -i and -o names so I can run the nth -i and save it with the nth -o name.
Passing both nth names together to run()
"""
if __name__ == "__main__":
    input = []
    for i in vars(args).get("i"): ##need vars() bc you wanna iterate over a namespace object
        input.append(i)
        #run(input)

    label = []
    if vars(args).get("l") == None:
        for i in vars(args).get("i"): ##need vars() bc you wanna iterate over a namespace object
            label.append(i[0:-4]) ##the label will be input file name without the extension
    else:
        for l in vars(args).get("l"):
            label.append(l)

    output = []
    if vars(args).get("o") == None:
        output = label ##note: these point to the same array now
    else:
        for o in vars(args).get("o"):
            output.append(o)
    print input, output, label


    for i in range(len(input)):
        print input[i], output[i]
        run(input[i], output[i], label[i])
