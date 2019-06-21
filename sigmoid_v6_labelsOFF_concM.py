#=================================================================================
# Programmer: Jenna Pellegrino
#
# inspired by this code:
# https://gist.github.com/andrewgiessel/5684769
#
# for the sigfit_4 equation:
# https://www.graphpad.com/guides/prism/7/curve-fitting/index.htm?reg_dr_inhibit_variable_2.htm
#
# Date Made:	16 October 2018
# Updated:      19 February 2019
# Language:     python
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
h = hill slope
i = IC50; as the last of the variables, it will be popt[-1]
"""
def sigfit_4(x, t, b, h, i):
    y = b + (t-b)/(1 + ((x**h)/(i**h)))
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

    df["x_data_in_M"] = df["[drug] uM"].multiply(0.000001)
    print df

    """
    Extract all x and y data. Feed to sigfit_4 functionself.
    """
    all_x_data_uM = np.concatenate([df["[drug] uM"][1:], df["[drug] uM"][1:], df["[drug] uM"][1:]])
    all_x_data_M = np.concatenate([df["x_data_in_M"][1:], df["x_data_in_M"][1:], df["x_data_in_M"][1:]])
    #print all_x_data_uM
    #print all_x_data_M

    all_y_data = np.concatenate([df["rep1-bkg"][1:], df["rep2-bkg"][1:], df["rep3-bkg"][1:]])
    #all_y_data_log10 = np.log10(all_y_data)

    #popt, pcov = curve_fit(sigfit_4, all_x_data_M, all_y_data)
    #popt, pcov = curve_fit(sigfit_4, all_x_data_M, all_y_data, p0=[1300, 2300, 1, 0.00000013]) ##for VM2
    popt, pcov = curve_fit(sigfit_4, all_x_data_M, all_y_data, p0=[5000, 800, 1, 0.0000003]) ##for all
    IC50 = popt[-1]
    ##popt = optimal values of the variables to sigfit_4: t, b, h, i (in that order)
    ##pcov = covariance of popt; the diagonals provide the variance of the parameter estimate
    stdev = np.sqrt(np.diagonal(pcov)) ##determines one standard deviation errors
    print "IC50 = ", IC50, " +- ", stdev[-1]
    print popt
    print pcov
    print stdev



    """
    Passing things to the sigfit_4 function
    Setting up min and max x in np.linspace
    If you have to change anything, only change args of linspace
    """
    x = np.linspace(0, 0.0001, 10000) #min x, max x, how many points in between 0 and 10
    y = sigfit_4(x, *popt)

    """
    Setting up the plot area (min and max y)
    Plotting your data and the fit curve
    """
    #plt.plot(all_x_data_uM, all_y_data, 'o', label=label)
    #plt.plot(all_x_data_uM, all_y_data, 'o')
    plt.plot(all_x_data_M, all_y_data, 'o')
    #plt.plot(x, y, '-', label='fit')
    plt.plot(x, y, '-')
    plt.xscale("log")
    #plt.xlabel("log"+r"$_1$"+r"$_0$"+" of concentration ("+u"\u03bc"+"M)")
    #plt.xlabel("concentration ("+u"\u03bc"+"M)")
    plt.xlabel("log"+r"$_1$"+r"$_0$"+" of concentration (M)")
    #plt.xlabel("concentration (M)")
    plt.ylabel("fluorescence at 535 nm (AU)")


    """
    Adding dashed vertical and horizontal lines to show where the IC50 (popt[-1]) is
    """
    vertical_x = [popt[-1],popt[-1]]
    vertical_y = [0, sigfit_4(popt[-1],*popt)]
    #plt.plot(vertical_x, vertical_y, "k--", label=round(float(IC50),2)) #rounding IC50 to 2 decimal places
    plt.plot(vertical_x, vertical_y, "k--") #rounding IC50 to 2 decimal places

    horizontal_x = [0, popt[-1]]
    horizontal_y = [sigfit_4(popt[-1],*popt), sigfit_4(popt[-1],*popt)]
    plt.plot(horizontal_x, horizontal_y, "k--")

    plt.legend(loc='best')

    rounded_IC50 = str(round(float(IC50),2))
    rounded_stdev = str(round(stdev[-1], 4))
    #informative_title = str(outFile+": IC50 = "+rounded_IC50+" +/- "+rounded_stdev)
    informative_title = str(outFile+": IC50 = "+str(IC50)+" +/- "+str(stdev[-1]))
    plt.title(informative_title)


    #title = raw_input("Give fig a title: ")
    #type(title)
    #plt.savefig(outFile+"_IC50="+str(round(float(IC50),3))+"_stdev="+str(rounded_stdev)+"_allData_labelsOFF.png", dpi=300)
    #plt.savefig(outFile+"_allData_labelsOFF_conc_in_log10_M.png", dpi=300)

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
