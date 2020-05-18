
"""
Written by: Jenna Pellegrino
Last updated: 20200427
Language: python 2
Execution: right now, you must edit the script with the file name(s) in order to execute
Purpose (and features): 
    1) to plot raw mouse data using matplotlib, where'x' are categories and 'y' are numbers, but with using
       a scatter plot with number-line x values (1,2,3,etc) and then changing the labels to be the categories
    2) to add mean and stdev error bars
    3) to color by category
    4) to use handles to achieve the color and category name requirements I want
Limitations:
- hardcoded input(s)
- input needs to already have a mean column
- input needs to already have a stdev column
- hardcoded output png name
"""

import pandas as pd
import numpy as np
from numpy import median
import matplotlib.pyplot as plt


df = pd.read_csv('MouseData.csv')
df = df.rename(columns={"SA name": "SA_name"})
df = df.iloc[::-1] #reverses row order, but keeps the original row index numbers
df.reset_index(drop=True,inplace=True) #gives new row index numbers from 0 and drops the old, scrambled one; updates df
print df


"""
SWAPPING 24H AND 2H ROWS
- i want plot to start with 2H then 24H, left to right
- swapping these two rows in dataframe to accomplish that
- requires making a temp row to hold the copy
"""

twentyfourH, twoH = df.iloc[0].copy(), df.iloc[1].copy()
df.iloc[1], df.iloc[0] = twentyfourH, twoH

"""
ADDING IN COLUMN TO HOLD X AXIS LABELS
- will make this column by concatenating two columns, adding a ( between the strings, and
"""

controls=df.iloc[0:2, :] #isolating control rows (first two rows) from df so I can edit members of other rows
df=df.drop(df.index[0:2]) #dropping those control rows from df so I can edit members of other rows

controls['SA_name_Dose'] = controls['Dose (mg/kg)'].str.cat(controls['SA_name'], sep =" ") 
controls['SA_name_Dose'] = "+" + controls['SA_name_Dose'].astype(str)
controls['xticklabels'] = controls['SA_name_Dose'] #duplicating a column so there aren't NaN when I put controls back into edited df

df['SA_name_Dose'] = df['SA_name'].str.cat(df['Dose (mg/kg)'], sep =" (") 
df['SA_name_Dose'] = df['SA_name_Dose'].astype(str) + " mg/kg)"
df['xticklabels'] = df['Dose (mg/kg)'].astype(str) + " mg/kg"
df = pd.concat([controls,df], axis=0) #putting control rows back in; axis=0 adds them as new rows

"DECIDED TO GO WITH A SIMPLER XAXIS LABEL FOR THE CONTROLS"
df.loc[0,"xticklabels"]="2 h"
df.loc[1,"xticklabels"]="24 h"


"""
FOLLOWING LINES DROP THE +2 HOUR INFECTION CONTROL
- originally were gonna drop the 2H infection control, but we later decided to keep if
- to remove it, then uncomment the below out, BUT YOU HAVE TO UPDATE XLIM AND STUFF FOR THE FIGURE LATER
"""

flip=df.transpose() 
flip.columns=flip.iloc[-2] #setting SA_name_Dose, which is in row -2, as the column headers; i believe headers should be unique
flip=flip.reset_index() #essential to move the now-row index names to be a normal column so i can drop columns that don't start with "Point"


flip = flip[flip['index'].str.startswith('Point')] #removes all rows that don't start with 'Point' under the 'index' col
flip

fake_x=[]
counter = 0
y=[]
for item in df["SA_name_Dose"]:
    #print item
    #print flip[item]
    n = len(flip[item])
    for i in range(n):
        fake_x.append(counter)
    counter+=1
    y.extend(flip[item])

print "fake_x", len(fake_x)
print "y", len(y)
# the lengths of those two should be the same


"""
IN ALL THE BELOW BLOCKS, I am changing xtick labels from # mg/kg to just be #.
"""
df["xticklabels_v2"] = df["Dose (mg/kg)"]
df.loc[0,"xticklabels_v2"]="2 h"
df.loc[1,"xticklabels_v2"]="24 h"


"""
https://stackoverflow.com/questions/26139423/plot-different-color-for-different-categorical-levels-using-matplotlib
"""
xkcd_colors = {'47':"xkcd:windows blue",'4':"xkcd:amber",'Infection Control':"xkcd:dusty purple"}

fig, ax = plt.subplots() #you can put figsize=(w,h) in the (); note that it stretches the axes
handles, labels = ax.get_legend_handles_labels()

"""
HARD-CODING THE PLOTTING IN ORDER TO CONTROL THE LEGEND
"""
handle_47 = plt.scatter(fake_x[35:],y[35:],label=r'$\bf{47}$',c="xkcd:windows blue",s=14)
handle_flopristin = plt.scatter(fake_x[10:35],y[10:35],label=r"flopristin ($\bf{4}$)",c="xkcd:amber",s=14)
handle_controls = plt.scatter(fake_x[0:10],y[0:10],label="no treatment",c="xkcd:dusty purple",s=14)

"""
ADDING ERROR BARS TO PLOT
modified from the below line of code, which came from the link below
- plt.errorbar(range(len(mean)), mean, yerr=std)
- https://stackoverflow.com/questions/44506136/plotting-errorbars-on-top-of-swarmplot
"""
handle_errorbar = plt.errorbar(range(len(df["Mean log10 CFU/Thigh"])), 
             df["Mean log10 CFU/Thigh"], 
             yerr=df["Standard Deviation"],
             ls="none", ##https://stackoverflow.com/questions/40020274/pyplot-errorbar-keeps-connecting-my-points-with-lines/40020492
             capsize=8,
             ecolor='black',
             marker='_', #adds marker for mean and makes it a horizontal line
             mec = 'black', #markeredgecolor; the caps and marker are edges
             mew = 1.5, #markeredgewidth
             ms=20, #markersize,
             zorder=0, #should force this to go in the background
             label=r'mean log$_{10}$ CFU/thigh' #can edit label name this way
            )

"""
CHANGING YLIM
- want min to be 0
- max can be whatever, but let's do 8
"""
ax.set_ylim(3,8)
 
"""
ADDING NAMES OF XAXIS AND YAXIS LABELS
- there are more than one way to do this
- included options below just to see them later
- plt.xlabel and plt.ylabel seems most straight-forward
"""
plt.xlabel('')
plt.ylabel(r'log$_{10}$ (CFU/thigh)', fontsize=12)

"""
FORMATING X AND Y TICK LABELS
"""
tick_locations = np.array([0,1,2,3,4,5,6,7,8,9,10,11])
ax.set_xticks(tick_locations)

myxticklabels=df["xticklabels_v2"].tolist()
ax.set_xticklabels(myxticklabels,
                   rotation=45,
                   #rotation=70,
                   #fontsize=10,
                   va="top",ha="center")

"""
REMOVE BORDER FROM TOP AND RIGHT OF RIGURE
These are called "spines"
"""
# Hide the top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.show()

fig.savefig("FileNameHere.png",bbox_inches="tight",dpi=300) #including bbox_inches="tight" ensures xtick labels aren't cut off



