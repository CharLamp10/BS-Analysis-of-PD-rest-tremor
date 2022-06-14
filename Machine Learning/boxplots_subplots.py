"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the boxplots of each feature are extracted, as shown in figure 2
of ( https://doi.org/10.36227/techrxiv.19589728.v1). 

-----------------------------------------------------------------------------------------------------------------
Authors: Charalampos Lamprou & Ioannis Ziogas
Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
-----------------------------------------------------------------------------------------------------------------

"""

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import seaborn as sb
import matplotlib.patches as mpatches

sb.set_theme()
save_plots = True
init_path = os.getcwd()

""" HAT/LAT"""
features = ["bicEn12h8","bspecEn22h8","totalArea"]
plot_labels = ["BicEnt1","BspecEnt2","TotalArea"]
dataset = pd.read_csv('characteristics_MedDbsOnOffFirst4.csv')

clas = dataset['class']
clas.where(clas == 0,'HAT',inplace = True)
clas.where(clas == 'HAT','LAT',inplace = True)
fig = plt.figure(figsize = (20,15),dpi = 350)
    
ax = fig.add_subplot(2, 1, 1)
sb.boxplot(x="class", y=features[0], data=dataset)
ax.set_ylabel(plot_labels[0],fontsize = 25)
ax.set_xlabel("class",fontsize = 25)
ax.tick_params(labelsize = 25)
    
ax = fig.add_subplot(2, 2, 3)
sb.boxplot(x="class", y=features[1], data=dataset)
ax.set_ylabel(plot_labels[1],fontsize = 25)
ax.set_xlabel("class",fontsize = 25)
ax.tick_params(labelsize = 25)
    
ax = fig.add_subplot(2, 2, 4)
sb.boxplot(x="class", y=features[2], data=dataset)
ax.set_ylabel(plot_labels[2],fontsize = 25)
ax.set_xlabel("class",fontsize = 25)
ax.tick_params(labelsize = 25)

if save_plots:
    sel_folder = os.path.join(init_path,'plots boxplots')
    if not os.path.exists(sel_folder):
        os.mkdir(sel_folder)
    fname = os.path.join(sel_folder,'HAT_LAT_boxplot.png')
    fig.savefig(fname,bbox_inches='tight')

""" Med On/Off"""
features = ["bicEn12h8","bspecEn12h8","totalArea","power46"]
plot_labels = ["BicEnt1","BspecEnt1","TotalArea","Power46"]
dataset = pd.read_csv('characteristics_MedDbsOnOffFirst4.csv')

clas = dataset['Med']
clas.where(clas == 0,'Med On',inplace = True)
clas.where(clas == 'Med On','Med Off',inplace = True)
fig = plt.figure(figsize = (20,15),dpi = 350)
    
ax = fig.add_subplot(2, 2, 1)
sb.boxplot(x="Med", y=features[0], data=dataset)
ax.set_ylabel(plot_labels[0],fontsize = 25)
ax.set_xlabel("Med",fontsize = 25)
ax.tick_params(labelsize = 25)
    
ax = fig.add_subplot(2, 2, 2)
sb.boxplot(x="Med", y=features[1], data=dataset)
ax.set_ylabel(plot_labels[1],fontsize = 25)
ax.set_xlabel("Med",fontsize = 25)
ax.tick_params(labelsize = 25)
    
ax = fig.add_subplot(2, 2, 3)
sb.boxplot(x="Med", y=features[2], data=dataset)
ax.set_ylabel(plot_labels[2],fontsize = 25)
ax.set_xlabel("Med",fontsize = 25)
ax.tick_params(labelsize = 25)

ax = fig.add_subplot(2, 2, 4)
sb.boxplot(x="Med", y=features[3], data=dataset)
ax.set_ylabel(plot_labels[3],fontsize = 25)
ax.set_xlabel("Med",fontsize = 25)
ax.tick_params(labelsize = 25)

if save_plots:
    sel_folder = os.path.join(init_path,'plots boxplots')
    if not os.path.exists(sel_folder):
        os.mkdir(sel_folder)
    fname = os.path.join(sel_folder,'med_on_off_boxplots.png')
    fig.savefig(fname,bbox_inches='tight')

""" DBS On/Off"""
features = ["bicEn22h8","totalBic2h8"]
plot_labels = ["BicEnt2","TotalBic"]
dataset = pd.read_csv('characteristics_MedDbsOnOffFirst4_HAT.csv')

clas = dataset['DBS']
clas.where(clas == 0,'DBS On',inplace = True)
clas.where(clas == 'DBS On','DBS Off',inplace = True)
fig = plt.figure(figsize = (20,15),dpi = 350)
    
ax = fig.add_subplot(2, 1, 1)
sb.boxplot(x="DBS", y=features[0], data=dataset)
ax.set_ylabel(plot_labels[0],fontsize = 25)
ax.set_xlabel("DBS",fontsize = 25)
ax.tick_params(labelsize = 25)
    
ax = fig.add_subplot(2, 1, 2)
sb.boxplot(x="DBS", y=features[1], data=dataset)
ax.set_ylabel(plot_labels[1],fontsize = 25)
ax.set_xlabel("DBS",fontsize = 25)
ax.tick_params(labelsize = 25)

if save_plots:
    sel_folder = os.path.join(init_path,'plots boxplots')
    if not os.path.exists(sel_folder):
        os.mkdir(sel_folder)
    fname = os.path.join(sel_folder,'dbs_on_off_boxplots.png')
    fig.savefig(fname,bbox_inches='tight')
    
    
        