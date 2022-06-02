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

    
        