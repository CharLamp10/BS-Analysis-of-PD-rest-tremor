"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the barplot with the selected features, as shown in figure 3
of ( https://doi.org/10.36227/techrxiv.19589728.v1).

-----------------------------------------------------------------------------------------------------------------
Authors: Charalampos Lamprou & Ioannis Ziogas
Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
-----------------------------------------------------------------------------------------------------------------

"""

import numpy as np
import pandas as pd
import numpy.matlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import seaborn as sb

sb.set_theme()
save_plots = True

""" Initializations """
path_init = os.getcwd()
path_results = os.path.join(path_init, 'classification_results')
path_chars = os.path.join(path_results, 'characteristics_MedDbsOnOffFirst4')
path_selected = os.path.join(path_chars, 'selected')


df = pd.DataFrame() 
dict_DBS = np.load(os.path.join(path_selected, 'DBS_k3.npy'),allow_pickle='TRUE').item()
dict_Med = np.load(os.path.join(path_selected, 'Med_k4.npy'),allow_pickle='TRUE').item()
dict_class = np.load(os.path.join(path_selected, 'class_k4.npy'),allow_pickle='TRUE').item()

selected_DBS = np.array(dict_DBS['k3']['selected'])
selected_Med = np.array(dict_Med['k4']['selected'])
selected_class = np.array(dict_class['k4']['selected'])

values_DBS, counts_DBS = np.unique(selected_DBS, return_counts=True)
values_Med, counts_Med = np.unique(selected_Med, return_counts=True)
values_class, counts_class = np.unique(selected_class, return_counts=True)

counts_DBS = counts_DBS/(sum(counts_DBS)/3)*100
counts_Med = counts_Med/(sum(counts_Med)/4)*100
counts_class = counts_class/(sum(counts_class)/4)*100

values_DBS = values_DBS[np.where(counts_DBS > 1)]
counts_DBS = pd.DataFrame(counts_DBS[np.where(counts_DBS > 1)])
values_Med = values_Med[np.where(counts_Med > 1)]
counts_Med = pd.DataFrame(counts_Med[np.where(counts_Med > 1)])
values_class = values_class[np.where(counts_class > 1)]
counts_class = pd.DataFrame(counts_class[np.where(counts_class > 1)])

counts_DBS = pd.DataFrame(counts_DBS)
counts_Med = pd.DataFrame(counts_Med)
counts_class = pd.DataFrame(counts_class)

varnames_DBS = np.array(dict_DBS['k3']['varnames'])
varnames_Med = np.array(dict_Med['k4']['varnames'])
varnames_class = np.array(dict_class['k4']['varnames'])

varnames_DBS = varnames_DBS[values_DBS.astype(int)]
varnames_Med = varnames_Med[values_Med.astype(int)]
varnames_class = varnames_class[values_class.astype(int)]

varnames = list(varnames_Med) + list(varnames_class) + list(varnames_DBS)

for i in range(len(varnames)):
    if varnames[i] == 'power46':
        varnames[i] = 'Power46'
    if varnames[i] == 'bspecEn12h8':
        varnames[i] = 'BspecEnt1'
    if varnames[i] == 'bspecEn22h8':
        varnames[i] = 'BspecEnt2'
    if varnames[i] == 'bicEn12h8':
        varnames[i] = 'BicEnt1'
    if varnames[i] == 'bicEn22h8':
        varnames[i] = 'BicEnt2'
    if varnames[i] == 'ellipse_similarity_ratio':
        varnames[i] = 'ESR'
    if varnames[i] == 'totalBic2h8':
        varnames[i] = 'totalBic'

varnames = pd.DataFrame(varnames)

hue = []
for i in range(len(values_Med)):
    hue.append('Med On/Off')
for i in range(len(values_class)):
    hue.append('HAT/LAT')
for i in range(len(values_DBS)):
    hue.append('DBS On/Off')

hue = pd.DataFrame(hue)

counts = pd.concat((counts_Med,counts_class,counts_DBS), axis = 0, ignore_index = 'FALSE')

df = pd.concat((varnames,counts,hue),axis = 1,ignore_index = True)

df = df.rename(columns = {0:'Features', 1:'Counts', 2:''})
plt.figure(dpi = 600)

plot = sb.barplot(x="Features", y = "Counts", hue = "", data = df, palette="Blues_d")    
plt.xticks(rotation=70)
plt.tick_params(labelsize = 7)
plt.ylabel('Percentage (%)',fontsize = 10)
plt.xlabel('',fontsize = 10)
plt.legend(loc='upper left', fontsize='7')

if save_plots:
    sel_folder = os.path.join(path_chars,'plots selected')
    if not os.path.exists(sel_folder):
        os.mkdir(sel_folder)
    fname = os.path.join(sel_folder,'plot_selected.png')
    plt.savefig(fname,bbox_inches='tight')
#plot.set_xticklabels(plot.get_xticks(), size=5)