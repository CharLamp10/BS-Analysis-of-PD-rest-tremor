"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the barplot with the selected features, as shown in table V
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
save = True
k = 8
all_vars = ['BicEnt1','BicEnt2','TotalBic','BspecEnt1','BspecEnt2',
            'MaxArea','TotalArea','Power46','ESR','Belongs46','NoPeaks']

""" Initializations """
names = ['Features_rof','Features_r15of','Features_r30of','Features_r45of','Features_r60of','total']
path_init = os.getcwd()
path_results1 = os.path.join(path_init, 'classification_results')
path_results2 = os.path.join(path_init, 'total_classification_results')
for j,name in enumerate(names):
    if name != 'total':
        path = os.path.join(path_results1, name)
    else:
        path = os.path.join(path_results2, name)
    
    path_selected = os.path.join(path, 'selected')

    dicti = np.load(os.path.join(path_selected, name + '_k' + str(k) +'.npy'),allow_pickle='TRUE').item()
    counts = np.array(dicti['kNN']['counts'])
    values = np.array(dicti['kNN']['values'])
    varnames = np.array(dicti['kNN']['varnames'])
    length = len(dicti['kNN']['selected'])

    counts = list(counts/(length)*100)
    varnames = varnames[values.astype(int)]


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
            varnames[i] = 'TotalBic'
        if varnames[i] == 'maxArea':
            varnames[i] = 'MaxArea'
        if varnames[i] == 'totalArea':
            varnames[i] = 'TotalArea'
        if varnames[i] == 'belongs46':
            varnames[i] = 'Belongs46'
        if varnames[i] == 'noPeaks':
            varnames[i] = 'NoPeaks'

    varnames = list(varnames)
    all_counts = np.zeros((11,))
    for i in range(len(varnames)):
        pos = all_vars.index(varnames[i])
        if pos != []:
            all_counts[pos,] = counts[i]
    all_counts = list(all_counts)
    all_counts = pd.DataFrame(all_counts)
    all_varss = pd.DataFrame(all_vars)
    if j == 0:
        df = pd.concat((all_varss,all_counts),axis = 1,)
    else:
        df = pd.concat((df,all_counts),axis = 1)
    
df.columns = ['Feature','rof','r15of','r30of','r45of','r60of','total']

if save:
    folder_path = os.path.join(path_init,'table selected')
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    fname = os.path.join(folder_path,'table_selected.csv')
    df.to_csv(fname,index = False)
