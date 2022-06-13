# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 10:55:39 2021

@author: HP
"""

import numpy as np
import pandas as pd
import numpy.matlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import seaborn as sb

sb.set_theme()

""" Initializations """
case = 'Bispec'#'Combined', 'Bispec'
cat = '_new'#'_new'
path_init = os.getcwd() #Ptes_Python/Features
method = 'total_classification'#'total_classification'
path_method = os.path.join(path_init,method + '_results')
conds = ["rof","total"] #['total'], ["rof","r15of","r30of","r45of","r60of"]
condsflip = np.flip(conds)
#k = 6
model_list = ["SVC","LinearSVC","kNN","LogRe","DecisionTree"]#,"RFC"]
save_plots = False

varnamesinit = np.array(["MaxArea","TotalArea","ESR","NoPeaks","Belongs46",
            "Power46","BspecEnt1","BspecEnt2","TotalBic","BicEnt1","BicEnt2"])

for k in range(8,9):
    df = pd.DataFrame()
    for cond in conds:    
        if method == 'total_classification':
            path_condition = os.path.join(path_method,cond + cat)
            name = cond + '_k' + str(k) + '.npy'
        else:
            path_condition = os.path.join(path_method, cond + cat)
            name = 'Features_' + cond + '_k' + str(k) + '.npy'
            
        path_selected = os.path.join(path_condition,'selected '+ case)
        
        dictionary = np.load(os.path.join(path_selected,name),allow_pickle='TRUE').item()
        selected = np.array(dictionary[model_list[0]]['selected'])
        values = dictionary[model_list[0]]['values'].reshape((-1,1))
        counts = dictionary[model_list[0]]['counts']
        #varnames = np.array(dictionary[model_list[0]]['varnames'])
        varnames = varnamesinit[values] 
        varnames_new = []
        for j in range(len(varnames)):
            temp = np.matlib.repmat(varnames[j,0],counts[j],1)
            
            if len(varnames_new) == 0:
                varnames_new = temp
            else:
                varnames_new = np.concatenate((varnames_new,temp),axis = 0)
        
        mod = np.matlib.repmat(cond,len(varnames_new),1)
        temp_df = pd.DataFrame(np.concatenate((varnames_new,mod),axis = 1))
        df = pd.concat((df,temp_df),axis = 0,ignore_index = True)
    
    
    df = df.rename(columns = {0:'Features',
                                  1:'Conditions'})
    plt.figure(dpi = 350)
    
    plot = sb.histplot(data = df,stat = "percent",common_norm = False,
                        x="Features", hue="Conditions", multiple="dodge", shrink=.8)    
    plt.xticks(rotation=45)
    plt.setp(plot.get_legend().get_texts(), fontsize='8')
    plt.setp(plot.get_legend().get_title(), fontsize='10')
    plt.title('Best Selected Features')
    #plt.title(case +' Case: k = ' + str(k) + ' Best Selected Features')
    plt.tick_params(labelsize = 7)
    plt.xlabel('Features',fontsize=10)
    plt.ylabel('Percentage (%)',fontsize = 10)
    plt.xlim([-1,12])
    if save_plots:
        sel_folder = os.path.join(path_method,'plots selected ' + case)
        if not os.path.exists(sel_folder):
            os.mkdir(sel_folder)
        fname = os.path.join(sel_folder,'selected_' + case + '_k' + str(k) +'.png')
        plt.savefig(fname)
    #plot.set_xticklabels(plot.get_xticks(), size=5)

        
    