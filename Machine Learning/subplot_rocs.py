"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the figure with the ROC curves, as in ( https://doi.org/10.36227/techrxiv.19589728.v1)
is created and saved

-----------------------------------------------------------------------------------------------------------------
Authors: Charalampos Lamprou & Ioannis Ziogas
Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
-----------------------------------------------------------------------------------------------------------------

"""

import os
import numpy as np
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.subplots import make_subplots

def plot_roc(fig,results,model,name,i,plot_labels,ksel,open_plots):
    tpr_orig = results['tpr_original']
    tpr_upper = results['tpr_upper']
    tpr_lower = results['tpr_lower']
    fpr_mean = results['fpr_mean']
    mean_auc = results['mean_auc']
    high_auc = results['high_auc']
    low_auc = results['low_auc']
    
    c_fill      = 'rgba(52, 152, 219, 0.2)'
    c_line      = 'rgba(52, 152, 219, 0.5)'
    c_line_main = 'rgba(41, 128, 185, 1.0)'
    c_grid      = 'rgba(189, 195, 199, 0.5)'
    i = i+1
    if i <= 3:
        row = 1
        col = i
    else:
        row = 2
        col = i - 3

    fig.add_trace(
        go.Scatter(        
        x          = fpr_mean,
        y          = tpr_upper,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'),row,col)
    fig.add_trace(
        go.Scatter(
            x          = fpr_mean,
            y          = tpr_lower,
            fill       = 'tonexty',
            fillcolor  = c_fill,
            line       = dict(color=c_line, width=1),
            hoverinfo  = "skip",
            showlegend = False,
            name       = 'lower'),row,col)
    fig.add_trace(
        go.Scatter(
            x          = fpr_mean,
            y          = tpr_orig,
            line       = dict(color=c_line_main, width=2),
            hoverinfo  = "skip",
            showlegend = False,
            name       = name + f' AUC: {mean_auc:.3f} [{low_auc:.3f} - {high_auc:.3f}]')
        ,row, col)
    """
    fig.add_shape(
        type ='line', 
        line =dict(dash='dash'),
        x0=0, x1=1, y0=0, y1=1
    )
    """
    
    fig.update_layout(
        #title_text = model,
        #title_font_size= 10,
        template    = 'plotly_white', 
        #title_x     = 0.5,
        xaxis1_title = "1-Specificity",
        yaxis1_title = "Sensitivity",
        #xaxis2_title = "1-Specificity",
        #yaxis2_title = "Sensitivity",
        #xaxis3_title = "1-Specificity",
        #yaxis3_title = "Sensitivity",
        xaxis4_title = "1-Specificity",
        yaxis4_title = "Sensitivity",
        #xaxis5_title = "1-Specificity",
        #yaxis5_title = "Sensitivity",
        #xaxis6_title = "1-Specificity",
        #yaxis6_title = "Sensitivity",
        font=dict(
        size=15),
        width       = 800,
        height      = 800,
        
        legend      = dict(
            #yanchor="bottom", 
            #xanchor="right", 
            x=0.95,
            y=0.01,
        ),
        
        yaxis1 = dict(range=[0, 1]),
        yaxis2 = dict(range=[0, 1]),
        yaxis3 = dict(range=[0, 1]),
        yaxis4 = dict(range=[0, 1]),
        yaxis5 = dict(range=[0, 1]),
        yaxis6 = dict(range=[0, 1])
    )
    """
    fig.update_yaxes(
        range       = [0, 1],
        gridcolor   = c_grid,
        scaleanchor = "x", 
        scaleratio  = 1,
        linecolor   = 'black')
    fig.update_xaxes(
        range       = [0, 1],
        gridcolor   = c_grid,
        constrain   = 'domain',
        linecolor   = 'black')
    """
    
    return fig



path_init = os.getcwd()#PTES_Python/Features
method = 'classification' #classification,total_classification
path_method1 = os.path.join(path_init,method + '_results')
path_method2 = os.path.join(path_init,'total_classification_results')

k = 8
save_plots = True
open_plots = False

names = ['Features_rof','Features_r15of','Features_r30of','Features_r45of','Features_r60of','total']
model_list = ["LinearSVC"]
plot_labels = ''
fields = ['tpr_original','tpr_upper','tpr_lower','mean_auc','high_auc',
          'low_auc','c_fill','c_line','c_line_main']

for model in model_list:
    all_results = {}
    all_results = {f:[] for f in fields}
    for i,name in enumerate(names):
        if name == 'total':
            fname = os.path.join(path_method2,model)
        else:
            fname = os.path.join(path_method1,model)
        if not os.path.exists(fname):               
            os.mkdir(fname)
        folderpath = os.path.join(path_init,'ROC subplots')
        if not os.path.exists(folderpath):               
            os.mkdir(folderpath)
       
        if name == 'total':
            filepath = os.path.join(path_method2,name,'boot_rocs',name+'_'+model+'_k'+str(k)+'.npy')
        else:
            filepath = os.path.join(path_method1,name,'boot_rocs',name+'_'+model+'_k'+str(k)+'.npy')
        btstrp_results =  np.load(filepath,allow_pickle='TRUE').item()
        ksel = btstrp_results['k']
        all_results['tpr_original'] = btstrp_results['tpr_original']
        all_results['tpr_upper'] = btstrp_results['tpr_upper']
        all_results['tpr_lower'] = btstrp_results['tpr_lower']
        all_results['mean_auc'] = btstrp_results['mean_auc']
        all_results['high_auc'] = btstrp_results['high_auc']
        all_results['low_auc'] = btstrp_results['low_auc']
        all_results['fpr_mean'] = btstrp_results['fpr_mean']
        if i == 0:
            fig = make_subplots(rows=2, cols=3,subplot_titles=("rof", "r15of", "r30of", "r45of","r60of","Grand Average"))
            fig.update_annotations(font_size=18)
        if name == 'total':
            name = 'Grand Average'
        else:
            name = name[9:]
        fig = plot_roc(fig,all_results,model,name,i,plot_labels,ksel,open_plots)
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=1,
    col=1
    )
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=1,
    col=2
    )
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=1,
    col=3
    )
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=2,
    col=1
    )
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=2,
    col=2
    )
    fig.add_shape(
    dict(type="line", x0=0, x1=1, y0=0, y1=1,line = dict(dash = 'dash')),
    row=2,
    col=3
    )
    plot(fig, auto_open = open_plots)
    if save_plots:
        img_path = os.path.join(folderpath, model + ' k' + str(ksel))
        fig.write_image(img_path + ".png", engine = "orca")