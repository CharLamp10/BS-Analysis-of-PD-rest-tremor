"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the nested LOSO classification scheme for the DBS, Med, HAT/LAT is implemented, as described
in ( https://doi.org/10.36227/techrxiv.19589728.v1).

-----------------------------------------------------------------------------------------------------------------
Authors: Charalampos Lamprou & Ioannis Ziogas
Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
-----------------------------------------------------------------------------------------------------------------

"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import f1_score,precision_score,recall_score,accuracy_score, roc_curve
from sklearn.model_selection import GridSearchCV, LeaveOneOut
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_auc_score,roc_curve
import scipy.stats
import os
from imblearn.over_sampling import SMOTE
import plotly.graph_objects as go
from plotly.offline import plot
from sklearn.pipeline import make_pipeline
import seaborn as sb
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a,axis = 0), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h

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
    
    fig = go.Figure([
        go.Scatter(        
        x          = fpr_mean,
        y          = tpr_upper,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'),
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower,
        fill       = 'tonexty',
        fillcolor  = c_fill,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_orig,
        line       = dict(color=c_line_main, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'AUC: {mean_auc:.3f} [{low_auc:.3f} - {high_auc:.3f}]')
    ])
    fig.add_shape(
        type ='line', 
        line =dict(dash='dash'),
        x0=0, x1=1, y0=0, y1=1
    )
    fig.update_layout(
        title_text = model + ' ' + name,
        title_font_size= 30,
        template    = 'plotly_white', 
        title_x     = 0.5,
        xaxis_title = "1 - Specificity",
        yaxis_title = "Sensitivity",
        width       = 800,
        height      = 800,
        legend      = dict(
            yanchor="bottom", 
            xanchor="right", 
            x=0.95,
            y=0.01,
        )
    )
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
    plot(fig, auto_open = False)
    
    return fig

def load_and_oversample(dataset,clas,colstart,colend,ovsmp,random_seed = 42,samp_ratio = 0.75):
    
        #clas = dataset.iloc[:,-1].values.reshape((-1,1))
        clas = dataset[clas].values.reshape((-1,1))
        class1 = np.where(clas == 1)[1]
        class0 = np.where(clas == 0)[1]
        no1 = len(class1)
        no0 = len(class0)
        dataset = dataset.iloc[:,colstart:colend]
        varnames = list(dataset.columns)
        dataset = dataset.values
        if ovsmp and ((no1/no0) < samp_ratio or (no0/no1) < samp_ratio) :
            ros = SMOTE(random_state = random_seed,sampling_strategy = samp_ratio,k_neighbors = no1 - 1)
            dataset, clas = ros.fit_resample(dataset, clas)
            clas = clas.reshape((-1,1))
        dataset = np.concatenate((dataset,clas), axis = 1)
        return dataset,clas,varnames
      
def construct_roc(score,y,step):
    # false positive rate
    fpr = []
    # true positive rate
    tpr = []
    # Iterate thresholds from 0.0, 0.01, ... 1.0
    thresholds = np.arange(0.0, 1+step, step)
    
    # get number of positive and negative examples in the dataset
    P = sum(y)
    N = len(y) - P
    
    # iterate through all thresholds and determine fraction of true positives
    # and false positives found at this threshold
    for thresh in thresholds:
        FP=0
        TP=0
        for i in range(len(score)):
            if (score[i] > thresh):
                if y[i] == 1:
                    TP = TP + 1
                if y[i] == 0:
                    FP = FP + 1
        fpr.append(FP/float(N))
        tpr.append(TP/float(P))
    return fpr,tpr

def bootstrap_roc(results,random_seed,model,name,B = 1000,threshold_length = 100):
    
    y_pred_proba = np.mean(np.squeeze(np.array(results['proba'])),0).tolist()
    y_true = np.array(results['class']).transpose()
    score_orig = roc_auc_score(y_true, y_pred_proba)
    fpr_orig,tpr_orig,_ = roc_curve(y_true,y_pred_proba,pos_label = 1)
    
    boot_fpr = []
    boot_tpr =[]
    boot_auc = []
    
    # control reproducibility    
    rng = np.random.RandomState(random_seed)
    for i in range(B):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.randint(0, len(y_pred_proba), len(y_pred_proba))
        if len(np.unique(y_true[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        score = roc_auc_score(y_true[indices], np.array(y_pred_proba)[indices])
        fpr,tpr,_ = roc_curve(y_true[indices],np.array(y_pred_proba)[indices],pos_label = 1)
        boot_auc.append(score)
        boot_fpr.append(fpr)
        boot_tpr.append(tpr)
    
    fpr_mean    = np.linspace(0, 1, threshold_length)
    tpr_orig    = np.interp(fpr_mean, fpr_orig, tpr_orig)
    tpr_orig[0] = 0.0
    tpr_orig[-1] = 1.0
    interp_tprs = []
    for i in range(len(boot_fpr)):
        fpr           = np.array(boot_fpr)[i]
        tpr           = np.array(boot_tpr)[i]
        interp_tpr    = np.interp(fpr_mean, fpr, tpr)
        interp_tpr[0] = 0.0
        interp_tprs.append(interp_tpr)
    
    mean_auc = np.mean(boot_auc,axis = 0)
    std_auc = 2*np.std(boot_auc,axis = 0)
    high_auc = np.clip(mean_auc + std_auc,0,1)
    low_auc = mean_auc - std_auc
    tpr_mean     = np.mean(interp_tprs, axis=0)
    tpr_mean[-1] = 1.0
    tpr_std      = 2*np.std(interp_tprs, axis=0)
    tpr_upper    = np.clip(tpr_mean+tpr_std, 0, 1)
    tpr_lower    = tpr_mean-tpr_std
    
    btstrp_results = {'tpr_original':tpr_orig,
               'tpr_upper':tpr_upper,
               'tpr_lower':tpr_lower,
               'tpr_mean':tpr_mean,
               'tpr_std':tpr_std,
               'tpr_orig':tpr_orig,
               'fpr_mean':fpr_mean,
               'mean_auc':mean_auc,
               'high_auc':high_auc,
               'low_auc':low_auc}
    return btstrp_results  
 
def plot_single_roc(results,model,name):

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
    
    fig = go.Figure([
        go.Scatter(        
        x          = fpr_mean,
        y          = tpr_upper,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'),
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower,
        fill       = 'tonexty',
        fillcolor  = c_fill,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_orig,
        line       = dict(color=c_line_main, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'AUC: {mean_auc:.3f} [{low_auc:.3f} - {high_auc:.3f}]')
    ])
    fig.add_shape(
        type ='line', 
        line =dict(dash='dash'),
        x0=0, x1=1, y0=0, y1=1
    )
    fig.update_layout(
        title_text = model + ' ' + name,
        title_font_size= 30,
        template    = 'plotly_white', 
        title_x     = 0.5,
        xaxis_title = "1 - Specificity",
        yaxis_title = "Sensitivity",
        width       = 800,
        height      = 800,
        legend      = dict(
            yanchor="bottom", 
            xanchor="right", 
            x=0.95,
            y=0.01,
        )
    )
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
    plot(fig, auto_open = False)
    
    return fig                             
"Initializations"

models = [SVC()]
kscore_func = f_classif

path_init = os.getcwd()
path_init = os.path.join(path_init, 'classification_results')

names = ['characteristics_MedDbsOnOffFirst4']
clas = 'class'#Med, DBS
rand = 100
model_list = ["rbfSVC"]
fields = ['auc','acc','sens','spec','proba','class']
sel_fields = ['selected','selected2','varnames','values','counts']

oversampling = 1
save_plots = 0
save_selected = 1
save_results = 1
save_boot_results = 1
n_jobs = -1

for ksel in range (4,5):
    for name in names:
    
        results = {}  
        selected = []
        
        for model in model_list:
            results[model] = {}
        
        predictions = {}
        
        for model in model_list:
            predictions[model] = {}
            
            fname = os.path.join(path_init,name)
            if not os.path.exists(fname):
                os.mkdir(fname)
            
            results = {}
            for model in model_list:
                results[model] = {f:[] for f in fields}
            results['k'] = ksel
            
        
        
        """Iterate through random states - 1.Different oversampling each time,
        2.Different grid initialization"""
        for m in range(rand):
            dataset = pd.read_csv(name + '.csv')
            colstart = 0
            colend = 11
            kselec = ksel
            dataset,_,varnames = load_and_oversample(dataset,clas,colstart,colend,oversampling,random_seed = m,samp_ratio = 0.5)
            cols = len(dataset[1,:])
            rows = len(dataset[:,1])              
                
            X = dataset[:,0:-1]
            y = dataset[:,-1]
            
 
            """Initialize Grids"""
        
            
            param_grid_linear_svm = {'C': [0.001,0.01,0.1,1,10,100],
            'gamma': [0.001,0.01,0.1,1,10,100],
            'kernel': ['rbf'],
            'probability':[True],
            'random_state':[m]}
            
            
            param_grid = [param_grid_linear_svm]
 

 
            """Leave-One-Out Pipeline""" 

            for j,model in enumerate(models):
                cv = LeaveOneOut()
                grid = GridSearchCV(model, param_grid[j], cv = cv,n_jobs = n_jobs)
                clf = make_pipeline(StandardScaler(), SelectKBest(f_classif, k = kselec), grid)
                
                clf.fit(X, y)
                y_pred = clf.predict(X)
                conf_mat = confusion_matrix(y,y_pred)
                results[model_list[j]]['sens'].append(conf_mat[0,0]/(conf_mat[0,0]+conf_mat[0,1]))
                results[model_list[j]]['spec'].append(conf_mat[1,1]/(conf_mat[1,1]+conf_mat[1,0]))
                results[model_list[j]]['acc'].append(clf.score(X,y))
                probs = clf.predict_proba(X)[:,1]
                results[model_list[j]]['proba'].append(probs)
                results[model_list[j]]['auc'].append(roc_auc_score(y,probs))
                
                if j == 0 and m == 0:
                    for train_index, test_index in cv.split(X):
                        selector = SelectKBest(f_classif,k = kselec)
                        selector.fit(X[train_index,:],y[train_index])
                        selected.append(selector.get_support(indices = True))
                
       
        
        for j,model_name in enumerate(model_list):
            results[model_list[j]]['class'].append(y)  
            
        "Save mean,+-CI of each feature in results dict"         
        metrics = ['auc','acc','sens','spec']
        for model in model_list:
            for mi in metrics:
                CI = []
                meanf,meanf_low,meanf_high = mean_confidence_interval(
                    results[model][mi], confidence=0.95)
                if meanf_high > 1:
                    meanf_high = 1
                CI.append(meanf_low)
                CI.append(meanf)
                CI.append(meanf_high)
                results[model][mi] = CI
                
       
        "Prepare folder and paths for saving plots"
        
        folderpath = os.path.join(fname,'plots')
        if not os.path.exists(folderpath):               
            os.mkdir(folderpath)
            
        "Bootstrap and plot ROC"
        for model in model_list:
            btstrp_results = bootstrap_roc(results[model],42,model,name,B = 1000,threshold_length = 100)
            btstrp_results['k'] = ksel
            fig = plot_single_roc(btstrp_results,model,name)
            if save_boot_results:
                boot_folder = os.path.join(fname,'boot_rocs')
                if not os.path.exists(boot_folder):
                    os.mkdir(boot_folder)
                boot_file = os.path.join(boot_folder,name + '_' + model + '_' + 'k' + str(ksel) + '_' + clas +'.npy')
                np.save(boot_file,btstrp_results)
            if save_plots:
                fig.write_image(os.path.join(folderpath,model+clas + ".png"))
                
        "Save Results"                        
        folderpath = os.path.join(fname,'results')
        if save_results:
            if not os.path.exists(folderpath):
                os.mkdir(folderpath)
            savepath = os.path.join(folderpath,name + '_' + 'k' + str(ksel) + '_' + clas + ' rbf sb.npy')
            np.save(savepath,results)       

            
        
                
                
