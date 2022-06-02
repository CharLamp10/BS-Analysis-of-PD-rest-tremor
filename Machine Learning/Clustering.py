"""
Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
and Classification Insights Pre-/Post-DBS and Medication Treatment

In this script the clustering that is described in the discussion section of
( https://doi.org/10.36227/techrxiv.19589728.v1) is implemented.

-----------------------------------------------------------------------------------------------------------------
Authors: Charalampos Lamprou & Ioannis Ziogas
Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
-----------------------------------------------------------------------------------------------------------------

"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn import manifold
from sklearn.metrics import silhouette_score
from collections import OrderedDict
from matplotlib.ticker import NullFormatter
from sklearn.neighbors import LocalOutlierFactor
from sklearn.ensemble import IsolationForest
from sklearn.cluster import MeanShift
from clusteval import clusteval
from itertools import combinations
import os

def outlier_detection(alg,data):
                    #IsolationForest(random_state,contamination,max_features)
                    #ISOForest.decision_function(data)
                    #LocalOutlierFactor(n_neighbors,contamination)
                    #LOF.negative_outlier_factor_
                    alg.fit(data)
                    preds = alg.predict(data)
                    outliers = np.where(preds == -1)[0]
                    cleaned_data = np.delete(data, (outliers), axis=0)
                    
                    return outliers,alg,cleaned_data

def mode_rows(a):
    a = np.ascontiguousarray(a)
    void_dt = np.dtype((np.void, a.dtype.itemsize * np.prod(a.shape[1:])))
    _,ids, count = np.unique(a.view(void_dt).ravel(), \
                                return_index=1,return_counts=1)
    largest_count_id = ids[count.argmax()]
    most_frequent_row = a[largest_count_id]
    return most_frequent_row,largest_count_id

def simple_matching_coeff(X,Y):
            X = np.array(X); Y = np.array(Y)
            m10 = 0; m01 = 0; m11 = 0; m00 = 0
            if not ((X==0) | (X==1)).all():
                raise Exception("Vector X is not binary, should contain only binary values")
            if not ((X==0) | (X==1)).all():
                raise Exception("Vector Y is not binary, should contain only binary values")
            if len(X) != len(Y):
                raise Exception("Vectors X,Y must have same length")
            for i in range(len(X)):
                if X[i] == Y[i] and X[i] == 1:
                    m11 += 1
                elif X[i] == Y[i] and X[i] == 0:
                    m00 += 1
                elif X[i] != Y[i] and X[i] == 1:
                    m10 += 1
                elif X[i] != Y[i] and X[i] == 0:
                    m01 += 1
            smc = (m11 + m00)/(m11 + m00 + m10 + m01)
            return smc    
  
def clustering_pipeline(dataset_np,clas,rand,n_components,max_clusters):
            allClasses = []        
            allScores = []
            "Scaling"    
            sc = StandardScaler()
            scData = sc.fit_transform(dataset_np)
                 
            "PCA"
            pca = PCA(n_components = n_components, svd_solver = 'full')
            pca.fit(scData)
            pcaData = pca.transform(scData)
            pcaData = pca.inverse_transform(pcaData)
        
            "Clustering"
            #Kmeans
            for m in range(rand):
                newClass = []
                scores = []
                for c in range(2,max_clusters):
                    kmeans = KMeans(n_clusters = c, init = 'k-means++', random_state = m+400)
                    cluster_labelsKmeans = kmeans.fit_predict(pcaData)
                    if c == 2 and cluster_labelsKmeans[1] != 1:
                        zeross = np.where(cluster_labelsKmeans == 0)[0]
                        oness = np.where(cluster_labelsKmeans == 1)[0]
                        cluster_labelsKmeans[zeross] = 1
                        cluster_labelsKmeans[oness]= 0
                        
                    silhouette_avgKmeans = silhouette_score(pcaData, cluster_labelsKmeans)
                    newClass.append(cluster_labelsKmeans)
                    penalty = 1/c
                    scores.append(penalty*silhouette_avgKmeans)
    
                max_score = max(scores)
                allScores.append(max_score)
                max_index = scores.index(max_score)
                allClasses.append(newClass[max_index])
            
            allClasses = np.array(allClasses)
            most_common_class,index = mode_rows(allClasses)
            score = allScores[index]*2
            smc = simple_matching_coeff(most_common_class,clas)
            
            return most_common_class,smc,score,scData

names = ['Features_rof','Features_r15of','Features_r30of','Features_r45of','Features_r60of']
rand = 100
max_clusters = 6
n_components = 1
save_new_dataset = True
path_init = os.getcwd()


remove = ['totalArea','belongs46']

    
for n in range(len(names)):
    dataset = pd.read_csv(names[n] + '.csv')
    clas = dataset.iloc[0:dataset.iloc[:,1].size,-1].values.reshape((-1,1)) 
        
    dataset_cut = dataset[dataset.columns.difference(remove,sort = False)]
    dataset_cut = dataset_cut.iloc[0:dataset_cut.iloc[:,1].size, 0:10-len(remove)]
    str1 = 'Bispec'
    dataset_np = dataset_cut.values    
    most_common_class,smc,silhouette,scData = clustering_pipeline(dataset_np,clas,rand,n_components,max_clusters)

    new_dataset = np.concatenate((dataset_np,most_common_class.reshape((-1,1))),axis = 1)
    col_names = list(dataset_cut.columns)
    col_names.append('Class')
    new_datasetDf = pd.DataFrame(data = new_dataset, columns = col_names)
    
    
    dictionary = {'data':new_datasetDf,
                  'old class':clas,
                  'silhouette':silhouette,
                  'SMC':smc}
    if save_new_dataset:
        df_folder = os.path.join(path_init,'clustering_results')
        if not os.path.exists(df_folder):
            os.mkdir(df_folder)
        df_path = os.path.join(df_folder,names[n] + ' ' + str1 + ' clustering dataset')            
        np.save(df_path+' dictionary.npy',dictionary)
        new_datasetDf.to_csv(df_path+' table.csv',index = False)

