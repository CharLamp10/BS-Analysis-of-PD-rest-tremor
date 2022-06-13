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

import pandas as pd
import os
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon

save = True
path_init = os.getcwd()
case = 'off' # off, first4

if case == 'off':
    names = ['Features_rof','Features_r15of','Features_r30of','Features_r45of',
         'Features_r60of','total_dataset']

    conds = ['Features','rof','r15of','r30of','r45of','r60of','total']

    for name in names:
        path_data = os.path.join(path_init, name + '.csv')
        data = pd.read_csv(path_data)
        data.rename(columns={"bspecEn12h8": "BspecEnt1", "power46": "Power46",
                             "bspecEn22h8": "BspecEnt2", "bicEn12h8": "BicEnt1",
                             "bicEn22h8": "BicEnt2", "ellipse_similarity_ratio": "ESR",
                             "totalBic2h8": "TotalBic", "maxArea": "MaxArea",
                             "totalArea": "TotalArea"})
        data.drop(columns = ['belongs46','noPeaks'], inplace = True)
        data1 = data[data['class'] == 1]
        data0 = data[data['class'] == 0]
        pval_ttest = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
        features = pd.DataFrame(data.iloc[:,0:9].columns)
        if name == 'Features_rof':
            ttest_df_Off = pd.concat((features,pd.DataFrame(pval_ttest[1])), axis = 1)
        else:
            ttest_df_Off = pd.concat((ttest_df_Off,pd.DataFrame(pval_ttest[1])), axis = 1)
        
        U = []
        Wx = []
        for i in range(len(features)):
            [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
            U.append(pvalU)
            """
            [_,pvalWx] = wilcoxon(x = data1.iloc[:,i].values, y = data0.iloc[:,i].values,
                                        zero_method='zsplit')
            Wx.append(pvalWx)
            """
            
        if name == 'Features_rof':
            U_df_Off = pd.concat((features,pd.DataFrame(U)), axis = 1)
            #Wx_df = pd.concat((features,pd.DataFrame(Wx)), axis = 1)
        else:
            U_df_Off = pd.concat((U_df_Off,pd.DataFrame(U)), axis = 1) 
            #Wx_df = pd.concat((Wx_df,pd.DataFrame(Wx)), axis = 1) 
            
        vwv = 1
            
    ttest_df_Off.columns = conds
    U_df_Off.columns = conds
    if save:
        save_path = os.path.join(path_init,'p_values')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        fname_ttest = os.path.join(save_path,'Off_ttest_pval.csv')
        fname_Utest = os.path.join(save_path,'Off_Utest_pval.csv')
        ttest_df_Off.to_csv(fname_ttest)
        U_df_Off.to_csv(fname_Utest)
    
else:
    names = ['characteristics_MedDbsOnOffFirst4','characteristics_MedDbsOnOffFirst4_HAT',
             'characteristics_MedDbsOnOffFirst4_LAT']
    
    for j,name in enumerate(names):
        path_data = os.path.join(path_init, name + '.csv')
        data = pd.read_csv(path_data)
        data.rename(columns={"bspecEn12h8": "BspecEnt1", "power46": "Power46",
                             "bspecEn22h8": "BspecEnt2", "bicEn12h8": "BicEnt1",
                             "bicEn22h8": "BicEnt2", "ellipse_similarity_ratio": "ESR",
                             "totalBic2h8": "TotalBic", "maxArea": "MaxArea",
                             "totalArea": "TotalArea"})
        data.drop(columns = ['belongs46','noPeaks'], inplace = True)
        if j == 0:
            conds = ['Features','HAT/LAT','Med','DBS']
            features = pd.DataFrame(data.iloc[:,0:9].columns)
            data1 = data[data['class'] == 1]
            data0 = data[data['class'] == 0]
            pval_ttestClass = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            UClass = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UClass.append(pvalU)
                
            data1 = data[data['Med'] == 1]
            data0 = data[data['Med'] == 0]
            pval_ttestMed = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            UMed = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UMed.append(pvalU)
                
            data1 = data[data['DBS'] == 1]
            data0 = data[data['DBS'] == 0]
            pval_ttestDBS = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            ttest_First4 = pd.concat((features,pd.DataFrame(pval_ttestClass[1]),
                                    pd.DataFrame(pval_ttestMed[1]),
                                    pd.DataFrame(pval_ttestDBS[1])),axis = 1)  
            UDBS = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UDBS.append(pvalU)
                
            U_First4 = pd.concat((features,pd.DataFrame(UClass),pd.DataFrame(UMed),
                                  pd.DataFrame(UDBS)), axis = 1)
             
            ttest_First4.columns = conds
            U_First4.columns = conds
        elif j == 1:
            conds = ['Features','Med','DBS']
            features = pd.DataFrame(data.iloc[:,0:9].columns)               
            data1 = data[data['Med'] == 1]
            data0 = data[data['Med'] == 0]
            pval_ttestMed = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            UMed = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UMed.append(pvalU)
                
            data1 = data[data['DBS'] == 1]
            data0 = data[data['DBS'] == 0]
            pval_ttestDBS = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            ttest_First4HAT = pd.concat((features,
                                    pd.DataFrame(pval_ttestMed[1]),
                                    pd.DataFrame(pval_ttestDBS[1])),axis = 1)  
            UDBS = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UDBS.append(pvalU)
                
            U_First4HAT = pd.concat((features,pd.DataFrame(UMed),
                                  pd.DataFrame(UDBS)), axis = 1)
             
            ttest_First4HAT.columns = conds
            U_First4HAT.columns = conds
        elif j == 2:
            conds = ['Features','Med','DBS']
            features = pd.DataFrame(data.iloc[:,0:9].columns)               
            data1 = data[data['Med'] == 1]
            data0 = data[data['Med'] == 0]
            pval_ttestMed = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            UMed = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UMed.append(pvalU)
                
            data1 = data[data['DBS'] == 1]
            data0 = data[data['DBS'] == 0]
            pval_ttestDBS = ttest_ind(data1.iloc[:,0:9].values, data0.iloc[:,0:9].values, equal_var=True)
            ttest_First4LAT = pd.concat((features,
                                    pd.DataFrame(pval_ttestMed[1]),
                                    pd.DataFrame(pval_ttestDBS[1])),axis = 1)  
            UDBS = []
            for i in range(len(features)):
                [_,pvalU] = mannwhitneyu(data1.iloc[:,i].values,data0.iloc[:,i].values)
                UDBS.append(pvalU)
                
            U_First4LAT = pd.concat((features,pd.DataFrame(UMed),
                                  pd.DataFrame(UDBS)), axis = 1)
             
            ttest_First4LAT.columns = conds
            U_First4LAT.columns = conds
    if save:
        save_path = os.path.join(path_init,'p_values')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        fname_ttest_first4 = os.path.join(save_path,'First4_ttest_pval.csv')
        fname_Utest_first4 = os.path.join(save_path,'First4_Utest_pval.csv')
        fname_ttest_first4HAT = os.path.join(save_path,'First4HAT_ttest_pval.csv')
        fname_Utest_first4HAT = os.path.join(save_path,'First4HAT_Utest_pval.csv')
        fname_ttest_first4LAT = os.path.join(save_path,'First4LAT_ttest_pval.csv')
        fname_Utest_first4LAT = os.path.join(save_path,'First4LAT_Utest_pval.csv')
        ttest_First4.to_csv(fname_ttest_first4)
        U_First4.to_csv(fname_Utest_first4)
        ttest_First4HAT.to_csv(fname_ttest_first4HAT)
        U_First4HAT.to_csv(fname_Utest_first4HAT)
        ttest_First4LAT.to_csv(fname_ttest_first4LAT)
        U_First4LAT.to_csv(fname_Utest_first4LAT)
    
    

        
        
        
        
        
        
        
    

