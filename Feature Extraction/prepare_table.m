%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment
%
% In this script the datasets for the statistical tests that presented in
% the tables II and III of (https://doi.org/10.36227/techrxiv.19589728.v1)
% are prepared

%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

clc
clear
close all

type = 'MedDbsOnOffFirst4';

ron = readtable('Features_ron.csv');
ref = readtable('Features_ref.csv');
ren = readtable('Features_ren.csv');
rof = readtable('Features_rof.csv');

ron(1,:) = [];%g1 out
ref(1,:) = [];%g1 out
ren(1,:) = [];%g1 out
rof([5,6],:) = []; %v4,v5 out 

new_class = [1,1,1,1,0,0,0,0,0,0,0,1];
new_class = repmat(new_class,1,4)';
new_class = array2table(new_class,'VariableNames',"new_class");

med = [ones(12,1); zeros(12,1); ones(12,1); zeros(12,1)];
med =  array2table(med,'VariableNames',"Med");

dbs = [zeros(12,1); ones(12,1); ones(12,1); zeros(12,1)];
dbs =  array2table(dbs,'VariableNames',"DBS");

names = rof.Properties.VariableNames;

all = [ron; ref; ren; rof]; 

all = [all,new_class,med,dbs];

allHAT = all(table2array(all(:,'class')) == 1,:);
allLAT = all(table2array(all(:,'class')) == 0,:);
writetable(allHAT,['characteristics_',type,'_HAT.csv']);  
writetable(allLAT,['characteristics_',type,'_LAT.csv']);  
writetable(all,['characteristics_',type,'.csv']);  
    





