%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment
%
% In this script the plots that illustrate the dynamic nature of the tremor as
% presented through the extracted features (as shown in figure 5 of 
% ( https://doi.org/10.36227/techrxiv.19589728.v1)), are extracted and saved.

%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

clc
clear
close all


figstate = 'off'; % 'on' to open figures and 'off' otherwise
savep = 0; % 1 to save figiures and 0 otherwise
alpha = 0.05; % The level of statistical significance
vars = {'maxArea','totalArea','ellipse_similarity_ratio','noPeaks',...
            'belongs46','power46','bspecEn12h8','bspecEn22h8','totalBic2h8',...
            'bicEn12h8','bicEn22h8','mean_FFT_X','max_FFT_X','std_FFT_X'}; % The wanted features
time_labels = {'mOn dOff','mOff dOn','mOff dOff','mOff dOff 15',...
    'mOff dOff 30','mOff dOff 45','mOff dOff 60','mOn dOn'};


for v = 1:length(vars)
    selFeat = vars{1,v};
    
    dir1Name = ['time_plots_'];
    if ~exist(dir1Name,'dir') 
        mkdir(dir1Name); 
        dir2Name = [dir1Name,'\',selFeat];
        if ~exist(dir2Name,'dir')
            mkdir(dir2Name);
        end
    else
        dir2Name = [dir1Name,'\',selFeat];
        if ~exist(dir2Name,'dir')
            mkdir(dir2Name);
        end
    end

    
    ron = readtable('Features_ron_.csv');
    ref = readtable('Features_ref_.csv');
    rof = readtable('Features_rof_.csv');
    r15of = readtable('Features_r15of_.csv');
    r30of = readtable('Features_r30of_.csv');
    r45of = readtable('Features_r45of_.csv');
    r60of = readtable('Features_r60of_.csv');
    ren = readtable('Features_ren_.csv');

    
    ron([1,4],:) = [];%g1,s7 out
    ref([1,4],:) = [];%g1,s7 out
    ren([1,4],:) = [];%g1,s7 out
    rof([3,5,6],:) = []; %s7,v4,v5 out 
    r30of(3,:) = []; %s7 out
    r60of(3,:) = []; %s7 out

    sf = [];
    sf = [sf ;table2array(ron(:,selFeat))];
    sf = [sf ;table2array(ref(:,selFeat))];
    sf = [sf ;table2array(rof(:,selFeat))];
    sf = [sf ;table2array(r15of(:,selFeat))];
    sf = [sf ;table2array(r30of(:,selFeat))];
    sf = [sf ;table2array(r45of(:,selFeat))];
    sf = [sf ;table2array(r60of(:,selFeat))];
    sf = [sf ;table2array(r60of(:,selFeat))];
    
    maxF = max(sf);
    minF = min(sf);
    
    names = rof.Properties.VariableNames;
    rowNames = {'g2','s6','s8','g10','g11','g12','g13','g9','s14','s15','s16'};
    ron.Properties.RowNames = rowNames; ref.Properties.RowNames = rowNames;
    rof.Properties.RowNames = rowNames; r15of.Properties.RowNames = rowNames;
    r30of.Properties.RowNames = rowNames; r45of.Properties.RowNames = rowNames;
    r60of.Properties.RowNames = rowNames; ren.Properties.RowNames = rowNames;

    condStruct.ron = ron; condStruct.ref = ref; condStruct.rof = rof;
    condStruct.r15of = r15of; condStruct.r30of = r30of; condStruct.r45of = r45of;
    condStruct.r60of = r60of; condStruct.ren = ren;

    fnames = fieldnames(condStruct);

    figure('Visible',figstate);
    LnLList = [1,2,3,11];
    msg = ['All LnL subjects, time evolution of feature ',selFeat];
    tetempFeat = zeros(1,8);
    for f = 1:length(fnames)
        tempFeat = zeros(1,4);
        for h = 1:length(LnLList)
            tempFeat(h) = table2array(condStruct.(fnames{f,1})(rowNames{1,LnLList(h)},selFeat));
        end    
        ci = bootci(100,@mean,tempFeat);
        cil(f) = ci(1);
        ciu(f) = ci(2);
        tetempFeat(f) = mean(tempFeat);
    end
    plot(tetempFeat,'LineStyle','--','LineWidth',1.5,'Marker','o')
    hold on 
    plot(cil,'LineStyle','--','LineWidth',1.5,'Color','r')
    legend('mean','95 % CI','AutoUpdate','off')
    hold on
    plot(ciu,'LineStyle','--','LineWidth',1.5,'Color','r')
    ylim([minF,maxF])
    title(msg)
    ylabel(selFeat)
    xlabel('Conditions Time Evolution')
    xticklabels(time_labels)
    xtickangle(45)
    if savep
        saveas(gcf,[dir2Name,'\',msg,'.jpg'])
    end
    hold off
    
    figure('Visible',figstate);
    HnLList = [4,5,6,7,8,9,10];
    msg = ['All HnL subjects, time evolution of feature ',selFeat];
    tetempFeat = zeros(1,8);
    for f = 1:length(fnames)
        tempFeat = zeros(1,7);
        for h = 1:length(HnLList)
            tempFeat(h) = table2array(condStruct.(fnames{f,1})(rowNames{1,HnLList(h)},selFeat));
        end    
        ci = bootci(100,@mean,tempFeat);
        cil(f) = ci(1);
        ciu(f) = ci(2);
        tetempFeat(f) = mean(tempFeat);
    end
    plot(tetempFeat,'LineStyle','--','LineWidth',1.5,'Marker','o')
    hold on 
    plot(cil,'LineStyle','--','LineWidth',1.5,'Color','r')
    legend('mean','95 % CI','AutoUpdate','off')
    hold on
    plot(ciu,'LineStyle','--','LineWidth',1.5,'Color','r')
    ylim([minF,maxF])
    title(msg)
    ylabel(selFeat)
    xlabel('Conditions Time Evolution')
    xticklabels(time_labels)
    xtickangle(45)
    if savep
        saveas(gcf,[dir2Name,'\',msg,'.jpg'])
    end
    hold off
    
    figure('Visible',figstate);
    msg = ['All HAT subjects, time evolution of feature ',selFeat];
    tetempFeat = zeros(1,8);
    for f = 1:length(fnames)
        tempFeat = zeros(1,3);
        for h = 1:3
            tempFeat(h) = table2array(condStruct.(fnames{f,1})(rowNames{1,h},selFeat));
        end  
        ci = bootci(100,@mean,tempFeat);
        cil(f) = ci(1);
        ciu(f) = ci(2);
        tetempFeat(f) = mean(tempFeat);
    end
    plot(tetempFeat,'LineStyle','--','LineWidth',1.5,'Marker','o')
    hold on 
    plot(cil,'LineStyle','--','LineWidth',1.5,'Color','r')
    legend('mean','95 % CI','AutoUpdate','off')
    hold on
    plot(ciu,'LineStyle','--','LineWidth',1.5,'Color','r')
    ylim([minF,maxF])
    title(msg)
    ylabel(selFeat)
    xlabel('Conditions Time Evolution')
    xticklabels(time_labels)
    xtickangle(45)
    if savep
        saveas(gcf,[dir2Name,'\',msg,'.jpg'])
    end
    hold off
    
    figure('Visible',figstate);
    msg = ['All LAT subjects, time evolution of feature ',selFeat];
    tetempFeat = zeros(1,8);
    for f = 1:length(fnames)
        tempFeat = zeros(1,8);
        for h = 4:11
            tempFeat(h) = table2array(condStruct.(fnames{f,1})(rowNames{1,h},selFeat));
        end    
        ci = bootci(100,@mean,tempFeat);
        cil(f) = ci(1);
        ciu(f) = ci(2);
        tetempFeat(f) = mean(tempFeat);
    end
    plot(tetempFeat,'LineStyle','--','LineWidth',1.5,'Marker','o')
    hold on 
    plot(cil,'LineStyle','--','LineWidth',1.5,'Color','r')
    legend('mean','95 % CI','AutoUpdate','off')
    hold on
    plot(ciu,'LineStyle','--','LineWidth',1.5,'Color','r')
    ylim([minF,maxF])
    title(msg)
    ylabel(selFeat)
    xlabel('Conditions Time Evolution')
    xticklabels(time_labels)
    xtickangle(45)
    if savep
        saveas(gcf,[dir2Name,'\',msg,'.jpg'])
    end
    hold off

end

