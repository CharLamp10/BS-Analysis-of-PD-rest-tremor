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


figstate = 'on';
savep = 1;
alpha = 0.05;
vars = {'bspecEn12h8','bicEn12h8'};%,{'totalBic2h8','power46'};{'bspecEn12h8','bicEn12h8'}
time_labels = {'mOn dOff','mOff dOn','mOff dOff','mOff dOff 15',...
    'mOff dOff 30','mOff dOff 45','mOff dOff 60','mOn dOn'};
fullpath = pwd;
varsplot = {'BspecEnt1','BicEnt1'};%{'BspecEnt1','BicEnt1','TotalBic','Power46'};

%sgtitle(msg)
fig = figure('Visible',figstate);

for v = 1:length(vars)
    selFeat = vars{1,v};

    ron = readtable([fullpath,'\','Features_ron.csv']);
    ref = readtable([fullpath,'\','Features_ref.csv']);
    rof = readtable([fullpath,'\','Features_rof.csv']);
    r15of = readtable([fullpath,'\','Features_r15of.csv']);
    r30of = readtable([fullpath,'\','Features_r30of.csv']);
    r45of = readtable([fullpath,'\','Features_r45of.csv']);
    r60of = readtable([fullpath,'\','Features_r60of.csv']);
    ren = readtable([fullpath,'\','Features_ren.csv']);

    
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
    ax = subplot(2,1,v);
    plot(tetempFeat,'LineStyle','--','LineWidth',2,'Marker','o')
    hold on 
    plot(cil,'LineStyle','--','LineWidth',1.5,'Color','r')
    hold on
    plot(ciu,'LineStyle','--','LineWidth',1.5,'Color','r')
    ylim([minF,maxF])
    ylabel(varsplot{1,v})
    if v ==1
       ylim([3,7])
       ylabel(varsplot{1,v})
    else
       ylim([3,7])
       ylabel([varsplot{1,v},' %'])
    end
    xticklabels(time_labels)
    xtickangle(45)
    hold off
end

if savep
    print(gcf,'time_paper1.png','-dpng','-r600'); % 300 dpi
end