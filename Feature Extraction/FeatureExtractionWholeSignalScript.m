%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

% The following script is part of the work of [1]. More specifically the
% feature extraction procedure implemented here as described in the
% manuscript.
%
% [1]:  https://doi.org/10.36227/techrxiv.19589728.v1
% 
% The following features are extracted: 
% BspecEnt1, BspecEnt2, MaxArea TotalArea, NoPeaks, ESR, BicEnt1, BicEnt2
% TotalBic, Power46, Belongs46
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

clc
clear
close all
warning('off','all')
tic

warning('off','MATLAB:MKDIR:DirectoryExists')
scriptDir = pwd;
pathSep = strfind(scriptDir,filesep);
projectDir = scriptDir(1:pathSep(end)-1);

disp(['Accessing ' projectDir])
fprintf('\n')
dataDir = dir(projectDir);

class = fullfile(dataDir(5).folder,dataDir(5).name);
disp(['Accessing ' class])
classDir = dir(class);

%% Initializations
Min = 2.5; % Lower frequecy for the bandpass filter
Max = 12; % Higher frequecy for the bandpass filter
Fs = 100; % Sampling rate
step = 2; % Subsampling step
nfft = 256; %  Number of lags for the Bispectrum
thres = [0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25]; % Parameter for NoPeaks calculation
color = ["yg","yg","yg","y","ygb","yg","yg","yg"]; % Parameter for NoPeaks calculation 


varnames = ["diagBicEn28","maxArea","totalArea",...
            "ellipse_similarity_ratio","noPeaks","belongs46","power46",...
            "bspecEn12h8","bspecEn22h8","totalBic2h8","bicEn12h8","bicEn22h8","class"];
        
situations = ["ren","ref","ron","rof","r15of","r30of","r45of","r60of"];

surr_thresh = [readtable('threshold_HAT.csv'), readtable('threshold_LAT.csv'),...
               readtable('threshold_HAT_ron.csv'), readtable('threshold_LAT_ron.csv')...
               readtable('threshold_HAT_ren.csv'), readtable('threshold_LAT_ren.csv')...
               readtable('threshold_HAT_ref.csv'), readtable('threshold_LAT_ref.csv')];% Import threshold values 
%            calculated by the surrogate data method, check
%            threshold_extraction.m
%% Feature Extraction Procedure
for i = 1:length(situations)
    %Access HATLAT_of
    for j = 1:2
        temp = NaN(8,1);
        matrix_mean_FFT = temp;
        matrix_max_FFT = temp;
        matrix_std_FFT = temp;
        matrixPower = temp;
        matrixArea = temp;
        matrixPeaks = temp;
        matrixbelongs_46 = temp;
        matrixbspecEn128 = temp;
        matrixtotalArea = temp;
        matrixbspecEn228 = temp;
        matrixellipse_similarity_ratio = temp;
        matrix_totalBic28 = temp;
        matrix_diagBicEn28 = temp;
        matrix_bicEn128 = temp;
        matrix_bicEn228 = temp;
        matrix_med_freq = temp;
        clas = temp;

        amplitude = fullfile(classDir(j+2).folder,classDir(j+2).name);
        disp(['Accessing ' amplitude])
        amplitudeDir = dir(amplitude);
        count = 1;
        for k = 1:length(amplitudeDir)-2
            if ~amplitudeDir(k+2).isdir %Ignore folders
                if contains(amplitudeDir(k+2).name ,situations(i))
                    sample = fullfile(amplitudeDir(k+2).folder,amplitudeDir(k+2).name); 
                    fprintf(['Loading sample ''' amplitudeDir(k+2).name ''' ...\n'])

                    X = load(sample); %Load the sample
                    if contains(amplitudeDir(k+2).name,'_b') 
                        X = 1000*X; % meters to milimeters
                    end     
                    X = bandpass(X,[Min Max],Fs);
                    X = subsampling(X,Fs,Max,step);

                    if length(X) > 3000
                        X = X(1:3000); % cropping to ensure equal duration
                    end
                    
                    Fstemp = Fs/step;
                    pos = strfind(convertCharsToStrings(amplitudeDir(k+2).name),'.');
                    [maxArea,totalArea,ellipse_similarity_ratio,number_of_peaks,belongs_46,...
                    power_per_4_6,bspecEn128,bspecEn228,totalBic28,diagBicEn28,bicEn128,bicEn228] ...
                                = FeatureExtractionWholeSignal(X,Fstemp,nfft,thres(i),...
                                    surr_thresh.(amplitudeDir(k+2).name(1:pos-1)),color(i),2.5,8);

                    matrixArea(count) = maxArea;
                    matrixtotalArea(count) = totalArea;
                    matrixellipse_similarity_ratio(count) = ellipse_similarity_ratio;
                    matrixPeaks(count) = number_of_peaks;
                    matrixbelongs_46(count) = belongs_46;
                    matrixPower(count) = power_per_4_6;
                    matrixbspecEn128(count) = bspecEn128;
                    matrixbspecEn228(count) = bspecEn228;                  
                    matrix_totalBic28(count) = totalBic28;
                    matrix_diagBicEn28(count) = diagBicEn28;
                    matrix_bicEn128(count) = bicEn128;
                    matrix_bicEn228(count) = bicEn228;
                    if j == 1
                        clas(count) = 1;
                    else
                        clas(count) = 0;
                    end
                    count = count + 1;
                end  
            end
        end
        Features1 = [matrix_diagBicEn28,matrixArea,matrixtotalArea,...
                    matrixellipse_similarity_ratio,matrixPeaks,...
                    matrixbelongs_46,matrixPower,matrixbspecEn128,matrixbspecEn228,...
                    matrix_totalBic28,matrix_bicEn128,matrix_bicEn228,clas]; 
        if j == 1
            Features_HAT = array2table(Features1, 'VariableNames', varnames);
        elseif j == 2
            Features_LAT = array2table(Features1, 'VariableNames', varnames);
        end
    end
    Features = [Features_HAT;Features_LAT];
    Features = Features(~any(ismissing(Features),2),:);
    writetable(Features,['Features_',convertStringsToChars(situations(i)),'.csv']);
end

