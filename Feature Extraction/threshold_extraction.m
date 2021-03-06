%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment
% 
% In the following script the calculation of the level of statistical significance 
% of the bispectrum of each signal is calculated.
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


clc
clear
close all
warning('off','all')

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
situations = ["ren","ref","ron","rof","r15of","r30of","r45of","r60of"];
Min = 2.5;
Max = 12;
Fs = 100;
step = 2;%Subsampling step
nfft = 256;
M = 100;
P = 100; 

%% Feature Extraction Procedure
for i = 1:length(situations)
    case_name = convertStringsToChars(situations(i));
    %Access HATLAT_of
    for j = 1:2
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
                    if j == 1
                        [~,~,threshold,~] = surrogateDataFunc(X,Fstemp,length(X),nfft,M,P);
                        matrixThreshold_HAT(count) = threshold;
                        names_HAT{count} = amplitudeDir(k+2).name(1:pos-1);
                    else
                        [~,~,threshold,~] = surrogateDataFunc(X,Fstemp,length(X),nfft,M,P);
                        matrixThreshold_LAT(count) = threshold;
                        names_LAT{count} = amplitudeDir(k+2).name(1:pos-1);
                    end
                    count = count + 1;
                end
            end 
        end             
    end


    writetable(array2table(matrixThreshold_HAT,'VariableNames', names_HAT),['threshold_HAT_',case_name,'.csv']);

    writetable(array2table(matrixThreshold_LAT,'VariableNames', names_LAT),['threshold_LAT_',case_name,'.csv']);
end
