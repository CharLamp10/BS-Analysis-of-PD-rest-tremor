%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function [maxArea,totalArea,ellipse_similarity_ratio,number_of_peaks,belongs_46,...
     power_per_4_6,bspecEn128,bspecEn228,totalBic28,bicEn128,bicEn228]...
     = FeatureExtractionWholeSignal(signal,Fs,nfft,thres,surr_thresh,color,lowLim,upLim)
 % In this function takes place the calculation of all features.
 %% Inputs:
% signal        -double array. The signal
% Fs            -double. Sampling rate
% nfft          -double. The number of lags for the bispectrum
% thres         -double. Parameter for NoPeaks calculation
% surr_thres    -double. The threshold of the bispectrum after surrogate
%                data
% color         -string. Parameter for NoPeaks calculation
% lowLim        -double. Lower frequency of interest
% upLim         -double. Higher frequency of interest
%% Outputs:
% All features as described in (https://doi.org/10.36227/techrxiv.19589728.v1)
%
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


signal = (signal - mean(signal))/std(signal);
L = length(signal);
bspec = abs(bispecd(signal,nfft,8,L,5)); 
per = surr_thresh/max(max(bspec));
bspec = abs(bispecd(signal,nfft,5,nfft,5));
bspec(bspec<per*max(max(bspec))) = 0;


int = Fs/nfft;
lastFreq = 20;
lastSamp = ceil(lastFreq/int);
up_vec = length(0:Fs/nfft:20) + nfft/2 + 1; %Crop image to have only 10 Hz band
upLim = ceil(upLim*nfft/Fs);
lowLim = ceil(lowLim*nfft/Fs);


if max(max(bspec)) ~= 0     
    new_bspec_crop = bspec(nfft/2+1:up_vec,nfft/2+1:up_vec);
    new_bspec_cut = bspec(nfft/2+1:end,nfft/2+1:end);
    bspecEn128 = bspecEntropy(new_bspec_cut,upLim,upLim,lowLim,lowLim,1);
    bspecEn228 = bspecEntropy(new_bspec_cut,upLim,upLim,lowLim,lowLim,2);
    [Area,number_of_peaks,NewPeaks,ellipse_similarity_ratio] = BispecArea(new_bspec_crop,Fs,nfft,thres,color,'off');
    freqPeaks = (NewPeaks-1)*Fs/nfft;
    bspecAmp = diag(bspec(nfft/2+NewPeaks(1,:),nfft/2+NewPeaks(2,:)));
    [~,ind] = max(bspecAmp);
    maxfreqPeaks = freqPeaks(:,ind);
    if maxfreqPeaks(1) > 3.9 && maxfreqPeaks(1) < 6.1 && maxfreqPeaks(2) > 3.9 && maxfreqPeaks(2)
        belongs_46 = 1;
    else
        belongs_46 = 0;
    end
    maxArea = Area(ind);
    totalArea = sum(Area);
    bic = bicoher(signal, nfft, 5, nfft,10);
    bic_cut = bic(nfft/2+1:end,nfft/2+1:end);
    bic_cut(bic_cut > 1) = 1;
    
    bicEn128 = bspecEntropy(bic_cut,upLim,upLim,lowLim,lowLim,1);
    bicEn228 = bspecEntropy(bic_cut,upLim,upLim,lowLim,lowLim,2);
       
    totalBic28 = sum(bic_cut(lowLim:upLim,lowLim:upLim),'all');
    M = upLim - lowLim + 1;
   
    power_per_4_6 = bandPower46(bspec,lastSamp,Fs,nfft);
end