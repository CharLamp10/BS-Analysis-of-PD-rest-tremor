%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function V = bandPower46(z,lastSamp,Fs,nfft)
% In this function the percentage of the bispectral power between 4 and 6
% Hz is calculated
%% Inputs:
% z         -double matrix. The bispectrum
% lastSamp  -double. The sample corresponding to the last frequecy of
%            interest
% Fs        -double. The sampling rate
% nfft      -double. The number of lags of the bispectrum
%% Outputs:
% V    -double. The resulted percentage
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

M = length(z);
x1 = M/2+1;
y1 = M/2+1;

x = x1:x1+lastSamp;
y = y1:y1+lastSamp;
z = z(x,y);

Z = abs(z);
Vtotal = trapz(y,trapz(x,Z.^2,2));

f4 = 3.5;
f6 = 6.5; %7.5
x = floor(f4*nfft/Fs):ceil(f6*nfft/Fs);
y = floor(f4*nfft/Fs):ceil(f6*nfft/Fs);
Z = abs(z(x,y));
V = trapz(y,trapz(x,Z.^2,2));
V = V/Vtotal*100;
end


