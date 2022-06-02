%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function [new_bspec,flag,threshold,threshold_per,max_orig,MeanMax] = surrogateDataFunc(x,Fs,L2,nfft,M,P)
% In this function the level of statistical significance of the bispectrum
% of a signal is calculated based on the method proposed by (DOI: 10.1109/TBME.2007.895751)
%% Inputs:
% x     -double array. The signal
% Fs    -double. Sampling rate
% L2    -double. Number of lags for the periodogram
% nfft  -double. Number of lags for the bispectrum
% M     -double. Number of surrogates in each P loop
% P     -double. Number of loops
%% Outputs:
% new_bspec      -double matric. The bispectrum after the thresholding
% flag           -binary. 1 if new_bspec ~= 0 and 0 otherwise
% threshold      -double. The calculated threshold
% threshold_per  -double. The threshold as a percentage of the max
%                 bispectral value
% max_orig       -double. The max bispectral value of x
% MeanMax        -double array. The max value for each P
%
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

N = length(x);
x = (x - mean(x))/std(x);
a = 2*sqrt(Fs/2);

[Gxx,fx] = periodogram(x,[],L2,Fs,'centered');
Gxx = Gxx(length(Gxx)/2+1:length(Gxx));
MeanMax = zeros(1,P);

Bispec = bispecd(x,nfft,8,nfft,5);
max_orig = max(max(abs(Bispec)));

sigma = sqrt((N/4).*Gxx)';
sigma = repmat(sigma,P,1,M);
sigma = shiftdim(sigma,1);
Xr = a*normrnd(0,sigma);
Xi = a*normrnd(0,sigma);
Xk1 = Xr - 1j*Xi;
Xk = cat(1,conj(flip(Xk1,1)),repmat(0+1i*0,1,M,P),Xk1(1:end-1,:,:));
X_surr = ifft(ifftshift(Xk,1),length(x),1,'symmetric');

for p=1:P 
    Bispec_surr = bispecd(X_surr(:,:,p),nfft,8,nfft,5);
    MeanMax(p) = max(max(abs(Bispec_surr))); 
    fprintf('%d \n',p)
end
threshold = BootstrapThresh(MeanMax,1000,0.95);
threshold_per = threshold/max_orig;
if threshold_per>=1
    flag = 0;
else
    flag = 1; %If bspec is not 0
end
new_bspec = abs(Bispec);
new_bspec(find(new_bspec<threshold)) = 0;
end
    
