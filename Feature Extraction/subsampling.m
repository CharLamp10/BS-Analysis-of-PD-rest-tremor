%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function x_under = subsampling(x,Fs,max_f,step)
%Subsampling Function
%% Inputs:
% x     -double array. The signal
% Fs    -double. Sampling rate
% max_f -double. Max frequency in x
% step  -double. The downsampling step
%% Outputs:
% x_under   -double array. The resulted signal
%
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


%Find acceptable steps
fn = 2*max_f;
if step > Fs/fn
    error('Step must be <= Fs/fn')
end
x_under = x(1:step:length(x));
end