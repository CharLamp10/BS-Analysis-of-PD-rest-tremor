%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function threshold = BootstrapThresh(data,btstrpsize,alpha)
% This function calculates the 95% percentile of the given thresholds with
% a 95% level of statistical significance
%% Inputs:
% data         -double array. The values that will be bootstraped
% btstrpsize   -double. The number of bootstraps
% alpha        -double. The level of statistical significance
%% Output:
% threshold    -double. The resulted threshold.
%
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


B = btstrpsize;
[~,bootsam] = bootstrp(B,@mean,data);
samples = data(bootsam);
samples = sort(samples);
samples = samples(ceil(alpha*length(data)),:);
samples = sort(samples);
threshold = samples(0.95*B);

end

