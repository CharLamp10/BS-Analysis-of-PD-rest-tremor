%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function Ent = bspecEntropy(bspec,upLim1,upLim2,lowLim1,lowLim2,q)
% Calculates the Entropy in a specific frequency block
%% Inputs: 
% bspec    -double matrix. The given bispectrum
% upLim1   -double. upper frequency for the x axis
% upLim2   -double. upper frequency for the y axis
% lowLim1  -double. lower frequency for the x axis
% lowLim2  -double. lower frequency for the y axis
% q        -double. 1 for entropy & 2 for squared entropy
%% Outputs:
% Ent      -double. The calculated entropy
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------




bspec = abs(bspec(lowLim1:upLim1,lowLim2:upLim2)).^q;
p = bspec./sum(bspec,'all');
if isnan(p)
    p = zeros(size(p,1),size(p,2));
    
end
p = p + 10e-10;
Ent = -sum(p.*log(p),'all');

end