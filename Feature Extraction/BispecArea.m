%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function [Area,number_of_peaks,NewPeaks,ellipse_similarity_ratio] = BispecArea(bspec,Fs,nfft,t,color,state)
% In this function the area, the number of peaks and the ESR are
% calculated.
%% Inputs:
% bspec    -double matrix. The given bispectrum
% Fs       -double. The sampling rate
% nfft     -double. The number of lags for the bispectrum
% t        -double. Parameter for NoPeaks calculation
% color    -string. Parameter for NoPeaks calculation
% state    -binary. check fit_ellipse.m
%% Outputs:
% Area                       -double. The calculated area
% number_of_peaks            -binary. 1 if more than one peak and 0
%                             otherwise
% NewPeaks                   -double matrix. The resulted peaks
% ellipse_similarity_ratio   -double. The calculated ESR
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


if (max(max(bspec))==0)
    Area = 0;
    number_of_peaks = 0;
    NewPeaks = [0 0];
else
    bspec(1:5,:) = 0;
    bspec(:,1:5) = 0;
    %waxis = [0:(nfft/2-1)]'/(nfft)*100;
    figure('Visible',state);
    [~,c] = contour(bspec,4);
    level_list = c.LevelList;
    yellow = level_list(end);
    green = level_list(3);
    blue = level_list(2);
    %Colour levels are obtained from contour plot and split into chromatic
    %clusters
    [xy,yy] = find(bspec>=yellow);
    [xg,yg] = find(bspec>=green & bspec<yellow);
    [xb,yb] = find(bspec>=blue & bspec<green);
    if strcmp(color,"ygb")
        xcoord = [xb; xg; xy]';
        ycoord = [yb; yg; yy]';
        
    elseif strcmp(color,"yg")
        xcoord = [xg; xy]';
        ycoord = [yg; yy]';

    elseif strcmp(color,"y")
        xcoord = xy';
        ycoord = yy';
        
    else
        error("Color must be: ygb,yg or y")
    end
    coord = [xcoord;ycoord];
    ellipse_similarity_ratio = fit_ellipse(bspec,Fs,nfft,state);

    k1 = 1;
    Peaks = findPeaksTest(xcoord, ycoord, bspec); %findPeaks is applied
    %%%Sort Peaks so they match their order in the coord matrix
    for i = 1:length(coord(1,:))
        for j = 1:length(Peaks(1,:))
            if coord(:,i) == Peaks(:,j)
                SortedPeaks(:,k1) = Peaks(:,j);
                k1 = k1+1;
            end
        end
    end


    for i = 1:length(SortedPeaks(1,:))
        for j = 1:length(SortedPeaks(1,:))
            if i~=j && sqrt((SortedPeaks(1,i)-SortedPeaks(1,j))^2 + (SortedPeaks(2,i)-SortedPeaks(2,j))^2)<2
                SortedPeaks(:,i) = 0;
            end
            if i~=j && SortedPeaks(1,i) == SortedPeaks(1,j) && SortedPeaks(2,i) == SortedPeaks(2,j)
                SortedPeaks(:,i) = 0;
            end
        end
    end
    SortedPeaks(:,all(SortedPeaks == 0 )) = [];
    NewPeaks = SortedPeaks;
    Area = zeros(1,length(NewPeaks(1,:)));
    for j=1:length(NewPeaks(1,:))
        [r,l,u,d,ur,ul,dr,dl] = find_limits(bspec,NewPeaks(1,j),NewPeaks(2,j),t);
        limx = [NewPeaks(1,j),NewPeaks(1,j) + ur,NewPeaks(1,j) + r, NewPeaks(1,j) + dr...
            NewPeaks(1,j), NewPeaks(1,j) - dl, NewPeaks(1,j) - l, NewPeaks(1,j) - ul];
        limy = [NewPeaks(2,j)+ u, NewPeaks(2,j) + ur,NewPeaks(2,j), NewPeaks(2,j) - dr...
            NewPeaks(2,j) - d, NewPeaks(2,j) - dl, NewPeaks(2,j), NewPeaks(2,j) + ul];
        % up,upright,right,downright,down,downleft,left,upleft
        Area(j) = polyarea(limx,limy);
    end
    number_of_peaks = length(Area);
    if number_of_peaks == 1
        number_of_peaks = 1;
    else
        number_of_peaks = 0;
    end
end
end
