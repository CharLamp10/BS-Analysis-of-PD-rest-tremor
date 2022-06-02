
%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment
function [r,l,u,d,ur,ul,dr,dl] = find_limits(bspec,peak_x,peak_y,t)
% This function calculates the area of a bispectral peak
%% Inputs:
% bspec    -double matrix. The bispectrum
% peak_x   -double. The x coordinate of the peak
% peak_y   -double. The y coordinate of the peak
%% Outputs:
% r    -double. The right limit of the area with respect to the coordinates
%       of the peak
% l    -double. The left limit of the area with respect to the coordinates
%       of the peak
% u    -double. The upper limit of the area with respect to the coordinates
%       of the peak
% d    -double. The down limit of the area with respect to the coordinates
%       of the peak
% ur   -double. The upper-right limit of the area with respect to the coordinates
%       of the peak
% ul   -double. The upper-left limit of the area with respect to the coordinates
%       of the peak
% dr   -double. The down-right limit of the area with respect to the coordinates
%       of the peak
% dl   -double. The down-left limit of the area with respect to the coordinates
%       of the peak
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------


k = 5;
b = 5;
q = 5;
w = 5;
o = 5;
y = 5;
g = 5;
p = 5;
thresh = t*bspec(peak_x,peak_y);
endd = 60;
for i=1:endd
    if (peak_y + i <= length(bspec)) 
        up = bspec(peak_x,peak_y+i);
        upbefore = bspec(peak_x,peak_y+i-1);
    end
    if (peak_y-i) >= 1
        down = bspec(peak_x,peak_y-i);
        downbefore = bspec(peak_x,peak_y-i+1);
    end
    if (peak_x-i) >= 1
        left = bspec(peak_x-i,peak_y);
        leftbefore = bspec(peak_x-i+1,peak_y);
    end
    if (peak_x+i) <= length(bspec)
        right = bspec(peak_x+i,peak_y);
        rightbefore = bspec(peak_x+i-1,peak_y);
    end
    if (peak_x-i) >=1  && (peak_y+i) <= length(bspec)
        upleft = bspec(peak_x-i,peak_y+i);
        upleftbefore = bspec(peak_x-i+1,peak_y+i-1);
    end
    if (peak_y-i) >= 1 && (peak_x-i) >=1 
        downleft = bspec(peak_x-i,peak_y-i);
        downleftbefore = bspec(peak_x-i+1,peak_y-i+1);
    end
    if (peak_x+i) <= length(bspec) && (peak_y+i) <=length(bspec)
        upright = bspec(peak_x+i,peak_y+i);
        uprightbefore = bspec(peak_x+i-1,peak_y+i-1);
    end
    if (peak_y-i) >= 1 && (peak_x+i) <= length(bspec)
        downright = bspec(peak_x+i,peak_y-i);
        downrightbefore = bspec(peak_x+i-1,peak_y-i+1);
    end   
    if (right<=thresh && k~=0 && rightbefore>=right)
        k = 0;
        r = i-1;
    elseif (right>thresh && k~=0 && rightbefore<right)
        k = 0;
        r = i-1;
    elseif k~= 0 && i == endd
        r = endd;
        k = 0;
    end 
    if (left<=thresh && b~=0 && leftbefore>=left)
        b = 0;
        l = i-1;
    elseif (left>thresh && b~=0 && leftbefore<left)
        b = 0;
        l = i-1;
    elseif b~=0 && i == endd
        l = endd;
        b = 0;
    end 
    if (up<=thresh && q~=0 && upbefore>=up)
        q = 0;
        u = i-1;
    elseif (up>thresh && q~=0 && upbefore<up)
        q = 0;
        u = i-1;
    elseif q~=0 && i == endd
        u = endd;
        q = 0;
    end  
    if (down<=thresh && w~=0 && downbefore>=down)
        w = 0;
        d = i-1;
    elseif (down>thresh && w~=0 && downbefore<down)
        w = 0;
        d = i-1;
    elseif w~=0 && i == endd
        d = endd;
        w = 0;
    end  
    if (upright<=thresh && g~=0 && uprightbefore>=upright)
        g = 0;
        ur = i-1;
    elseif (upright>thresh && g~=0 && uprightbefore<upright)
        g = 0;
        ur = i-1;
    elseif g~=0 && i == endd
        ur = endd;
        g = 0;
    end  
    if (downright<=thresh && y~=0 && downrightbefore>=downright)
        y = 0;
        dr = i-1;
    elseif (downright>thresh && y~=0 && downrightbefore<downright)
        y = 0;
        dr = i-1;
    elseif y~=0 && i == endd
        dr = endd;
        y = 0;
    end 
    if (upleft<=thresh && o~=0 && upleftbefore>=upleft)
        o = 0;
        ul = i-1;
    elseif (upleft>thresh && o~=0 && upleftbefore<upleft)
        o = 0;
        ul = i-1;
    elseif o~=0 && i == endd
        ul = endd;
        o = 0;
    end  
    if (downleft<=thresh && p~=0 && downleftbefore>=downleft)
        p = 0;
        dl = i-1;
    elseif (downleft>thresh && p~=0 && downleftbefore<downleft)
        p = 0;
        dl = i-1;
    elseif p~=0 && i == endd
        dl = endd;
        p = 0;
    end  
end

end