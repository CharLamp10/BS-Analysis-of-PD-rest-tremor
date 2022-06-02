%% Bispectral Analysis of Parkinsonian Rest Tremor: New Characterization
%% and Classification Insights Pre-/Post-DBS and Medication Treatment

function mean_ratio = fit_ellipse(bspec,Fs,nfft,state)
% In this function an ellipse is fitted to the given bispectrum and the
% deviation between this ellipse and the bispectrum is calculated
%% Inputs:
% bspec    -double matrix. The bispectrum
% Fs       -double. The sampling rate
% nfft     -double. The Number of lags for the bispectrum
% state    -binary. 'on' to open figure and 'off' otherwise
%% Outputs:
% mean_ratio  -double. The calculated ratio
% 
%-----------------------------------------------------------------------------------------------------------------
% Authors: Charalampos Lamprou & Ioannis Ziogas
% Copyright (C) 2022 Ioannis Ziogas and Charalampos Lamprou,SPBTU,ECE,AUTh
%-----------------------------------------------------------------------------------------------------------------

p = 5;
bspec(1:5,:) = 0;
bspec(:,1:5) = 0;
[Xold,Yold] = meshgrid(0:Fs/nfft:20+Fs/nfft);

[Xnew,Ynew] = meshgrid(0:1/p*Fs/nfft:20+Fs/nfft);
bspec_interp = interp2(Xold,Yold,bspec,Xnew,Ynew);

figure('Visible',state);
[~,c] = contour(bspec_interp,4);
level_list = c.LevelList;
close
yellow = level_list(end);
green = level_list(end-1);
blue = level_list(end-2);
purple = level_list(end-3);

%Colour levels are obtained from contour plot and split into chromatic
%clusters
[xy,yy] = find(bspec_interp>=yellow);
[xg,yg] = find(bspec_interp>=green & bspec_interp<yellow);
[xb,yb] = find(bspec_interp>=blue & bspec_interp<green);
[xp,yp] = find(bspec_interp>=purple & bspec_interp<blue);


ellipse_similarity_ratio1 = 1*find_ratio([xy;xg;xb;xp],[yy;yg;yb;yp],state);
ellipse_similarity_ratio2 = 1*find_ratio([xy;xg;xb],[yy;yg;yb],state);
ellipse_similarity_ratio3 = 1*find_ratio([xy;xg],[yy;yg],state);
ellipse_similarity_ratio4 = 1*find_ratio([xy],[yy],state);

ellipse_ratios = [ellipse_similarity_ratio1 ellipse_similarity_ratio2 ...
    ellipse_similarity_ratio3 ellipse_similarity_ratio4];

mean_ratio = mean(ellipse_ratios);

    function ellipse_similarity_ratio = find_ratio(xbound,ybound,state)
        figure('Visible',state)
        plot(xbound,ybound,'.')
        xlim([min(xbound)-5,max(xbound+5)])
        ylim([min(ybound)-5,max(ybound+5)])
        k = boundary(xbound,ybound,1);
        hold on
        plot(xbound(k),ybound(k));

        geo_xbound = geomean(xbound(k));
        geo_ybound = geomean(ybound(k));

        Ereal = polyarea(xbound(k),ybound(k));
        Eshape = length(xbound);

        weight1 = Eshape/Ereal;
        if weight1 > 1
            weight1 = 1;
        end

        %Find big axis of ellipse
        maxdist = 0;
        for i = 1:length(xbound)
            xx1 = xbound(i);
            yy1 = ybound(i);
            for j = 1:length(xbound)
                if j ~= i
                   xx2 = xbound(j);
                   yy2 = ybound(j);
                   dist = sqrt((xx1-xx2)^2 + (yy1-yy2)^2);
                   if dist > maxdist
                       maxdist = dist;
                       x1 = xx1;
                       x2 = xx2;
                       y1 = yy1;
                       y2 = yy2;
                   end
                end
            end
        end

        a = maxdist/2;
        b = Ereal/(pi*a);
        t = linspace(0,2*pi,500);
        X = a*cos(t);
        Y = b*sin(t);
        w = atan2(y2-y1,x2-x1);
        x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
        y = (y1+y2)/2 + X*sin(w) + Y*cos(w);

        geox = geomean(x);
        geoy = geomean(y);

        distx = geo_xbound - geox;
        disty = geo_ybound - geoy;

        x = x + distx;
        y = y + disty;

        hold on
        plot(x,y,'b-')
        axis equal
        hold off

        figure('Visible',state)
        x_union = [xbound;x'];
        y_union = [ybound;y'];
        k = boundary(x_union,y_union);
        Eunion = polyarea(x_union(k),y_union(k));
        plot(x_union(k),y_union(k),'b-')
        hold on

        [in,on] = inpolygon(x,y,xbound,ybound);
        [in2,on2] = inpolygon(xbound,ybound,x,y);
        ind = in + on;
        ind2 = in2 + on2;
        ind(find(ind == 2)) = 1;
        ind2(find(ind2 == 2)) = 1;
        ind = ind >= 1 ;
        ind2 = ind2 >= 1 ;

        x_temp = x(ind);
        y_temp = y(ind);
        x_temp2 = xbound(ind2);
        y_temp2 = ybound(ind2);

        x_temp = [x_temp x_temp2'];
        y_temp = [y_temp y_temp2'];

        k = boundary(x_temp',y_temp',0.7);

        Esection = polyarea(x_temp(k),y_temp(k));
        plot(x_temp(k),y_temp(k),'r-')

        if rad2deg(w) > 135
            weight2 = 135/rad2deg(w);
        else
            weight2 = rad2deg(w)/135;
        end
        ellipse_similarity_ratio = weight1*weight2*Esection/Eunion;
        end
end