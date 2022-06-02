% Advanced DSP Project
% Effect of Deep Brain Stimulation on Parkinsonian Tremor

clc
clear
% close all

cat = 1;
Min = 2.5;
Max = 12;
Fs = 100;
step = 2;%Subsampling step
nfft = 256;
xt = [0 2 4 6 8 10 12];

surr_thresh = [readtable('threshold_HAT.csv'), readtable('threshold_LAT.csv')];

sample_path1 = 'D:\newLaptop\Project_Final_2.0_29_3_21\Project 2.0\PTES_MATLAB\HATLAT_of\LAT_of\g13r15of.ritF50_s';
sample_path2 = 'D:\newLaptop\Project_Final_2.0_29_3_21\Project 2.0\PTES_MATLAB\HATLAT_of\HAT_of\s8r15of.letM64_b';

sample_path{1} = sample_path1;
sample_path{2} = sample_path2;

names = {'g13r15of','s8r15of'};

fig = figure;
for i = 1:length(sample_path)
    x = load(sample_path{i});
    if contains(sample_path{i},'_b')
        x = 1000*x;
    end     
    x = filtering(x,Fs,Min,Max);
    x = subsampling(x,Fs,Max,step);
    newFs = Fs/step;

    if length(x) > 3000
        x = x(1:3000);
    end
    x = x - mean(x);
    xtime = x./max(abs(x));
    x = x/std(x);
%     x_time = x_time - mean(x_time);
    L = length(x);
    
    if cat == 1
        %% Time
        subplot(2,2,i)
        time = 1:length(xtime);
        time = time/newFs;
        plot(time,xtime)
        ylim([-1.1 1.1])
        xlabel('Time (sec)')
        ylabel('Normalized Amplitude')
        title('Time Series')
        if i == 1
            text(25,1.48,'LAT');
        else
            text(25,1.48,'HAT');
        end

        %% Power Spectrum
        subplot(2,2,i+2)
        Xf = fft(x,2*nfft);
        Xf = abs(Xf(1:nfft)).^2;
        f = 1:nfft;
        f = f*newFs/(2*nfft);
        plot(f,Xf)
        ylim([0 3.5*10^4])
        xlim([0 xt(end)+1])
        xticks(xt)
        xlabel('Frequency (Hz)')
        ylabel('Magnitude')
        title('Power Spectrum')
    else    
        %% Bispectrum
        thres = surr_thresh.(names{i});
        bspec = abs(bispecd(x,nfft,8,L,5)); 
        per = thres/max(max(bspec));
        bspec = abs(bispecd(x,nfft,5,nfft,5));
        bspec(bspec<per*max(max(bspec))) = 0;

        subplot(2,2,i)
        if i == 1
            up_vec = length(0:newFs/nfft:xt(end)+1) + nfft/2 + 1;
            bspec_crop = bspec(nfft/2 + 1:up_vec,nfft/2 + 1:up_vec);
            waxis = linspace(0,xt(end)+1,length(bspec_crop));
            contour(waxis,waxis,bspec_crop)
            caxis(subplot(2,2,i),[0,7]);
            hold on
            plot(waxis,waxis,'--')
            xlabel('f1 (Hz)')
            ylabel('f2 (Hz)')
            title('Bispectrum')
            xticks(xt)
            yticks(xt)
        elseif i == 2
%             low_vec = length(0:newFs/nfft:3) + nfft/2 + 1;
%             up_vec = length(0:newFs/nfft:6) + nfft/2 + 1;
            up_vec = length(0:newFs/nfft:xt(end)+1) + nfft/2 + 1;
            bspec_crop = bspec(nfft/2+1:up_vec,nfft/2+1:up_vec);
%             waxis = linspace(3,6,length(bspec_crop));
            waxis = linspace(0,xt(end)+1,length(bspec_crop));
            contour(waxis,waxis,bspec_crop)
            caxis(subplot(2,2,i),[0,7]); 
            hold on
            plot(waxis,waxis,'--')
            xlabel('f1 (Hz)')
            ylabel('f2 (Hz)')
            xticks(xt)
            yticks(xt)
            title('Bispectrum')
        end
        if i == 2
            b=colorbar; colormap(subplot(2,2,2),jet); cl = caxis;
%             b.Limits = [0,7];
            set(b, 'Position', [0.913 0.58 0.022 0.35]) 
        end
%         h = axes(fig,'visible','off'); 
%         c1 = colorbar(h,'Position',[0.93 0.58 0.022 0.35]);  % attach colorbar to h
%         colormap(c1,'jet')
%         caxis(subplot(2,2,2),[0,7]); 
%         display(max(max(bspec)))
    
        %% Bicoherence
        bic = bicoher(x, nfft, 5, nfft,10);

        subplot(2,2,i+2)
        if i == 1
            up_vec = length(0:newFs/nfft:xt(end)+1) + nfft/2 + 1;
            bic_crop = bic(nfft/2 + 1:up_vec,nfft/2 + 1:up_vec);
            waxis = linspace(0,xt(end)+1,length(bic_crop));
            contour(waxis,waxis,bic_crop)
            caxis(subplot(2,2,i+2),[0,1]);
            hold on
            plot(waxis,waxis,'--')
            xlabel('f1 (Hz)')
            ylabel('f2 (Hz)')
            xticks(xt)
            yticks(xt)
            title('Bicoherency Index')
        elseif i == 2
%             low_vec = length(0:newFs/nfft:3) + nfft/2 + 1;
%             up_vec = length(0:newFs/nfft:6) + nfft/2 + 1;
            up_vec = length(0:newFs/nfft:xt(end)+1) + nfft/2 + 1;
            bic_crop = bic(nfft/2+1:up_vec,nfft/2+1:up_vec);
            bic_crop(1:30,1:5) = 0;
            bic_crop(1:5,1:30) = 0;
%             waxis = linspace(3,6,length(bic_crop));
            waxis = linspace(0,xt(end)+1,length(bic_crop));
            contour(waxis,waxis,bic_crop)
            caxis(subplot(2,2,i+2),[0,1]); 
            hold on
            plot(waxis,waxis,'--')
            xlabel('f1 (Hz)')
            ylabel('f2 (Hz)')
            xticks(xt)
            yticks(xt)
            title('Bicoherency Index')
        end
%         h = axes(fig,'visible','off'); 
%         c2 = colorbar(h,'Position',[0.93 0 0.022 0.35]);  % attach colorbar to h
%         colormap(c2,'jet')
%         caxis(h,[0,1]); 
        if i == 2
            c=colorbar; colormap(subplot(2,2,4),jet);
            c.Limits = [0,1];
            set(c, 'Position', [0.913 0.1 0.022 0.35]) 
        end
        display(max(max(bic)))
    end   
end

f = gcf;
if cat == 1
%     saveas(fig,'C:\Users\xaral\Desktop\Time_PS.jpg')
    exportgraphics(f,'C:\Users\xaral\Desktop\Time_PS.jpg','Resolution',600,'BackgroundColor','none')
else
    exportgraphics(f,'C:\Users\xaral\Desktop\BS_BIC.jpg','Resolution',600,'BackgroundColor','none')
end
    
