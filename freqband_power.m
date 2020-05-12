clear all;close all;
load neuro_act.mat;
x = neuro_act(:,1);
fs = 1000;
fcutlow = 1;
fcuthigh = 15;
n1ms = fs/1000;
frsize = 500*n1ms;
frshift = 25*n1ms;
frovlap = frsize - frshift;

nfft = 1024;
nfftby2 = round(nfft/2 + 1);
hfhz = linspace(0,fs/2,nfftby2); % frequency scale
flidx = round(fcutlow * (nfft/fs)); % index of fcutlow in fft return values
fhidx = round(fcuthigh * (nfft/fs)); % index of fcuthigh in fft return values

% make frames and window the signal
xb = buffer(x,frsize); % divide the signal into non-overlapping frames
xbw = bsxfun(@times,xb,hamming(frsize));
nof = size(xb,2);

% nif = 1000;
% xb = xb(:,1:nif);
% xbw = xbw(:,1:nif);

% compute spectrum 
magy_xb = abs(fft(xb,nfft));
magy_xbw = abs(fft(xbw,nfft));

% half magnitude spectrum
hmagy_xb = magy_xb(1:nfftby2,:);
hmagy_xbw = magy_xbw(1:nfftby2,:);

%     figure();
%     for i = 1:nof
%         subplot(411);
%         plot(xb(:,i));
%         subplot(412);
%         plot(hfhz,hmagy_xb(:,i));  
%         subplot(413);
%         plot(xbw(:,i));
%         subplot(414);
%         plot(hfhz,hmagy_xbw(:,i));
%         pause
%     end

avgpwrb = sum(sum(magy_xb(flidx:fhidx,:).^2));
avgpwrb = avgpwr / nof;
avgpwrbw = sum(sum(magy_xbw(flidx:fhidx,:).^2));
avgpwrbw = avgpwr / nof;
    