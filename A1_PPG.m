% Script that introduces spectral analysis using a PPG trace sampled at @ 125 Hz

close all; clear all;
load bidmc_data.mat			   % load variable 'data'

rec = 4;					   % Recording index
ppg = data(rec).ppg.v;         % PPG signal
Fs  = data(rec).ppg.fs;        % Sampling frequency (125 Hz)

time = (1:numel(ppg))*(1/Fs);

% time-domain representation
subplot(3,1,1),plot(time,ppg),xlabel('Time(s)'),grid,title('PPG trace'),xlim([0 16])

% frequency-domain representation via FFT
L=numel(ppg);
Y = fft(ppg); %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);P1 = P2(1:round(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
%Define the frequency domain f and plot the single-sided amplitude spectrum P1.
f = Fs*(0:round(L/2))/L;
subplot(3,2,3),plot(f,P1),title('Single-Sided Amplitude Spectrum'),xlabel('f (Hz)'),ylabel('|P(f)|'),xlim([0 1])

% frequency-domain representation via periodogram
[pxx,faxis] = periodogram(ppg,hamming(L),L,Fs,'psd') ;
subplot(3,2,4),plot(faxis,pxx),title('periodogram'),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 1])

% WELCH-technique for estimating PSD
WINDOW=256;NOVERLAP=32;NFFT=256;[pxx,faxis] = pwelch(ppg,WINDOW,NOVERLAP,NFFT,Fs,'onesided');
subplot(3,2,5),plot(faxis,pxx),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 10]),title('Welch-based PSD'),legend('32samples-overlap ')

WINDOW=256;NOVERLAP=200;NFFT=256;[pxx,faxis] = pwelch(ppg,WINDOW,NOVERLAP,NFFT,Fs,'onesided');
subplot(3,2,6),plot(faxis,pxx),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 10]),title('Welch-based PSD'),legend('200samples-overlap')


