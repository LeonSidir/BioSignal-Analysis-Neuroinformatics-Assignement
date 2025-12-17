% Script that introduces spectral analysis using a respiration impedance trace sampled at @ 125 Hz

close all; clear all;
load bidmc_data.mat			   			 % load variable 'data'

rec = 4;								 % Recording index
resp = data(rec).ref.resp_sig.imp.v;     % Respiration impedance signal 
Fs  = data(rec).ref.resp_sig.imp.fs;     % Sampling frequency (125 Hz)

time = (1:numel(resp))*(1/Fs);

% time-domain representation
subplot(3,1,1),plot(time,resp),xlabel('Time(s)'),grid,title('Respiration impedance trace'),xlim([0 16])

% frequency-domain representation via FFT
L=numel(resp);
Y = fft(resp); %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);P1 = P2(1:round(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);
%Define the frequency domain f and plot the single-sided amplitude spectrum P1.
f = Fs*(0:round(L/2))/L;
subplot(3,2,3),plot(f,10*log10(P1)),title('Single-Sided Amplitude Spectrum'),xlabel('f (Hz)'),ylabel('|P(f)|'),xlim([0 1])

% frequency-domain representation via periodogram
[pxx,faxis] = periodogram(resp,hamming(L),L,Fs,'psd') ;
subplot(3,2,4),plot(faxis,10*log10(pxx)),title('periodogram'),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 1])

% WELCH-technique for estimating PSD
WINDOW=256;NOVERLAP=32;NFFT=256;[pxx,faxis] = pwelch(resp,WINDOW,NOVERLAP,NFFT,Fs,'onesided');
subplot(3,2,5),plot(faxis,10*log10(pxx)),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 10]),title('Welch-based PSD'),legend('32samples-overlap ')

WINDOW=256;NOVERLAP=200;NFFT=256;[pxx,faxis] = pwelch(resp,WINDOW,NOVERLAP,NFFT,Fs,'onesided');
subplot(3,2,6),plot(faxis,10*log10(pxx)),xlabel('f (Hz)'),ylabel('PSD'),xlim([0 10]),title('Welch-based PSD'),legend('200samples-overlap')
