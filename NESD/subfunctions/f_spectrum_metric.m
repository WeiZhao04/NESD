function [spec,fft_regressor] = f_spectrum_metric(X,TR)
%  To perform FFT and compute the spectrum-related metrics including ultra
% low spectrum and ectras.
% X:   (matrice)the iuput signals, should be reformed into timepoints x
%      samples, samples can be 1
% TR:  (int)TR in seconds, needed for FFT
%
%  Outputs
% Spec:(struct) a struct containing metrics, only ul is used for
%      denoise and others are used quality metrics for now
% fft_regressor:(matrice) the full-band fft of signals, should be fed into
%      later denoise functions       


%% basic parameters
fs = 1/TR;             % Sampling frequency
[L,N] = size(X);         % Length of signal
%% basic parameters
if L < 1024
    NFFT = 1024;
else
    NFFT = 2^nextpow2(L);     % Number of the point of FFT
end
freqBin = fs/NFFT;        % Resolution of the frequency                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       /
f = freqBin*(0:(NFFT/2)); % Index of frequency figure
%%
% paddning zeros first
Nzeros = NFFT-L;
X1 = [X; zeros(Nzeros,N)];
Y = fft(X1,NFFT);
power = 2*abs(Y(1:NFFT/2+1,:))/L;
%% ultra low/high power spectrum (UL:drift etc., UH: heartbeat rate etc.)
Lind = find(f>0.01,1,'first')-1;
fALFFind = find(f>0.08,1,'first');
Hind = find(f>0.15,1,'first');
fft_ul = complex(zeros(NFFT,N));
fft_ul(1:Lind,:) = Y(1:Lind,:);
fft_uh = zeros(NFFT,N);
fft_uh(Hind:end,:) = Y(Hind:end,:);
fft_falff = zeros(NFFT,N);
fft_falff(Lind:fALFFind,:) = Y(Lind:fALFFind,:);
%% abnormol power in specific band (CR:cardiac rate, CO2, VA:vasomotion)
CR_Lind = find(f>0.025,1,'first');
CR_Hind = find(f>0.035,1,'first');
CO2_Hind = find(f>0.045,1,'first');
VA_Lind = find(f>0.095,1,'first');
VA_Hind = find(f>0.105,1,'first');

% fft_cr = zeros(NFFT,N);
% fft_cr(CR_Lind:CR_Hind,:) = Y(CR_Lind:CR_Hind,:);
% fft_co2 = zeros(NFFT,N);
% fft_co2(CR_Hind+1:CO2_Hind,:) = Y(CR_Hind+1:CO2_Hind,:);
% fft_va = zeros(NFFT,N);
% fft_va(VA_Lind:VA_Hind,:) = Y(VA_Lind:VA_Hind,:);
%% compute the ratio
% power_fullband = sum(power);
spec.UL = sum(power(1:Lind,:));
spec.UH = sum(power(Hind:end,:));
spec.fALFF = sum(power(Lind:fALFFind,:))./sum(power(Lind:VA_Hind,:));
spec.CR = sum(power(CR_Lind:CR_Hind,:));
spec.CO2 = sum(power(CR_Hind+1:CO2_Hind,:));
spec.VA = sum(power(VA_Lind:VA_Hind,:));
%%
fft_regressor = cat(3,fft_ul);%,fft_cr,fft_co2,fft_va,fft_falff,fft_uh
%% plot figure
% figure('visible','on')
% subplot(2,1,1)
% plot(X,'linewidth',2)
% title('Temporal Course','fontsize',10)
% set(gca,'fontsize',10)
% xlabel('Time')
% ylabel('Amplitude')
% subplot(2,1,2)
% plot(f,Y,'linewidth',2)
% set(gca,'fontsize',10)
% title('Power Spectrum')
% xlabel('Frequency Hz')
% ylabel('Amplitude')

%%

