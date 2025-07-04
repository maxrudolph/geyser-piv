clear;
close all

load PIV_results_32216_F10_57_Y20241122H095919.mat
%% metadata:
fs = 500;

% extract the velocity in the middle of the timeseries
umed = cellfun(@(x) median(reshape(x(:,60:80),[],1)), u_original);
t = (0:length(umed)-1)/fs;

figure, plot(t,umed);
%% plot a spectrum
y = fft(detrend(umed));
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power

plot(f0,power0)
xlabel('Frequency')
ylabel('Power')

%%
spectrogram(umed,1500,500,fs);
% figure, pcolor(s,f,t)