clear;
close all

% load video and set scale:
prefix = '32216_F1_47_Y20241121H131006'
switch prefix
    case '32216_F1_47_Y20241121H131006'
        lscale = .02/98;
        fs=1000;
    case '32216_F2_45_Y20241121H131348'
        lscale = .02/98;
        fs=1000;
    case '32216_F3_46_Y20241121H131931'   
        lscale = .02/98;
        fs=1000;
    case '32216_F4_51_Y20241122H084427'
        lscale = .02/98;
        fs=500;
    case '32216_F5_52_Y20241122H085754'
        lscale=.02/98;
        fs=500;
    case '32216_F6_53_Y20241122H090903'
        lscale=.02/98;
        fs= 400;    
    case '32216_F7_54_Y20241122H093634'
        fs = 500;
        lscale = 0.02/110;
    case '32216_F8_55_Y20241122H094618'
        fs = 400;
        lscale = 0.02/111;
    case '32216_F9_56_Y20241122H095501'
        fs = 500;
        lscale = 0.02/111;
    case '32216_F10_57_Y20241122H095919'
        fs = 500;
        lscale = .02/111; % 2 cm = 111 pixel for 32216_F10_57_Y20241122H095919
    case '32216_F11_59_Y20241122H102242'
        lscale = .02/111;
        fs=400;
end
load(['PIV_results_' prefix '.mat']);

if lscale == 1.0
    warning('lscale not set!!')
end
vscale = lscale*fs;%(pix/frame)*(m/pix)*(frame/s)

%% metadata:

% extract the velocity in the middle of the timeseries
umed = cellfun(@(x) median(reshape(x(:,60:80),[],1)), u_original)*vscale;
ubar = cellfun(@(x) mean(reshape(x(:,60:80),[],1)), u_original)*vscale;
umax = cellfun(@(x) max(reshape(x(:,60:80),[],1)), u_original)*vscale;
umin = cellfun(@(x) min(reshape(x(:,60:80),[],1)), u_original)*vscale;
t = (0:length(umed)-1)/fs;

figure, plot(t,umed);
%% plot a spectrum
y = fft(detrend(ubar));
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
y0 = fftshift(y);         % shift y values
f0 = (-n/2:n/2-1)*(fs/n); % 0-centered frequency range
power0 = abs(y0).^2/n;    % 0-centered power
figure()
plot(f0,power0)
xlabel('Frequency')
ylabel('Power')

%%
window_length = 3*fs;
overlap = 0.95*window_length;

[all_Pxx,frequency,times] = mt_spectrogram(ubar,window_length,overlap,fs);
%% 
f=figure();
f.Position(3:4)=[1664 420];
subplot(2,1,2);
pcolor(times,frequency,all_Pxx);
set(gca,'ColorScale','log')
grid on
colorbar();
shading interp
set(gca,'YScale','log');
ylabel('Frequency (Hz)');
xlabel('Time (s)')
h1 = gca();
subplot(2,1,1)
plot(t,ubar);
hold on
plot(t,umed);
plot(t,umax);
% plot(t,umin);
ylabel('Velocity (m/s)')
legend('average','median','maximum')%,'minimum')
h2 = gca();
linkaxes([h1,h2],'x')
h1.Position(3) = h2.Position(3)
% figure, pcolor(s,f,t)
figname = [prefix '_velocities.pdf'];
exportgraphics(gcf,figname)
savefig([prefix '_velocities.fig']);

uavg = ubar;

save([prefix '_velocities.mat'],'t','ubar','umed','umax');