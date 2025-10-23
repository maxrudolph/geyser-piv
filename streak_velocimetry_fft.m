clear; close all;
addpath SoundSpeedScript_Karlstrometal2013/
% do simple streak velocimetry
% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000';
% video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000.mp4';
% roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000_roi.mat';
% location = 1;
% load(roi_file);
% roi(1) = 233;
% roi(3) = 150;

% % second video
% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/';
filepath = '~/Box/Geyser Field Experiments/Old Faithful/High Speed Video/04-11-2025/'
% video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000.mp4';
% roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000_roi.mat';
location = 1; % camera location - corresponds to locations in field notebook
load(roi_file);
roi(3) = 150;

% note that this varies video by video depending on location and framerate.
switch location
    case 1
        lscale = 0.0310;        % (meters) = (lscale)*(pixels)
        fs = 120;               % framerate, 1/s
        vscale = lscale*fs;    %(pixel/frame)*(meter/pixel)*(frame/s)
    case 2
        % needs to be updated...

    otherwise
        error('NO!!!');
end

% start by extracting slice from movie
video = VideoReader(video_file);
amount = video.NumFrames;

streak_roi = zeros(roi(3),amount,3,'uint8');
streak = zeros(video.Width,amount,3,'uint8');
row = fix(roi(2)+0.5*roi(4));
for i=1:amount
    img = video.readFrame();
    streak_roi(:,i,:) = img(row,roi(1):roi(1)+roi(3)-1,:);
    streak(:,i,:) = img(row,:,:);
end

figure, imagesc(streak);
figure, imagesc(streak_roi);
%%
bwstreak = rgb2gray(streak_roi);
figure, imshow(bwstreak);
streak1 = adapthisteq(bwstreak);
figure, imshow(streak1);
f = fspecial('gaussian',20,10);
streak2 = streak1;% - imfilter(streak1,f,'symmetric');
streak2 = adapthisteq(streak2);
figure, imagesc(streak2)

% streak2=streak1;% this actually gets used for processing, so make sure it's preprocessed adequately.
%%
% define chunk length in pixels
chunk_length = 120; % number of frames
chunk_stride = 10;%fix(chunk_length/2);
chunk_start = 1:chunk_stride:(amount-chunk_length-1);
% chunk_start = 2000;
% chunk_vel = zeros(size(chunk_start));
% velocities = fliplr([1./linspace(1/vmax,1/0.01,nvel-1) 0])
chunk_velocity = zeros(size(chunk_start));

for i=1:5:length(chunk_start)
    % select the chunk
    % norms = zeros(1,nvel);
    k=1;
    % for vel=velocities
    nrow = size(streak2,1);
    chunk = streak2(:,chunk_start(i):chunk_start(i)+chunk_length-1);
    % take 2D FFT of the chunk.
    % 2D FFT
    F = fftshift(fft2(chunk-mean(chunk(:))));
    P = abs(F).^2; % power spectrum

    % Frequency axes
    f_t = (-chunk_length/2:chunk_length/2-1) / (chunk_length / fs); % Hz
    num_x = size(chunk,1);
    dx = lscale;
    f_x = (-num_x/2:num_x/2-1) / (num_x * dx);                   % cycles/m
    [ftt,fxx] = meshgrid(f_t,f_x);
        
    % remove values near zero frequency
    % P(:,abs(f_t)<0.5)= 0;
    % P(abs(f_x)<0.5,:)  = 0;

    % Find peak power location
    [~, idx_max] = max(P(:));
    [idx_fx, idx_ft] = ind2sub(size(P), idx_max);

    dom_ft = f_t(idx_ft);
    dom_fx = f_x(idx_fx);
    v = dom_fx / dom_ft; % m/s

    % compute a covariance matrix of P.
    Pmean = mean(log10(P(:)));
    Pstd  = std( log10(P(:)));
    Pmask = log10(P) > Pmean%+Pstd; % exclude power values except for top ~5% assuming normal dist.

    weights = P(Pmask)';

    angles = atan2(fxx(Pmask)' ,ftt(Pmask)');
    weights(angles>pi/2 | angles < -pi/2) = [];
    angles(angles>pi/2 | angles < -pi/2) = [];

    % weighted average
    sin_sum = sum( weights.*sin(angles) )/sum(weights);
    cos_sum = sum( weights.*cos(angles) )/sum(weights);
    % unweighted average performed beter empirically:
    sin_sum = sum( sin(angles)/length(angles) );
    cos_sum = sum( cos(angles)/length(angles) );
    average_angle = pi/2-atan2(sin_sum,cos_sum);

    %diagnostic plot
    figure(101);clf;
    subplot(2,1,1);
    P1 = P;
    P1(~Pmask)=0;
    pcolor(ftt,fxx,P1); hold on; shading flat; colorbar(); set(gca,'ColorScale','log');
    plot([0 60*cos(pi/2-average_angle)],[0 60*sin(pi/2-average_angle)],'r')

    subplot(2,1,2);
    pcolor((0:chunk_length-1)/fs,lscale*(0:size(chunk,1)-1),chunk); shading flat;
    colorbar();
    hold on
    for k=0:0.2:1.6
        xdata = [0 2];
        plot(xdata+k,xdata*tan(-average_angle),'r');
    end
    chunk_velocity(i) = tan(-average_angle);
end
%% smoothing / postprocessing of velocity data
v_cutoff = 40; %m/s, maximum allowable velocity
velocities = vscale*chunk_vel(1:i-1);
times = chunk_start(1:i-1);
mask = velocities>v_cutoff;
velocities( mask ) = interp1(times(~mask),velocities( ~mask),times(mask));
velocities = movmedian(velocities,15);
figure, plot(velocities);

%% Thermodynamics
vmax = max(velocities);
xx = logspace(-4,0,100);
P = 0.75;
for ix=1:length(xx)
    [cne(ix),ce(ix)] = SoundSpeedKieffer(xx(ix),P); % YNP is at 0.75 Bar elevation
end
% compute density corresponding to maximum velocity
x_vmax = interp1(ce,xx,max(velocities));
h = XSteam('h_px',P,x_vmax);
rho = XSteam('rho_ph',P,h);
figure, plot(xx,ce,'--');
hold on
plot(xx,cne,'-');
legend('Eqlm','Non Eqlm');
set(gca,'XScale','log','YScale','log');
xlabel('mass fraction vapor')
ylabel('sound speed (m/s');
plot(x_vmax,max(velocities),'r+','MarkerSize',10);



%%
% volume estimate
A = pi*0.4^2; % just a guess of the orifice size
% rho = 16;% an estimate - given thermodynamic constraints...
Q = A*velocities;
volume = cumtrapz(times/fs,Q)*rho/1000;

% height = 0.5*(1./chunk_vel(1:i-1)*vscale).^2/9.81;% 1/2*v^2 = gh
height = 0.5*(velocities).^2/9.81;
figure();
subplot(3,1,1);
imagesc(((1:amount)-1)/fs,((1:video.Width)-1)*lscale,streak);
% imagesc(streak,'XData',[0 amount-1]/fs,'YData',[0 video.Width-1]);
hold on
plot((chunk_start(1:i-1)+chunk_length/2)/fs,height+6.5,'r')
set(gca,'YDir','normal');
ylabel('Height (m)')
colormap gray
h1 = gca;
subplot(3,1,2);
plot(chunk_start(1:i-1)/fs,vscale*chunk_vel(1:i-1));
hold on
plot(times/fs,velocities)
ylabel('Velocity (m/s)')
h2 = gca;
set(gca,'YLim',[0 40])
subplot(3,1,3)
plot(times/fs,volume)
ylabel('Discharge (m^3)')
xlabel('Time (s)')
h3 = gca;
linkaxes([h1 h2 h3],'x')
% set(gca,'XLim',[0 110])

%%
% another idea - assume ballistic trajectory
% dv/dt = -g
% dh/dt = v
% d2h/dt2 = -g
% dh/dt = -g*t + v0
% h(t) = -g*t^2/2 + v0*t + h0
% let h0=0, v0 unknown
% (g/2)*t^2 - v0*t + (h-h0)=0
% t=((v0) +- sqrt(v0^2-4*(g/2)*(h)))/(2*g/2)


