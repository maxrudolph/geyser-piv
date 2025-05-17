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
filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/';
video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000.mp4';
roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000_roi.mat';
location = 1;
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
streak2 = streak1 - imfilter(streak1,f,'symmetric');
streak2 = adapthisteq(streak2);
figure, imagesc(streak2)

% streak2=streak1;% this actually gets used for processing, so make sure it's preprocessed adequately.
%%
% define chunk length in pixels
chunk_length = 30; % at 100 fps...
% define chunk overlap in pixels
vmin = 0.33;
vmax = 10; % frame/pixel
nvel = 100;

chunk_stride = 10;%fix(chunk_length/2);
chunk_start = 1:chunk_stride:(amount-chunk_length-1);
% chunk_start = 2000;
chunk_vel = zeros(size(chunk_start));
% velocities = fliplr([1./linspace(1/vmax,1/0.01,nvel-1) 0])
velocities = linspace(vmin,vmax,nvel);

for i=1:length(chunk_start)
    % select the chunk
    norms = zeros(1,nvel);
    k=1;
    for vel=velocities
        nrow = size(streak2,1);
        chunk = zeros(size(streak2,1),chunk_length);
        for j=1:nrow
            % g = 9.81;% g in units of pixel/frame^2
            % use a quadratic shift
            % if vel ~= 0 && (j-1) <= 1/2*(1/vel)^2% h=1/2*v^2
            % time_shift = ((1/(vel*vscale)) - sqrt((1/(vel*vscale))^2 - 4*(g/2)*(j-1)*lscale ))/(g)/fs;
            % mychunk_start = fix(chunk_start(i) + time_shift);% vel is the velocity in frames/pixel
            % chunk(j,:) = streak1(j,mychunk_start:mychunk_start + chunk_length-1);
            % else
            % chunk(j,:) = 0;
            % end
            % use a linear shift
            time_shift = (1/vel)*(j-1);
            mychunk_start = fix(chunk_start(i) + time_shift);% vel is the velocity in frames/pixel
            chunk(j,:) = streak2(j,mychunk_start:mychunk_start + chunk_length-1);
        end
        % tmp = corrcoef(chunk);
        imsum = sum(chunk,1);

        imdiff = diff(chunk,1);
        % norms(k) = norm(imdiff(:));
        norms(k) = var(imsum); % an image with slanting lines will look smoother in the column-sum. maximize the variance in the column sum to determine optimal shift.
        k = k+1;
        % figure, imagesc(chunk), set(gca,'YDir','normal'), title([num2str(i) ' ' num2str(vel)])
    end

    [n,ind] = max(norms);
    chunk_vel(i) = velocities(ind);
    % make a summary figure
    % tmp = 1./(velocities(ind)*vscale);
    % tshift = ((tmp) - sqrt(tmp^2-4*(g/2)*(0:nrow-1)*lscale))/(2*g/2);

    % figure(); imagesc( streak1(:,chunk_start(i):chunk_start(i)+chunk_length-1) );
    % hold on;
    % plot( ((1:nrow)-1)*velocities(ind),1:nrow,'r');
    % plot( ((1:nrow)-1)*velocities(ind)+125,1:nrow,'r');
    % plot( ((1:nrow)-1)*velocities(ind)+250,1:nrow,'r');
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


