clear; close all;
addpath SoundSpeedScript_Karlstrometal2013/
% do simple streak velocimetry
% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000';
% video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000.mp4';
% roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000_roi.mat';
% fs = 120;               % framerate, 1/s
% location = 1;
% load(roi_file);
% roi(1) = 233;

% % second video
% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/';
filepath = '~/Box/Geyser Field Experiments/Old Faithful/High Speed Video/04-11-2025/'
% video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000.mp4';
video_file = '~/Box/Geyser Field Experiments/Old Faithful/High Speed Video/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000.mp4';
% roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000_roi.mat';
roi_file = '~/Box/Geyser Field Experiments/Old Faithful/High Speed Video/04-11-2025/old-faithful-4112028_32216_F2_103_Y20250411H101013.692794000_roi.mat';
fs = 120;
location = 1; % camera location - corresponds to locations in field notebook
load(roi_file);
roi(3) = 300;
g=9.81;

% note that this varies video by video depending on location and framerate.
switch location
    case 1
        lscale = 0.0310;        % (meters) = (lscale)*(pixels)        
        vscale = lscale*fs;     %(pixel/frame)*(meter/pixel)*(frame/s)
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

streak_bw = rgb2gray(streak);
streak_roi_bw = rgb2gray(streak_roi);

streak_bw = adapthisteq(streak_bw);
streak_roi_bw = adapthisteq(streak_roi_bw);

figure, imagesc(streak_bw);
figure, imagesc(streak_roi_bw);

%% 2. Do a color-based segmentation using jet pixels

figure
imshow(streak);
poly=drawpolygon();
mask = poly.createMask;
maskind = find(mask);
[i,j]= ind2sub(size(mask),maskind); % i,j indices
imr = streak(:,:,1);
img = streak(:,:,2);
imb = streak(:,:,3);
pix_r = imr(maskind);
pix_g = img(maskind);
pix_b = imb(maskind);
pix = [pix_r pix_g pix_b];
centroid = mean(pix,1);
[Cpix] = cov(double(pix));
%% compute mahalanobis distance for all pixels in input image
impix = double(reshape(streak,[],3));
impix2 = impix-ones(length(impix),1)*centroid ;
imdist = sqrt(dot(impix2',(inv(Cpix)*impix2')));
imdist = reshape(imdist,size(mask));
figure, histogram(imdist,1000);
threshold = 7.5; % note that this is a tunable parameter.
im_mask = imdist < threshold;
se = strel('disk',10);
im_mask = imopen(im_mask,se);

B = imoverlay(streak,im_mask);
figure, imshow(B);

%% optionally - mask out unwanted areas
B = imoverlay(streak,im_mask);
figure, imshow(B);
title('select mask region to remove')
poly2=drawpolygon();
mask2 = poly2.createMask;
maskind2 = find(mask2);
[i,j]= ind2sub(size(mask2),maskind2); % i,j indices
im_mask2 = im_mask;
im_mask2(maskind2) = false; % remove these pixels from the mask.

B = imoverlay(streak,im_mask2);
figure, imshow(B);
%% estimate plume height from the binary mask
nt = size(streak,2);
height_estimate = zeros(1,nt)*NaN;
vent_estimate = zeros(1,nt)*NaN;
for i=1:nt
    v = im_mask2(:,i);
    if any(v)
        ind = find(v,1,'last');
        height_estimate(i) = ind;
        ind = find(v,1,'first');
        vent_estimate(i) = ind;
    end
end
figure,
imshow(streak);
hold on
plot(height_estimate,'r');
plot(vent_estimate,'g');

%% From the binary mask estimate, compute the exit velocity
% compute the median vent position
median_vent_position = median(vent_estimate(~isnan(vent_estimate)));
% at each time, measure the height
jet_height = (height_estimate - median_vent_position)*lscale;
dt = 1/fs;
jet_time = (0:nt-1)*dt;

% compute the vent velocity at the time of ejection
ind=1;
jet_v0 = [];
jet_time0 = [];
jet_time1 = [];
for i=1:nt
    if ~isnan(jet_height(i))
        jet_v0(ind) = sqrt(2*jet_height(i)*g);
        jet_time0(ind) = jet_time(i) - jet_v0(ind)/g;
        jet_time1(ind) = jet_time(i);
        ind=ind+1;
    end
end
% [jet_time0,i] = sort(jet_time0);
% jet_v0 = jet_v0(i);
figure();
plot(jet_time,jet_height);
title('Jet height vs. time')
figure();
plot(jet_time0,jet_v0,'.');

%% cross-correlation based image analysis
window_length = 240; % number of samples in time direction
window_shift = 1;%window_length/2;
window_height = 128;
window_roi_min = 250;
vertical_spacing = 4;
nslice = 7;
clear v0_values c_values;
for i=1:window_shift:nt-window_length
    window_height_max = max([height_estimate(i),window_roi_min+window_height]);
    slice_heights = fix(linspace(window_roi_min,window_height_max,nslice));
    coefs = zeros(nslice-1,2*window_length+1);
    for j=1:nslice-1
        v1 = streak_bw(slice_heights(j),i:i+window_length);
        v2 = streak_bw(slice_heights(j+1),i:i+window_length);
        v1f = detrend(double(v1));
        v2f = detrend(double(v2));
        lags = -window_length:window_length;
        coefs(j,:) = xcorr(v1f,v2f,window_length,'normalized');        
    end
    % find maximum correlation at -ve lag   
    v = zeros(nslice-1,1)*NaN;
    c=v;
    v_heights = ((slice_heights(1:end-1)+slice_heights(2:end))/2-min(vent_estimate))*lscale;
    for j=1:nslice-1
        [cvals,pos] = findpeaks(coefs(j,:),'MinPeakProminence',0.25);
        clags = lags(pos);
        
        ind1 = find(clags<=0,1,'first');
        lag1 = clags(ind1);
        amp1 = cvals(ind1);
        if ~isempty(lag1)
            % the lag is the dt in the image corresponding to the dy between
            % slices
            dy = (slice_heights(j+1)-slice_heights(j))*lscale; % dy in meters
            v(j) = dy/(-lag1/fs);% dydt is velocity
            c(j) = amp1;
        end
    end
    v0_values(:,i) = (v.^2+2*9.81*v_heights').^(1/2);
    c_values(:,i) = c;
end

%%
vmed = zeros(1,length(v0_values))*NaN;
vmax = zeros(1,length(v0_values))*NaN;
for i=1:length(vmed)
    if any(v0_values(:,i) > 0)
        vmed(i) = median(v0_values(:,i));
        vmax(i) = max(v0_values(:,i));
    end
end

% effective jet height
height_from_velocimetry_med = vmed.^2/2/9.81; %1/2 v^2 = h
height_from_velocimetry_max = vmax.^2/2/9.81; %1/2 v^2 = h

figure, plot(max(vmed,1))

%% Plot various estimates together
vent_height_estimate = median(vent_estimate(~isnan(vent_estimate)))*lscale; %in units of meters

t=1/fs*((1:size(streak,2))-1);
t1=1/fs*((1:length(height_from_velocimetry_med))-1);
y = lscale*(0:size(streak,1)-1);
figure();
subplot(2,1,1)
imagesc(t,y,streak)
hold on
plot(t1,height_from_velocimetry_med+vent_height_estimate);
plot(t1,height_from_velocimetry_max+vent_height_estimate);
set(gca,'YDir','normal');
subplot(2,1,2)
plot(t1,vmed);