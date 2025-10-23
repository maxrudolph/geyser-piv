clear; close all;
addpath SoundSpeedScript_Karlstrometal2013/
% do simple streak velocimetry
% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000';
% video_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000.mp4';
% roi_file = '/Volumes/GeyserData/Old Faithful/04-11-2025/32216_1_102_Y20250411H090557.221485000_roi.mat';
fs = 120;               % framerate, 1/s
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
%% cross correlation analysis
step=3;
chunk_size = 100;
chunk_step = chunk_size/2;
ncol = size(streak_roi_bw,1);
clear corrs;
ind1=1;
for i=1:amount-step   % outer loop over frames in the video (columns in streak image)
    ind=1;
    for j=1:chunk_step:ncol-chunk_size % loop over chunks that correspond to position above the bottom of the video.
        chunk_start = j;
        
        v1 = streak_roi_bw(chunk_start:chunk_start+chunk_size-1,i);
        v2 = streak_roi_bw(chunk_start:chunk_start+chunk_size-1,i+step);      
        
        v1f = detrend(double(v1));
        v2f = detrend(double(v2));
        % v1f = v1;
        % v2f = v2;
        % NOTE - a correlation at negative lag corresponds to an upward
        % velocity.
        corr = xcorr(v1f,v2f,chunk_size,'normalized');
        lags=-chunk_size:chunk_size;
        corrs(:,ind,ind1) = corr; % 1st index is across lag, 2nd index is across spatial positions, 3rd index is across time within the series.
        ind=ind+1;
    end
    ind1 = ind1+1;
end
%% in order to extract velocities, we need to stack correlations over X
% frames (temporal averaging):
stack_length = 10;
clear avg_corrs;
for i=1:(size(corrs,3)-stack_length)
    avg_corrs(:,:,i) = sum(corrs(:,:,i:i+stack_length-1),3)/stack_length;
end
figure, imagesc(avg_corrs(:,:,1));
% now loop over time. At each time, find the peak corresponding to the max.
% corr at each vertical position. compute the velocity.
clear vy;
for i=1:size(avg_corrs,3)
    for j=1:size(avg_corrs,2) % loop over spatial position
        [cor,pos] = findpeaks(avg_corrs(:,j,i),"MinPeakProminence",0.3);
        % sort the peaks by correlation        
        [cor1,ind1] = sort(cor);
        pos1 = pos(ind1);
        % cor1 and pos1 contain the sorted peak correlations and indices
        lags1 = lags(pos1); % lags1 is the sorted lags corresponding to the peak positions
        mask = lags1 < 0;
        lags1 = lags1(mask); % sorted lags (by correlation), for nonpositive lags
        cor1 = cor1(mask); % sorted correlations at nonpositive lag (+vy)
        if ~isempty(lags1)
            vmax = lags1(end)/step; % pixels/frame
            vy(j,i) = vmax; % pixels/frame
        else
            vy(j,i) = NaN;
        end
        
    end
end
figure, imagesc(vy); colorbar;
%%
h = (1:chunk_step:ncol-chunk_size)*lscale;
v_avg = zeros(amount-step-stack_length,1)*NaN;
for i=1:amount-step-stack_length
    values = vy(:,i);
    % correct the velocities back to an exit velocity
    v0_values = ((vy(:,i)*vscale).^2+2*9.81*h).^(1/2);
    mask = ~isnan(values);
    if any(mask)
        v_avg(i) = mean(v0_values(mask));
    end
end
figure, plot(v_avg);
% v_avg_ms = v_avg*vscale;
h_jet = 1/2*v_avg.^2/9.81;% 1/2 V^2/g = h

% plot on streak image
figure, 
% subplot(2,1,1);
imagesc(1:amount,(1:size(streak_bw,1)-1)*lscale,streak_bw); colormap("bone"); colorbar();
set(gca,'YDir','normal')
hold on
plot(h_jet+roi(1)*lscale)


%% try optical flow?

% verdict - didn't work.

% % Create an optical flow object
% opticFlow = opticalFlowLK();
% 
% % Estimate the flow field for the entire kymograph
% flow = estimateFlow(opticFlow, streak_bw_roi);
% 
% % The 'flow' object contains the velocity field.
% % We can visualize the vertical component of the flow.
% figure;
% imagesc(flow.Vy); % Vy is the vertical velocity component
% title('Vertical Velocity from Optical Flow');
% xlabel('Time (frames)');
% ylabel('Position (pixels)');
% colorbar;
% h = colorbar;
% ylabel(h, 'Velocity (pixels/frame)');