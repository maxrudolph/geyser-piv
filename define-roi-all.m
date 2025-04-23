clear;
close all;

filepath = '/Volumes/GeyserData/Old\ Faithful/04-11-2025/';
% find files on path plus dates with suffix .mp4
[filelist] = dir(filepath);

% Loop over the video files

% Manually define a ROI for this file

% Save the ROI to disk as filename+'_roi.mat'