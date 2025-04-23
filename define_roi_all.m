clear;
close all;

filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/';
% find files on path plus dates with suffix .mp4
[filelist] = dir([filepath '*.mp4']);

% Loop over the video files
nfile = length(filelist);
for i=1:1%nfile
    success = false;
    while ~success
        nplot = 15;
        ff = fullfile(filelist(i).folder,filelist(i).name);
        v = VideoReader(ff);
        nframes = v.NumFrames;
        % seek to middle of file
        img1 = [];
        for j=1:nplot
            ind = floor((j-1)*nframes/nplot)+1;
            img = read(v,ind);
            if j==1
                img1=img;
            else
                img1 = cat(1,img1,img);
            end
        end
        img1 = permute(img1,[2 1 3]);
        figure();
        imshow(img1);
        set(gca,'YDir','normal');
        axis image;
        % Manually define a ROI for this file
        roi = drawrectangle();
        nx = size(img,1);
        ny = size(img,2);
        roi_pix = fix(roi.Position);% x,y,width,height
        roi_pix(1) = mod(roi_pix(1),nx);
        for j=1:nplot
            ind = floor((j-1)*nframes/nplot)+1;
            img = read(v,ind);
            img = img( roi_pix(1):roi_pix(1)+roi_pix(3) , roi_pix(2):roi_pix(2)+roi_pix(4), : );
            if j==1
                img1=img;
            else
                img1 = cat(1,img1,img);
            end
        end
        clf();

        img1 = permute(img1,[2 1 3]);
        imshow(img1);
        set(gca,'YDir','normal');
        yn = menu('Looks good?','Yes','No');
        if yn == 1
            success = true;
        end
    end % end while ~success loop
    % Save the ROI to disk as filename+'_roi.mat'
    roi_file = [filepath filelist(i).name(1:end-4) '_roi.mat']
    roi = roi_pix([2 1 4 3]); % permute back to the original orientation of the video.
    save(roi_file,"roi");
end