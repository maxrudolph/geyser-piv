clear;
close all;

% filepath = '/Volumes/GeyserData/Old Faithful/04-11-2025/';
filepath = '~/Box/Geyser Field Experiments/Old Faithful/High Speed Video/04-11-2025/'

% find files on path plus dates with suffix .mp4
[filelist] = dir([filepath '*.mp4']);
filelist = filelist( arrayfun(@(x) x.name(1) ~= '.',filelist));

%% Loop over the video files
nfile = length(filelist);
for i=2:nfile
    success = false;
    while ~success
        nplot = 15;
        ff = fullfile(filelist(i).folder,filelist(i).name);
        v = VideoReader(ff);
        nframes = v.NumFrames;
        % seek to nplot many locations in the file
        img1 = [];
        for j=1:nplot
            ind = floor((j-1)*nframes/nplot)+1;
            img = read(v,ind); % assemble a composite image from nplot-many frames of the video
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
        % selct a sub-ROI for pixel statistics used in masking
        % mask_roi = drawrectangle();
        % mask_roi_pix = fix(mask_roi.Position);
        % mask_roi_pix0 = mask_roi_pix;
        % mask_roi_pix(1) = mod(mask_roi_pix(1),nx);
        % % make a vector of pixel values inside the selected region
        % pixval = [];
        % for j=1:nplot
        %     ind = floor((j-1)*nframes/nplot)+1;
        %     img = read(v,ind);
        %     img = img( mask_roi_pix(1):mask_roi_pix(1)+mask_roi_pix(3) , mask_roi_pix(2):mask_roi_pix(2)+mask_roi_pix(4), : );
        %     if j==1
        %         img2=img;
        %     else
        %         img2 = cat(1,img2,img);
        %     end
        % end
        % figure(3);
        % img2 = permute(img2,[2 1 3]);
        % imshow(img2);

        figure(2);
        img1 = permute(img1,[2 1 3]);
        imshow(img1);
        title('ROI image')
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