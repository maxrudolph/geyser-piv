clear
close all

% make a synthetic streak pattern
nt=400;
ny=202;

t=linspace(0,2,nt);
y=linspace(0,40,ny);
[tt,yy]=meshgrid(t,y);
v=20;
f=10;
fst = 1/(t(2)-t(1));
fsy = 1/(y(2)-y(1));

zz = 128+64*cos(2*pi*f*(tt-1/v *yy) )+rand(size(tt))*64;

figure, imagesc(t,y,zz); colorbar; set(gca,'YDir','normal')

zz = uint8(zz);


gray = double(zz); % Ensure it's in double format

% Compute 2D FFT and shift zero-frequency to center
F = fftshift(fft2(gray));
tfreq = fst/nt*(-nt/2:nt/2-1);
yfreq = fsy/ny*(-ny/2:ny/2-1);

% fx = linspace(-0.5, 0.5, nt);
% fy = linspace(-0.5, 0.5, ny);
% [FX, FY] = meshgrid(fx, fy);
% freq_radius = sqrt(FX.^2 + FY.^2);

% Compute magnitude spectrum
magF = abs(F);
logMag = log(1 + magF);

% Normalize and display spectrum
figure;
% imshow(logMag, [], 'InitialMagnification', 'fit');
imagesc(tfreq,yfreq,logMag); colorbar();
title('Log Magnitude Spectrum');

% Threshold to find peaks in frequency domain
bw = imbinarize(mat2gray(logMag), 0.8); % adjust threshold as needed
% bw = bwareaopen(bw, 20); % remove small noise

% Get peak locations (pixels)
[iy, ix] = find(bw);
center = size(gray)/2;

% Analyze peaks
angles = atan2d(yfreq(iy) - yfreq(center(1)), tfreq(ix) - tfreq(center(2)));
radii = sqrt((iy - center(1)).^2 + (ix - center(2)).^2);

% Filter out central low-frequency region
valid = radii > 10; % skip low-frequency noise near center
angles = angles(valid);
radii = radii(valid);

% Display results
disp('Detected feature angles (degrees):');
disp(angles);

% Estimate spacing: spacing = image_size / frequency_peak
image_size = mean(size(gray)); % average dimension
frequencies = radii / image_size;
spacings = 1 ./ frequencies;

disp('Estimated spacings (pixels):');
disp(spacings);

% Optional: show peaks over frequency image
hold on;
plot(tfreq(ix(valid)), yfreq(iy(valid)), 'r+');