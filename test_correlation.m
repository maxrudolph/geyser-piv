% correlation test

y = linspace(0,5,101);
x = linspace(0,10,200);

[xx,yy] = meshgrid(x,y);
zz = sin(2*pi*(xx-yy/5));%+0.5*rand(size(xx));
figure, imagesc(zz);

Y = fft2(zz);
figure()
imagesc(abs(fftshift(Y)));