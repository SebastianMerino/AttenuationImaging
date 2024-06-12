
function attFinal = getIdealAcsInc(center,rInc,acsBack,acsInc,x_ACS,z_ACS,windowSize)
% Inputs 
%   center      2-array of xz coordinates in cm
%   rInc        inclusion radius in cm
%   windowSize  2-array of window size in cm

% cx = 1.85; cz = 1.93;
% rInc = 0.95;
% acsBack = 0.5;
% acsInc = 1;

% Ideal ACS dsitribution
xUp = linspace(x_ACS(1), x_ACS(end),1000);
zUp = linspace(z_ACS(1), z_ACS(end),1000);
[Xup,Zup] = meshgrid(xUp,zUp);
attUp = ones(size(Zup)) * acsBack;
circle_inc = ((Xup-center(1)).^2 + (Zup-center(2)).^2)<= rInc^2;
attUp(circle_inc) = acsInc;

% Smoothing with blocksize
nxDown = round(windowSize(1)/(xUp(2) - xUp(1)));                 % Window size
nzDown = round(windowSize(2)/(zUp(2) - zUp(1))); % Window size
attSmooth = imfilter(attUp,ones([nzDown nxDown])/nzDown/nxDown, ...
    "symmetric","same");

% Downsampling to original resolution
[X,Z] = meshgrid(x_ACS,z_ACS);
attFinal = interp2(Xup,Zup,attSmooth,X,Z);

% figure, tiledlayout(1,3)
% nexttile
% imagesc(xUp,zUp, attUp, [0.3,1.2])
% colormap turbo
% colorbar
% axis image
% 
% nexttile
% imagesc(xUp,zUp, attSmooth, [0.3,1.2])
% colormap turbo
% colorbar
% axis image
% 
% nexttile
% imagesc(x_ACS,z_ACS, attFinal, [0.3,1.2])
% colormap turbo
% colorbar
% axis image

end