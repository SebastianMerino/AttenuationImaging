function [hF,hB,hColor] = imoverlay2(B,SWS,climB,clim,alpha,x,z,ROI,xBm,zBm)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm
%   ROI:    Region of interest
B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(xBm,zBm,B);%
axis image on;
% xlabel('\bf x [mm]')
% ylabel('\bf z [mm]')
colormap(gray)

hColor = colorbar; colormap turbo;
hold on;
hF = imagesc(x,z,SWS,clim);

% Make the foreground image transparent
alphadata = alpha.*(ROI);
set(hF,'AlphaData',alphadata);
hold off
axis image
