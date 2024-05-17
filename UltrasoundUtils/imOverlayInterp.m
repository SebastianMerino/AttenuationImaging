function [hF,hB,hColor] = imOverlayInterp(B,SWS,climB,clim,alpha,x,z,ROI,xBm,zBm)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm
%   ROI:    Region of interest

[X,Z] = meshgrid(x,z);
[Xq,Zq] = meshgrid(xBm,zBm);
imgInterp = interp2(X,Z,SWS,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & ROI;
% [hF,hB,hColor] = imoverlay2(B,imgInterp,climB,clim,alpha,x,z,newRoi,xBm,zBm);


B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(xBm,zBm,B);%
axis image on;
% xlabel('\bf x [mm]')
% ylabel('\bf z [mm]')
colormap(gray)

hColor = colorbar; colormap turbo;
hold on;
hF = imagesc(x,z,imgInterp,clim);

% If images are different sizes, map the front image to back coordinates
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))

% Make the foreground image transparent
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
axis image
