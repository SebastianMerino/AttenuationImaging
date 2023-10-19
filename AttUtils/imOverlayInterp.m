function [hF,hB,hColor] = imOverlayInterp(B,SWS,climB,clim,alpha,x,z,ROI,xBm,zBm)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm
%   ROI:    Region of interest

[X,Z] = meshgrid(x,z);
[Xq,Zq] = meshgrid(xBm,zBm);
imgInterp = interp2(X,Z,SWS,Xq,Zq);
[hF,hB,hColor] = imoverlay2(B,imgInterp,climB,clim,alpha,x,z,ROI,xBm,zBm);