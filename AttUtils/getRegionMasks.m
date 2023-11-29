function [back,inc] = getRegionMasks(x,z,c1x,c1z,L,d,Lz)
% ========================================================================
% Function that gives data from a square in the center of an image 
% (inclusion) and two rectangles, one at each side (background)
%
% INPUTS:
%   x,z:        Coordinates of the image. Must be the same size of img.
%   c1x,c1z:    Coordinates of the center of the square
%   L:          Square length
%   d:          Distance from between regions
%
% OUTPUTS:
%   back:       Binary mask of background as a column vector
%   inc:        Binary mask of inclusion as a column vector
% ========================================================================

% Upper left corner of inclusion
x0 = c1x - L/2; 
z0 = c1z - Lz/2;

% Mask inclusion
[X,Z] = meshgrid(x,z);
maskZ = Z>z0 & Z<z0+Lz;
maskInc = X>x0 & X<x0+L & maskZ;
inc = maskInc;

% Upper left corner of each background rectangle
xb1 = x0 - d - L/2;
xb2 = x0 + L + d;

% Mask background
maskBack1 = X>xb1 & X<xb1+L/2 & maskZ;
maskBack2 = X>xb2 & X<xb2+L/2 & maskZ;
back = maskBack1 | maskBack2;


% Vertical squares
% [X,Z] = meshgrid(x,z);
% maskX = X>x0 & X<x0+L;
% maskInc = Z>z0 & Z<z0+L & maskX;
% inc = maskInc;
% zb1 = z0 - d - L/2;
% zb2 = z0 + L + d;
% maskBack1 = Z>zb1 & Z<zb1+L/2 & maskX;
% maskBack2 = Z>zb2 & Z<zb2+L/2 & maskX;
% back = maskBack1 | maskBack2;

% ========================================================================
end