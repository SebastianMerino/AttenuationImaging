function [xLeft,xRight] = findFreqBand(x, y, ratio)
% ====================================================================== %
% Function that finds the left and right points in the x axis such that 
% they cut the y axis at ratio*max(y). Used to find the frequency band.
% Inputs:
%   x       Vector of frequencies
%   y       Vector of magnitudes, must be the same size of x
%   ratio   Scalar
%   
% Created on Dec, 2023
% ====================================================================== %
y = y / max(y);
dx = x(2)-x(1);
N = length(y);

[~,imax]=max(y);
if imax == 1 || imax == N
    disp('Error');
    ix0 = 0; ixf = 0;
    return;
end

for iLeft = imax:-1:1
    if y(iLeft) < ratio
        xLeft = x(iLeft) + dx*(ratio-y(iLeft))/(y(iLeft+1) - y(iLeft));
        break;
    end
    if (iLeft == 1)
        disp('Error');
    end
end

for iRight = imax:N
    if y(iRight) < ratio
        xRight = x(iRight) - dx*(y(iRight) - ratio)/(y(iRight) - y(iRight-1));
        break;
    end
    if (iRight == N)
        disp('Error');
    end
end