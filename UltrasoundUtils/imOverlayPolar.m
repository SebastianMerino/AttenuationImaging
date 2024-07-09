function [ax1,ax2] = imOverlayPolar(Bm,ACS,dynRange,attRange,transparency,xPolar,zPolar,xPolarACS,zPolarACS)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm
%   ROI:    Region of interest

omitLines = 0;
% Create axes for Bmode
ax1 = axes;
pcolor(ax1, xPolar(1:end-omitLines,omitLines+1:end-omitLines)*1e2, ...
    zPolar(1:end-omitLines,omitLines+1:end-omitLines)*1e2, ...
    Bm(1:end-omitLines,omitLines+1:end-omitLines))
shading interp
axis equal ij tight

% Create axes for Attenuation image
view(2)
ax2 = axes;
s2 = pcolor(ax2, xPolarACS*100,zPolarACS*100,ACS);
shading interp
axis equal ij tight

% Link them together
linkaxes([ax1,ax2])

% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

% Give each one its own colormap
colormap(ax1,'gray')
clim(ax1,dynRange)
colormap(ax2,'turbo')
clim(ax2,attRange)

% Transparency
ax1.Color = [0,0,0];
s2.FaceAlpha = transparency;

% Then add colorbars and get everything lined up
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);
% cb1 = colorbar(ax1,'Position',[.05 .11 .0675 .815]);
% cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);

end