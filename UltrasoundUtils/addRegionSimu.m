function medium = addRegionSimu(medium, c0, rho0, scatterStd, alpha, mask)
% =======================================================================
% Adds a region to medium with the specified properties, given the mask.
%
% INPUTS:
%   medium      struct with properties from the background. If not
%               created, use []
%   c0          mean SoS
%   rho0        mean density
%   alpha       mean attenuation
%   mask        mask to indicate where to locate the nodule/layer
%
% OUTPUTS:
%   medium     new medium
% =======================================================================

[Nx,Ny] = size(mask);

% Define ROI properties for each region 
sound_speed_region = c0 * ones(Nx,Ny) .* (1 + scatterStd * randn(Nx,Ny));
density_region = rho0 .* (1 + scatterStd * randn(Nx,Ny));
alpha_region = alpha + zeros(Nx,Ny);

if isempty(medium)
    medium.sound_speed = sound_speed_region;
    medium.density = density_region;
    medium.alpha_coeff = alpha_region;
else
    medium.sound_speed(mask) = sound_speed_region(mask);
    medium.density(mask) = density_region(mask);
    medium.alpha_coeff(mask) = alpha_region(mask);
end

end