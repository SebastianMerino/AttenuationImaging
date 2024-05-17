%% Saran Wrap Layer
function [T]=saran_wrap(band)

c_sw=2400; % [m/s] Sound Velocity of Saran Wrap
c_p=1540;  % [m/s] Sound Velocity of Phantom (Equal to Water)
rho_sw=1690; % [Kg/m^3] Density of Saran Wrap
rho_p=1000;  % [Kg/m^3] Density of Phantom (Equal to Water)
z_sw=c_sw*rho_sw; % Impendance of Saran Wrap
z_p=c_p*rho_p;    % Impendance of Phantom (Equal to Water)
L_saran=25*1e-6;   % [m]

alpha_sw=5*(band*1e-6).^1.5; % [Np/m/MHz^1.5]
% Considering phantom with equal properties like water
t=abs(2*z_p./( 2*z_p*cos((2*pi/c_sw*(band)-1j*alpha_sw)*L_saran) + ...
    1j*(z_sw+z_p^2/z_sw)*sin((2*pi/c_sw*(band)-1j*alpha_sw)*L_saran) )).^2;
T=t.^2; % Intensity correction for Intensity of Spectrum

% figure (50); plot(band*1e-6, T); 
% title('Correction factor T(f) for phantom'); 
% xlabel('Frequency [MHz]'); ylabel('Coefficient Factor');

end