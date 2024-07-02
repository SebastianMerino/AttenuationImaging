%%
clear; clc;
close all;

%% Extract structures from .mat
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation\' ...
%     'Liver_24_06_28\set2'];
baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\Attenuation\' ...
    'Liver\24_06_28'];
input_f = fullfile(baseDir,'sm02.mat');
load(input_f, 'Trans','P','Receive');

%%  Input and output file
input_f = fullfile(baseDir,'rm4.mat');
output_f = fullfile(baseDir,'ref','r4.mat');

saved_data = load(input_f);
RcvData = cell2mat(saved_data.RcvData);
n_frame = size(RcvData,3); % Frame selector
RcvData = RcvData(:, :, n_frame); % Select frame of RF Data


%% Additional variables
central_freq = Receive(1).demodFrequency*1e6; % Central frequency of pulse
sample_freq = Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
%central_freq = 7.6*1e6;
%sample_freq = 4*central_freq;
n_pulses = P.numRays; % number of pulses
n_elements = Trans.numelements; % number of transducer elements
num_samples = Receive(1).endSample - Receive(1).startSample +1; % number of samples per channel
sound_speed = 1540; %Resource.Parameters.speedOfSound; % [m/s] Sound speed defined for acquisition
wvl = sound_speed/central_freq; % [m] wavelength in meters
scalemm2wvl = 1/wvl;
% Initialize variables
rf_channel = zeros(num_samples , n_elements, n_pulses); 
rx_apods = zeros(1, n_elements, n_pulses);
rf_data = zeros(num_samples, n_pulses);
%tx_apods = zeros(1, n_elements, n_pulses);
%steering_angles = zeros(1, n_pulses);

%% Organize data
for n = 1:n_pulses % Iterate through pulses
    rf_channel(:, :, n) = RcvData(Receive(n).startSample:Receive(n).endSample, :); % RF Data from Buffer
    %rx_apods(:, :, n) = TX(n).Apod; % Reception apodizations the same as transmission
    %tx_apods(:, :, n) = TX(n).Apod; % Transmission apodizations
end
%rx_apods = reshape(rx_apods, [size(rx_apods, 2), size(rx_apods, 3)]);
%tx_apods = reshape(tx_apods, [size(tx_apods, 2), size(tx_apods, 3)]);

%% Acommodate to time delays and rf signals
focus = 20/1000;

t = (0:(num_samples-1))/sample_freq; % [sec.] time domain 0:T:(N_sample-1)*T
[rx_delays] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl);

%% Dynamic Aperture

f_num = 4;
z = sound_speed*t/2;
elem_pitch = Trans.spacingMm*1e-3;
maxAprSz = 32;

dyn_aperture = zeros(length(z), n_elements, n_pulses);

for n = 1:n_pulses
    for z_i = 1:length(z)
        a = z(z_i)/(2*f_num);
        hlfAprSz = floor(a / elem_pitch);
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        a_i = -hlfAprSz: hlfAprSz;    % aperture indices
        fulAprSz = 2*hlfAprSz + 1;

        aper_center = n;
        aper = aper_center + a_i;

        aper = aper(aper>=1);
        aper = aper(aper<=128);

        
        dyn_aperture(z_i, aper, n) = 1;
    end
end


%% Beamforming
% Delay-and-sum
for n = 1:n_pulses
    % Delaying
    for e = 1:n_elements
        rf_channel(:, e, n) = interp1(t, rf_channel(:, e, n), rx_delays(:,e, n), 'linear', 0);%rx_apods(e, n);
        %rx_samples = round(rx_delays(:,e, n)*sample_freq);
        % rx_samples(rx_samples<1) = 1;
        % rx_samples(rx_samples>num_samples) = num_samples;
        % rf_channel(:, e, n) = rf_channel(rx_samples, e, n);
        
        %delay_channel = round((2*focus - rx_delays(e,n))/sound_speed*sample_freq);
        %rf_channel(:, e, n) = circshift(rf_channel(:, e, n), delay_channel) * rx_apods(e, n);
    end
    rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
    % Summing
    rf_data(:, n) = sum(rf_channel(:, :, n), 2);
end

%% B-Mode
b_mode = 20*log10(abs(hilbert(rf_data)));
b_mode = b_mode-max(b_mode(:));

param = getparam('C5-2v');
siz = size(rf_data);
z = sound_speed*t/2;
zmax = z(end);
R = param.radius;
p = param.pitch;
N = param.Nelements;
L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
d = sqrt(R^2-L^2/4); % apothem
z0 = -d;

th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
r = linspace(R+p,-z0+zmax,siz(1));

clim_min = -60;
figure();
x = linspace(1, 128, 128); % Same as verasonics image
imagesc(th, r*100, b_mode);
%ylim([0 80])
colormap("gray");
colorbar;
clim([clim_min, 0]);

%% To Polar Coordinates
[xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);

figure();
pcolor(xPolar*1e2,zPolar*1e2,b_mode)
colorbar;
clim([clim_min, 0]);
colormap gray
title('Bmode image')
ylabel('[cm]')
shading interp
axis equal ij tight
%ylim([0 80])
%xlim([-50 50])
%set(gca,'XColor','none','box','off')

%% Saving
rf = rf_data;
fs = sample_freq;
save(output_f,'rf','th','r','xPolar','zPolar',"z0Polar",'fs');

%% Auxiliary functions

% Get delays

function [t_delay] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)
    
    t_delay = zeros(length(t), n_elements, n_pulses);

    element_pos_x = Trans.ElementPos(:, 1)*wvl; % (x, z) [m] Obtain positions of center of every element
    element_pos_z = Trans.ElementPos(:, 3)*wvl;
    phi = Trans.ElementPos(:, 4);

    for n = 1:n_pulses
        for e = 1:n_elements
            focus = sound_speed*t(:)/2;
            xfocus = element_pos_x(n) + sin(phi(n)) * focus;
            zfocus = element_pos_z(n) + cos(phi(n)) * focus;
          %  keyboard 
            %focus = (20/1000)*ones(size(t));
            %t_delay(:,e,n) = (focus + sqrt((focus - element_pos_z(e)).^2 + (element_pos_x(e)-element_pos_x(n))^2))/sound_speed; 
            t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
            
        end
    end

end