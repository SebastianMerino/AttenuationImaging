% Script to show how to setup a linear array transducer using the
% kWaveArray class.
%

clearvars;

% =========================================================================
% DEFINE LITERALS
% =========================================================================
%%    
% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
model           = 1;

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 6e6;      % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 3.5;      % number of toneburst cycles
source_focus    = 20e-3;    % focal length [m]
element_pitch   = 0.3e-3;   % pitch [m]
element_width   = 0.25e-3;  % width [m]
focal_number    = 2;
nLines          = 96;

% grid parameters
grid_size_x     = 50e-3;    % [m]
grid_size_y     = 40e-3;    % [m]

% transducer position
translation     = [-20e-3, 0];
rotation        = 0;

% computational parameters
DATA_CAST       = 'single'; % set to 'single' or 'gpuArray-single'
ppw             = 4;        % number of points per wavelength, 4 to 8
depth           = 40e-3;    % imaging depth [m]
cfl             = 0.5;      % CFL number

% =====================================================================
% RUN SIMULATION
% =====================================================================
%% GRID


% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(grid_size_x / dx);
Ny = roundEven(grid_size_y / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% create the time array
t_end           = depth*2/c0;     % [s];    % total compute time [s]
kgrid.makeTime(c0, cfl, t_end);

% ------------------------------------------------------------
% ------------------------------------------------------------
%% MEDIUM
rz = kgrid.x - translation(1);
rx = kgrid.y;

cz = 20e-3; cx = 0;
r = 8e-3;
maskLayer = (rz-cz).^2 + (rx-cx).^2 < r^2;

background_map_std = 0.008;
inclusion_map_std = background_map_std/2;
background_alpha = 0.6;       % [dB/(MHz^y cm)]
inclusion_alpha = 1.2;            % [dB/(MHz^y cm)]

% Define properties for each region
sound_speed_map = c0 * ones(Nx,Ny) .* (1 + background_map_std * randn(Nx,Ny));
density_map = rho0 * ones(Nx,Ny) .* (1 + background_map_std * randn(Nx,Ny));
alpha_map = background_alpha + zeros(Nx,Ny);

sound_speed_layer = c0 * ones(Nx,Ny) .* (1 + inclusion_map_std * randn(Nx,Ny));
density_layer = rho0 * ones(Nx,Ny) .* (1 + inclusion_map_std * randn(Nx,Ny));
alpha_layer = inclusion_alpha + zeros(Nx,Ny);      % [dB/(MHz^y cm)]

% assign region
sound_speed_map(maskLayer) = sound_speed_layer(maskLayer);
density_map(maskLayer) = density_layer(maskLayer);
alpha_map(maskLayer) = alpha_layer(maskLayer);

% medium.sound_speed = sound_speed_map;
medium.sound_speed = sound_speed_map;
medium.density = density_map;
medium.alpha_coeff = alpha_map;

medium.alpha_power = 1.05;
medium.alpha_mode = 'no_dispersion';
medium.sound_speed_ref = 1540;
% save(fullfile(BaseDir,folderNames{iSim},['medium_',folderNames{iSim},'.mat']),...
% 'medium');


figure('Units','centimeters', 'Position',[5 5 25 10]),
tiledlayout(1,3),
nexttile,
imagesc(100*rx(1,:),100*rz(:,1),medium.sound_speed)
xlabel('x [cm]'), ylabel('z [cm]')
title('Sound speed')
c = colorbar; ylabel(c,'m/s')
axis image

nexttile,
imagesc(100*rx(1,:),100*rz(:,1),medium.density)
xlabel('x [cm]'), ylabel('z [cm]')
title('Density')
c = colorbar; ylabel(c,'kg/m^3')
colorbar,
axis image

nexttile,
imagesc(100*rx(1,:),100*rz(:,1),medium.alpha_coeff)
xlabel('x [cm]'), ylabel('z [cm]')
title('Absorption')
c = colorbar; ylabel(c,'dB/cm/MHz')
axis image

% ------------------------------------------------------------
% ------------------------------------------------------------
%% SOURCE
aperture = source_focus/focal_number;
element_num = floor(aperture/element_pitch);

% set indices for each element
ids = (0:element_num-1) - (element_num-1)/2;

% set time delays for each element to focus at source_focus
time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
time_delays = time_delays - min(time_delays);

% create time varying source signals (one for each physical element)
source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% add rectangular elements
for ind = 1:element_num
    
    % set element y position
    y_pos = 0 - (element_num * element_pitch / 2 - element_pitch / 2) + (ind - 1) * element_pitch;
    
    % add element (set rotation angle to match the global rotation angle)
    karray.addRectElement([0, y_pos], element_width/4, element_width, rotation);
    
end

% ------------------------------------------------------------
% ------------------------------------------------------------
%% For Looping
bf_data_final = zeros(kgrid.Nt ,nLines);
yCords = ( (0:nLines-1) - (nLines-1)/2 )* element_pitch; % Lateral cord of each element

for iLine = 1:nLines
    %% SOURCE
    % move the array
    translation(2) = yCords(iLine);
    karray.setArrayPosition(translation, rotation)

    % assign binary mask
    source.p_mask = karray.getArrayBinaryMask(kgrid);
    
    % plot the off-grid source mask
    figure, 
    imagesc(kgrid.y(1,:)*1000,kgrid.x(:,1)*1000,single(source.p_mask));
    xlabel('Axial[mm]')
    ylabel('Lateral[mm]')
    axis image
    title('Off-grid source mask');

    % assign source signals
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

    % ------------------------------------------------------------
    % ------------------------------------------------------------
    %% SENSOR
    
    % set sensor mask to record central plane
    sensor.mask = karray.getArrayBinaryMask(kgrid);
    sensor.directivity_size = 10*kgrid.dx;
    sensor.directivity_angle = zeros(size(sensor.mask));
    
    % ------------------------------------------------------------
    % ------------------------------------------------------------
    %% SIMULATION
    
    % set input options
    input_args = {...
        'PMLSize', 'auto', ...
        'PMLInside', false, ...
        'PlotPML', true, ...
        'DisplayMask', 'off'};
    
    % run code
    switch model
        case 1
            % MATLAB CPU
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'single', ...
                'PlotScale', [-0.5, 0.5] * source_amp);
        case 2
            % MATLAB GPU
            sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
                input_args{:}, ...
                'DataCast', 'gpuArray-single', ...
                'PlotScale', [-0.5, 0.5] * source_amp);
            
        case 3
            % C++
            sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
        case 4
            % C++/CUDA GPU
            sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
            
    end
    
    combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);
    %% BEAMFORMING
    fs = 1/kgrid.dt;
    bf_data = BFangle(combined_sensor_data',element_num,fs,c0,element_pitch,...
        'rect',focal_number,0);
    index = floor(element_num/2)+1;
    bf_data_final(:,iLine) = bf_data(:,index);

    close all
end
% =========================================================================
% =========================================================================
%% VISUALISATION
offset = 60;

axAxis = 0:size(bf_data_final,1)-1; axAxis = axAxis*kgrid.dt*c0/2;
latAxis = 0:nLines-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *element_pitch;

rf = bf_data_final(offset:end,:);
axAxis = axAxis(offset:end);
Bmode = db(hilbert(rf));
Bmode = Bmode - max(Bmode(:));

% plot the pressure field 
figure;
imagesc(1e3 * latAxis, 1e3 * axAxis, Bmode, [-50,-5]);
xlabel('z-position [mm]');
ylabel('x-position [mm]');
axis image;
title('Intensity');
cb = colorbar;
title(cb, '[dB]');
colormap gray
