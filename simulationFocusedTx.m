
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.

clear; close all; clc; rng shuffle;
addpath(genpath(pwd))

% save parameters
BaseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets\Attenuation\Simulation_24_01_03';
mkdir(BaseDir)
folderNames = {'homogeneous1','homogeneous2'};

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 6.66e6;   % source frequency [Hz]
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
ppw             = 6;        % number of points per wavelength, 4 to 8
depth           = 40e-3;    % imaging depth [m]
cfl             = 0.5;      % CFL number, could be 0.5

%% For looping simulations

for iSim = 1:length(folderNames)
    if (~exist(fullfile(BaseDir,folderNames{iSim},"input"),"dir"))
        mkdir(fullfile(BaseDir,folderNames{iSim},"input"));
    end
    if (~exist(fullfile(BaseDir,folderNames{iSim},"output"),"dir"))
        mkdir(fullfile(BaseDir,folderNames{iSim},"output"));
    end

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

    %% MEDIUM
    rz = kgrid.x - translation(1);
    rx = kgrid.y;
    
    background_map_std = 0.008;
    % background_alpha = 0.6;       % [dB/(MHz^y cm)]
    
    switch iSim
        case 1
            background_alpha = 1.2;       % [dB/(MHz^y cm)]
            c0mean = 1.05;
        case 2
            background_alpha = 1.2;
            c0mean = 0.95;
    end
    % Define background properties for each region
    sound_speed_map = c0 * ones(Nx,Ny) .* (c0mean + ...
        background_map_std * randn(Nx,Ny));
    density_map = rho0 * ones(Nx,Ny) .* (1 + background_map_std * randn(Nx,Ny));

    alpha_map = background_alpha + zeros(Nx,Ny);


    medium.sound_speed = sound_speed_map;
    medium.density = density_map;
    medium.alpha_coeff = alpha_map;

    medium.alpha_power = 1.05;
    medium.alpha_mode = 'no_dispersion';
    medium.sound_speed_ref = 1540;

    save(fullfile(BaseDir,folderNames{iSim},['medium_',folderNames{iSim},'.mat']),...
    'medium');
    
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
    
    %% SOURCE
    aperture = source_focus/focal_number;
    element_num = floor(aperture/element_pitch);
    
    % set indices for each element
    ids = (0:element_num-1) - (element_num-1)/2;
    
    % set time delays for each element to focus at source_focus
    time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
    time_delays = time_delays - min(time_delays);
    
    % create time varying source signals (one for each physical element)
    source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, ...
        source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));
    
    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);
    
    % add rectangular elements
    for ind = 1:element_num
        
        % set element y position
        y_pos = 0 - (element_num * element_pitch/2 - element_pitch/2) + (ind-1) * element_pitch;
        
        % add element (set rotation angle to match the global rotation angle)
        karray.addRectElement([0, y_pos], element_width/4, element_width, rotation);
    end

    %% For Looping
    yCords = ( (0:nLines-1) - (nLines-1)/2 )* element_pitch; % Lateral cord of each element
    
    for iLine = 1:nLines
        %% SOURCE  
        % move the array
        translation(2) = yCords(iLine);
        karray.setArrayPosition(translation, rotation)
    
        % assign binary mask
        source.p_mask = karray.getArrayBinaryMask(kgrid);
            
        % assign source signals
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        clc
        disp(['Line: ',num2str(iLine),' of ',num2str(nLines)]);
    
        %% SENSOR
        
        % set sensor mask to record central plane
        sensor.mask = karray.getArrayBinaryMask(kgrid);
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));
        
        %% SIMULATION
        
        % set input options
        input_args = {...
            'PMLInside', false, ...
            'PMLSize', 'auto', ... 
            'DataCast', 'single', ...
            'DataRecast', true, ...
            'PlotSim',false,...
            'SaveToDisk', fullfile(BaseDir,folderNames{iSim},'input',[num2str(iLine),'.h5'])};
        
        % MATLAB CPU
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
            input_args{:});
    end
end
