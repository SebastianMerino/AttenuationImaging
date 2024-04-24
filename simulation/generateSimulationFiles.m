
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.

clear; close all; clc; rng shuffle;
addpath(genpath(pwd))
DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
% DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

BaseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets\Attenuation\Simulation_23_12_18';
folderNames = {'target1','target2','target3'};

%% For looping simulations

for iSim = 1:length(folderNames)
    if (~exist(fullfile(BaseDir,folderNames{iSim},"input"),"dir"))
        mkdir(fullfile(BaseDir,folderNames{iSim},"input"));
    end
    if (~exist(fullfile(BaseDir,folderNames{iSim},"output"),"dir"))
        mkdir(fullfile(BaseDir,folderNames{iSim},"output"));
    end

    %% Generating grid
    elem_pitch = 0.30e-3;

    pml_x_size = 2*40;                % [grid points]
    pml_y_size = 2*14;                % [grid points]

    % set total number of grid points not including the PML
    Nx = 2*810 - 2 * pml_x_size;      % [grid points]
    Ny = 2*540 - 2 * pml_y_size;      % [grid points]

    PML_size = [pml_x_size pml_y_size];   % size of the PML in grid points TONY

    ratio = 8;
    dy = elem_pitch/ratio;          % grid point spacing in the y direction [m]
    dx = dy;

    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    offset = 5; % to prevent echo top
    %% Medium properties 
    axAxis = (0:Nx-1)*dx;
    latAxis = (0:Ny-1)*dy; latAxis = latAxis-mean(latAxis);
    [rx,rz] = meshgrid(latAxis,axAxis);

    c0 = 1540;
    rho0 = 1000;

    cz = 20e-3 + offset*dx; cx = 0;
    r = 8e-3;
    maskLayer = (rz-cz).^2 + (rx-cx).^2 < r^2;

    % define properties of the layer
    switch iSim
        case 1
            background_map_std = 0.008;
            inclusion_map_std = background_map_std;
            background_alpha = 0.6;       % [dB/(MHz^y cm)]
            inclusion_alpha = 1.2;            % [dB/(MHz^y cm)]
        case 2
            background_map_std = 0.008;
            inclusion_map_std = background_map_std*4;
            background_alpha = 0.6;
            inclusion_alpha = 1.2;
        case 3
            background_map_std = 0.008;
            inclusion_map_std = background_map_std/4;
            background_alpha = 0.6;
            inclusion_alpha = 1.2;
    end

    % Multiplicative maps for background and layer
    % background_map = 1 + background_map_std * randn(Nx,Ny);
    % layer_map = 1 + layer_map_std * randn(Nx,Ny);

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

    %%
    c0 = 1540;
    t_end = (Nx*dx)*2/c0;     % [s]
    kgrid.makeTime(c0, [], t_end);
    fs = 1/kgrid.dt;
    
    center_depth = 25e-3;
    focal_distance = center_depth;   % center of circle
    focal_number = 2;
    nAperture = (focal_distance/focal_number)/dy;
    nAperture = ratio*floor(nAperture/ratio);
    nApertureEle =  nAperture/ratio;

    nAperture = nApertureEle*ratio;

    nLines = floor(Ny/ratio); % vary slightly plm_y to get nLines=128

    bf_data_final = nan(kgrid.Nt ,nLines);

    for ii = 1:nLines
        jj = ratio*ii;
        axis_x = rz(:,1);
        axis_y = rx(1,:);
        disp(['Lines: ',num2str(ii),' de ',num2str(nLines)]);
        src_ini = max(1,jj - nAperture/2) ;
        src_fin = min(Ny,jj + nAperture/2-1) ;

        [temp,pos] = min(abs(axis_x-focal_distance));
        focus_point = [axis_y(jj) axis_x(pos)];

        aperture_point_src = [axis_y(src_ini:src_fin)' axis_x(offset)*ones(src_fin-src_ini+1,1)];

%         figure (6); plot(aperture_point_src(:,1),aperture_point_src(:,2),'sb');
%         hold on; plot(focus_point(:,1),focus_point(:,1),'sr'); hold off;

        %%
        sensor.mask = zeros(Nx, Ny);
        % Need a slight offset here to prevent backwards propagating pulse
        sensor.mask(offset, src_ini:ratio:src_fin) = 1;
        %sensor.directivity_size = 0.2698e-3;
        %sensor.directivity_angle = zeros(size(sensor.mask));
        %sensor.directivity_size = 0.4e-3;
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));

        %%
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(offset,src_ini:src_fin) = 1;
        apod = boxcar(nnz(source.p_mask));

        source_strength = 1e6;
        tone_burst_freq = 6.66e6;        % [Hz]
        tone_burst_cycles = 3.5;

        rho0 = 1000;
        input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
        input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

        %excitation = amp*toneBurst(1/kgrid.dt,6.66*1e6,5);
        %src_exc = apod(:)*excitation(:)';
        src_exc = apod(:)*input_signal(:)';

        angle = 0;
        source.p = src_exc;
        %figure; imagesc(source.p2);
        %source.p = source.p2;

        medium.sound_speed_ref = 1540;
        %%

        input_args = {'PMLInside', false, 'PMLSize', PML_size, ... 
            'DataCast',DATA_CAST, 'DataRecast', true, 'PlotSim',false,...
            'SaveToDisk',fullfile(BaseDir,folderNames{iSim},'input',[num2str(ii),'.h5'])};

        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    end
end
