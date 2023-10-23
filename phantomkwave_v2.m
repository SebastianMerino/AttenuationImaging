
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.

clear; close all; clc; rng shuffle;
addpath(genpath(pwd))
DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
% DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

BaseDir = '.';
fileNames = {'layered1.mat','layered2.mat'};
%% For looping simulations

for iSim = 1:2
    %% Generating grid
    elem_pitch = 0.30e-3;

    pml_x_size = 2*40;                % [grid points]
    pml_y_size = 2*14;                % [grid points]

    % set total number of grid points not including the PML
    Nx = 2*810 - 2 * pml_x_size;      % [grid points]
    Ny = 2*540 - 2 * pml_y_size;      % [grid points]

    PML_size = [pml_x_size pml_y_size];       % size of the PML in grid points TONY

    ratio = 8;
    dy = elem_pitch/ratio;          % grid point spacing in the y direction [m]
    dx = dy;

    %kgrid = makeGrid(Nx, dx, Ny, dy);
    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    offset = 5; % to prevent echo top
    %% Medium properties 
    rx = ones(Nx,1)*linspace(-Ny*dy/2,Ny*dy/2,Ny);
    rz = linspace(0,Nx*dx,Nx)'*ones(1,Ny);

    % ACTIVATE THIS SECTION FOR INCLUSION
    % radius_disk = 10e-3;
    % center_depth = 25e-3;

    maskLayer = rz > 15e-3 & rz < 30e-3;

    sd1 = 0.02; sd2 = 0.05;
    medium.sound_speed = 1540;
    medium.density = 1000 + zeros(Nx,Ny);
    medium.density = medium.density + 1000*sd1*randn(size(medium.density)).*maskLayer;
    medium.density = medium.density + 1000*sd2*randn(size(medium.density)).*~maskLayer;
    medium.density = single(medium.density);

    switch iSim
        case 1
            medium.alpha_coeff = 0.5 + zeros(Nx,Ny);
            medium.alpha_coeff = medium.alpha_coeff + maskLayer*0.5;
            medium.alpha_coeff = single(medium.alpha_coeff);
        case 2
            medium.alpha_coeff = 0.5 + zeros(Nx,Ny);
            medium.alpha_coeff = single(medium.alpha_coeff);
    end
    figure, tiledlayout(1,2),
    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.density)
    colorbar,
    axis image

    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.alpha_coeff)
    colorbar,
    axis image

    medium.alpha_power = 1.05;
    medium.alpha_mode = 'no_dispersion';
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

    %nApertureEle = 32;
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

        figure (6); plot(aperture_point_src(:,1),aperture_point_src(:,2),'sb');
        hold on; plot(focus_point(:,1),focus_point(:,1),'sr'); hold off;

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
        %source.p1 = steer_delay(src_exc,angle,dy,dt,1540);
        %figure; imagesc(src_exc);
        %figure; imagesc(source.p1);
        %source = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,focus_point);
        %source3 = steer_delay_focus(src_exc,angle,dy,kgrid.dt,1540,aperture_point_src,[0 Inf]);
        source.p = src_exc;
        %figure; imagesc(source.p2);
        %source.p = source.p2;

        medium.sound_speed_ref = 1540;
        %%
        %PML_alpha = 2;   % Default is 2
        %DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
        %input_args = {'PMLInside', false, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size, 'PlotPML', false,...
        %    'Smooth', false,'PlotSim',false, 'DataCast',DATA_CAST, 'DataRecast', true};
        input_args = {'PMLInside', false, 'PMLSize', PML_size, ... %'PMLAlpha', PML_alpha, ...
            'DataCast',DATA_CAST, 'DataRecast', true, 'PlotSim',false};

        sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});

        %%
        max_apert = 64; % elements
        f_number = 2;
        c_bf = 1540;
        bf_data = BFangle(sensor_data,max_apert,fs,c_bf,elem_pitch,'rect',f_number,0);

        if src_ini <= 1
            index = size(bf_data,2) - floor(nApertureEle/2);
        elseif src_fin>= Ny
            index = floor(nApertureEle/2)+1;
        else
            index = floor(nApertureEle/2)+1;
        end

        bf_data_final(:,ii) = bf_data(:,index);

    end

    %% Plotting functions
    normZero = @(x) x-max(x(:));
    rf2Bmode = @(x) 20*log10(abs(hilbert(x)));

    %% Plotting
    axAxis = 0:size(bf_data_final,1)-1; axAxis = axAxis*1/fs*c0/2;
    latAxis = 0:size(bf_data_final,2)-1; latAxis = latAxis-mean(latAxis); latAxis = latAxis *elem_pitch;

    rf = bf_data_final(100:end,:);
    z = axAxis(100:end);
    x = latAxis;
    figure; imagesc(x*1e3, z*1e3, normZero(rf2Bmode(rf)),[-60 0])
    title('B-mode'); colormap gray
    xlabel('Lateral distance (mm)')
    ylabel('Depth (mm)')
    axis image

    density_map = medium.density;
    attenuation_map = medium.alpha_coeff;
    save([BaseDir,'\',fileNames{iSim}],'rf','x','z','fs','density_map',...
        'attenuation_map','kgrid','maskLayer');
end
