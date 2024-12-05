
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.
function [] = simulateLiver(BaseDir)
close all; clc; rng shuffle;
addpath(genpath(pwd))
% DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

normZero = @(x) x-max(x(:));
rf2Bmode = @(x) 20*log10(abs(hilbert(x)));
%iSim =10;
%%%
for iSim = 10:12
    elem_pitch = 0.30e-3;

    pml_x_size = 2*44;                % [grid points]
    pml_y_size = 2*18;                % [grid points]

    % set total number of grid points not including the PML
    Nx = 2*828 - 2 * pml_x_size;      % [grid points]
    Ny = 2*508 - 2 * pml_y_size;      % [grid points]

    PML_size = [pml_x_size pml_y_size];       % size of the PML in grid points TONY

    ratio = 8;
    dy = elem_pitch/ratio;          % grid point spacing in the y direction [m]
    dx = dy;

    %kgrid = makeGrid(Nx, dx, Ny, dy);
    kgrid = kWaveGrid(Nx, dx, Ny, dy);



    %%

    sd = 0.02;
    medium.sound_speed = 1540;
    %-0*makeDisc(Nx,Ny,0,0,10e-3/dy);
    %medium.sound_speed = medium.sound_speed+1540*sd*randn(size(medium.sound_speed));
    medium.density = 1000-0*makeDisc(Nx,Ny,0,0,10e-3/dy);
    medium.density = medium.density + 1000*sd*randn(size(medium.density));
    load(fullfile(BaseDir,sprintf("dens_%i_2024b.mat",iSim)),"liver_map")
    for j =1:Ny
        for i=1:Nx
            if (liver_map(i,j) == 1)   % Hypoechoic 2nd layer
                medium.density(i,j)= 1060+1000*0.008*randn(1);
            elseif (liver_map(i,j) == 3) % Hyperechoic 1st layer circles
                medium.density(i,j)= 700+1000*0.01*randn(1);
            elseif (liver_map(i,j) == 13) % Interfaces
                medium.density(i,j)= 1100+1000*0.13*randn(1);
            elseif (liver_map(i,j) == 10)
                medium.density(i,j)= 900+1000*0.08*randn(1);
            elseif (liver_map(i,j) == 6) % Liver & first layer
                medium.density(i,j)= 1000+1000*0.01*randn(1);
            else
                medium.density(i,j)= 1000+1000*0.03*randn(1);
            end
        end
    end

    %medium.density = liver_map;
    Nx_tot = Nx;
    Ny_tot = Ny;
    rx = ones(Nx_tot,1)*linspace(-Ny_tot*dy/2,Ny_tot*dy/2,Ny_tot);
    rz = linspace(0,Nx_tot*dx,Nx)'*ones(1,Ny_tot);

    %medium.alpha_coeff = 0.5;


    load(fullfile(BaseDir,sprintf("att_%i_2024b.mat",iSim)),"liver_map")
    medium.alpha_coeff = liver_map;

    %medium.alpha_power = 1.05;
    medium.alpha_power = 1;
    medium.alpha_mode = 'no_dispersion';


    figure('Units','centimeters', 'Position',[5 5 25 10]),
    tiledlayout(1,2),

    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.density)
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Density')
    c = colorbar; ylabel(c,'kg/m^3')
    colorbar,
    axis image

    t3 = nexttile;
    imagesc(100*rx(1,:),100*rz(:,1),medium.alpha_coeff, [0.4 1.7])
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Absorption')
    c = colorbar; ylabel(c,'dB/cm/MHz')
    axis image
    colormap(t3,"turbo")
    %%
    c0 =1540;
    t_end = (Nx*dx)*2/c0;     % [s]
    kgrid.makeTime(c0, [], t_end);
    fs = 1/kgrid.dt;
    %[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed,0.40875,1.2*Nx*dx*2/1540);

    center_depth = 25e-3;
    offset = 5; % to prevent echo top

    focal_distance = center_depth;   % center of circle
    focal_number = 2;
    nAperture = (focal_distance/focal_number)/dy;
    nAperture = ratio*floor(nAperture/ratio);
    nApertureEle =  nAperture/ratio;

    %nApertureEle = 32;
    nAperture = nApertureEle*ratio;

    nLines = floor(Ny/ratio); % vary slightly plm_y to get nLines=128

    bf_data_final = nan(kgrid.Nt ,nLines);
    %%
    for ii = 1:nLines

        jj = ratio*ii;
        disp(['Lines: ',num2str(ii),' de ',num2str(nLines)]);
        src_ini = max(1,jj - nAperture/2) ;
        src_fin = min(Ny,jj + nAperture/2-1) ;
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
        %amp = 100000; % [au]
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(offset,src_ini:src_fin) = 1;
        apod = boxcar(nnz(source.p_mask));

        source_strength = 1e6;
        tone_burst_freq = 6.66e6;        % [Hz]
        tone_burst_cycles = 3.5;

        rho0 = 1000;
        input_signal_norm = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
        input_signal = (source_strength ./ (c0 * rho0)) .* input_signal_norm;

        src_exc = apod(:)*input_signal(:)';

        source.p = src_exc;
        medium.sound_speed_ref = 1540;


        %%
        %PML_alpha = 2;   % Default is 2
        %DATA_CAST       = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
        %input_args = {'PMLInside', false, 'PMLAlpha', PML_alpha, 'PMLSize', PML_size, 'PlotPML', false,...
        %    'Smooth', false,'PlotSim',false, 'DataCast',DATA_CAST, 'DataRecast', true};
        input_args = {'PMLInside', false, 'PMLSize', PML_size, ... %'PMLAlpha', PML_alpha, ...
            'DataCast',DATA_CAST, 'DataRecast', true, 'PlotSim',false};

        % run the simulation
        %colormap gray
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

        sensor_data=sensor_data';

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

    %%
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

    mkdir(BaseDir);
    save(fullfile(BaseDir,sprintf("rf_liver%i.mat",iSim)),...
        "rf","x","z","fs","density_map","attenuation_map","kgrid")
    %return

end
end