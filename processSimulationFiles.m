
% This example demonstrates the use of k-Wave for the reconstruction of a
% two-dimensional photoacoustic wave-field recorded  over a linear array of
% sensor elements  The sensor data is simulated using kspaceFirstOrder2D
% and reconstructed using kspaceLineRecon. It builds on the Homogeneous
% Propagation Medium and Heterogeneous Propagation Medium examples.

clear; close all; clc; rng shuffle;
addpath(genpath(pwd))
DATA_CAST = 'single';     % set to 'single' or 'gpuArray-single' to speed up computations
% DATA_CAST = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

BaseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets\Attenuation\Simulation_23_11_11';
folderNames = {'twoLayers1','twoLayers2','twoLayers3'};
%% For looping simulations

for iSim = 1:length(folderNames)
    %% Loading medium
    load([BaseDir,'\',folderNames{iSim},'\medium_',folderNames{iSim}])

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

    kgrid = kWaveGrid(Nx, dx, Ny, dy);

    %% cuec
    % ======
    c0 = 1540;
    t_end = (Nx*dx)*2/c0;     % [s]
    kgrid.makeTime(c0, [], t_end);
    fs = 1/kgrid.dt;
    % ======
    
    center_depth = 25e-3;
    focal_distance = center_depth;   % center of circle
    focal_number = 2;
    nAperture = (focal_distance/focal_number)/dy;
    nAperture = ratio*floor(nAperture/ratio);
    nApertureEle =  nAperture/ratio;
    nAperture = nApertureEle*ratio;
    % ======

    nLines = floor(Ny/ratio); % vary slightly plm_y to get nLines=128

    bf_data_final = nan(kgrid.Nt ,nLines);
    fprintf("\nProcessing file: %s\n",folderNames{iSim})

    for ii = 1:nLines
        fprintf('Processing line %i out of %i\n',ii,nLines)
        jj = ratio*ii;
        src_ini = max(1,jj - nAperture/2) ;
        src_fin = min(Ny,jj + nAperture/2-1) ;
        sensor_data = h5read(fullfile(BaseDir,folderNames{iSim},...
            'output',[num2str(ii),'.h5']), '/p');
        sensor_data=sensor_data';

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

    %density_map = medium.density;
    %attenuation_map = medium.alpha_coeff;
    save(fullfile(BaseDir,['rf_',folderNames{iSim},'.mat']),...
        'rf','x','z','fs','medium');
end
