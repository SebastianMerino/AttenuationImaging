%%
clear,clc
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\', ...
    'Attenuation\DataQUS_4_Merino'];
refDir = [baseDir,'\References\P4-CUELLO-3'];
refFiles = dir([refDir,'\*.mat']);

%% LOADING PARAMETERS
for iRef = 1:length(refFiles)
    load([refDir,'\',refFiles(iRef).name])
    dx = x(2)-x(1);
    dz = z(2)-z(1);
    x = x*1e2; % [cm]
    z = z*1e2; % [cm]
    dynRange = [-60,-10];
    
    blocksize = 20;     % Block size in wavelengths
    c0 = 1540;
    freq_L = 3; freq_H = 9;
    overlap_pc      = 0.8;
    winsize         = 0.5;
    window_type     = 5;
    saran_layer     = 0;
    
    % Region for attenuation imaging
    x_inf = 0.5; x_sup = 3;
    z_inf = 0.5; z_sup = 3;
    
    sam1 = RF;
    
    % Limits for ACS estimation
    ind_x = x_inf <= x & x <= x_sup;
    x = x(ind_x);
    ind_z = z_inf <= z & z <= z_sup;
    z = z(ind_z);
    sam1 = sam1(ind_z,ind_x);
    
    % Freq limits in Hz
    freq_L = freq_L*1e6;   % (Hz)
    freq_H = freq_H*1e6;   % (Hz)
    
    wl = c0/mean([freq_L freq_H]);   % Wavelength (m)
    
    % Number of repetitions of datablocks within a single datablock
    rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
    rpt = round(rpt);
    overlap = 1 - (1/rpt);
    
    
    % Lateral samples between windows
    wx  = round((blocksize*wl)/(dx*rpt));
    
    % Number of lines laterally = repetitions * block without repetition
    nx  = rpt*wx;
    new_blocksize = round(nx*dx/(wl));
    
    L2   = size(sam1,2);                % RF data columns
    ncol = floor((L2-(rpt-1)*wx)/wx);   % Blocksize colums
    sam1 = sam1(:,1:wx*(ncol+rpt-1));
    L2 = size(sam1,2);                  % Actual RF data columns
    x  = x(1:L2);
    xi = 1; xf = L2;
    
    x0 = (xi:wx:xf+1-nx);
    x_ACS = x(1,x0+round(nx/2));
    n  = length(x0);
    
    
    % Axial samples between windows
    wz = floor(nx*dx/(dz*rpt));
    nz = rpt*wz;
    
    % winsize: Percentage of window (max 0.5)
    % nw: Samples of each window axially
    nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;
    L = (nz - nw)*dz*100;   % (cm)
    
    NFFT = 2^(nextpow2(nw)+2);
    band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
    rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
    f  = band(rang)*1e-6; % [MHz]
    p = length(rang);
    
    L1   = size(sam1,1);                % RF data: rows
    nrow = floor((L1-(rpt-1)*wz)/wz);   % Blocksize rows
    sam1 = sam1(1:wz*(nrow+rpt-1),:);
    L1   = size(sam1,1);                % RF data: rows
    z    = z(1:L1);
    zi = 1; zf = L1;
    
    z0 = (zi:wz:zf+1-nz);
    z_ACS = z(z0+round(nz/2));
    m  = length(z0);
    
    z0p = z0 + (nw-1)/2;
    z0d = z0 + (nz-1) - (nw-1)/2;
    
    
    % Plot region of interest B-mode image
    Im=abs(hilbert(sam1));   % envelope calculation
    Im_db=20*log10(Im/max(Im(:)));   % log scale
    figure,
    imagesc(x,z,Im_db); axis image; colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    
    ref{iRef} = sam1;
end

%% Diffraction compensation
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
windowing = windowing*ones(1,nx);

att_ref = attenuation_phantoms_Np(f, 3, []);

% Attenuation reference
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Four 'samples' of the reference phantom
Sp_ref = zeros(m,n,p);
Sd_ref = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        % Reference 1
        sub_block_p = ref{1}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{1}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp1,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd1,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 2
        sub_block_p = ref{2}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{2}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp2,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd2,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 3
        sub_block_p = ref{3}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{3}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp3,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd3,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        % Reference 4
        sub_block_p = ref{4}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{4}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp4,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd4,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);

        tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
        tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));
    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;

save([baseDir,'\References\P4_CUELLO_3.mat'],"diffraction_compensation");