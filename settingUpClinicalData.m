clear,clc
addpath('./functions_v7');

% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets' ...
%     '\Attenuation\DataQUS_4_Merino'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\DataQUS_4_Merino'];

targetDir = [baseDir,'\Carcinoma'];
refDir = [baseDir,'\References\P4-CUELLO-3'];

% targetDir = [baseDir ,'\TD_31456\31456'];
% refDir = [targetDir,'\ref'];

croppedDir = [targetDir,'\cropped'];
if (~exist(croppedDir,"dir")), mkdir(croppedDir); end

rawFiles = dir([targetDir,'\*.rf']);
refFiles = dir([refDir,'\*.rf']);

%% Generating .mat data
for iAcq = 1:length(rawFiles)
    [out]=lectura_OK([targetDir,'\',rawFiles(iAcq).name]);
    save([targetDir,'\',rawFiles(iAcq).name(1:end-3),'.mat'],'-struct','out');
end
for iAcq = 1:length(refFiles)
    [out]=lectura_OK([refDir,'\',refFiles(iAcq).name]);
    save([refDir,'\T',num2str(iAcq),'.mat'],'-struct','out');
end
%% Cropping data
targetFiles = dir([targetDir,'\*.mat']);
for iAcq = 1:length(rawFiles)
    load(fullfile(targetDir,targetFiles(iAcq).name));
    disp(['Target: ', rawFiles(iAcq).name]);
    attRange = [0.3,1.2];

    dx = x(2)-x(1);
    dz = z(2)-z(1);
    x = x*1e2; % [cm]
    z = z*1e2; % [cm]

    sam1 = RF(:,:,1);
    dynRange = [-50,0];

    % Plotting
    % [pSpectrum,f] = pwelch(RF,hamming(500),300,1024,fs);
    % meanSpectrum = mean(pSpectrum,2);
    % figure, plot(f,meanSpectrum)
    % xlim([0 fs/2])
    % title('Mean power spectrum')
    %
    % Bmode = db(hilbert(sam1));
    % Bmode = Bmode - max(Bmode(:));
    % figure, imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
    % hb2=colorbar; ylabel(hb2,'dB')
    % xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


    %% SETTING PARAMETERS
    blocksize = 15;     % Block size in wavelengths
    c0 = 1540;
    %freq_L = 2; freq_H = 9;
    freq_L = 3; freq_H = 8;
    overlap_pc      = 0.8;
    winsize         = 0.5;
    % Region for attenuation imaging
    % x_inf = 0.1; x_sup = 3.8;
    % z_inf = 0.2; z_sup = 3.1; % 3.5 if discarting T3
    x_inf = 0.5; x_sup = 3;
    z_inf = 0.5; z_sup = 2.1; % 3.5 if discarting T3

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
    Bmode = db(hilbert(sam1));
    Bmode = Bmode - max(Bmode(:));
    figure, imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
    % hb2=colorbar; ylabel(hb2,'dB')
    % xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
        num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
    disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
    disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
    disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
    disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);
    disp(' ')

    save(fullfile(croppedDir,targetFiles(iAcq).name));
end

%% Generating references
x1 = x; z1 = z;
for iRef = 1:length(refFiles)
    load([refDir,'\T',num2str(iRef),'.mat'],'RF','x','z');
    x = x*1e2; z = z*1e2; % [cm]
    ind_x = x_inf <= x & x <= x_sup;
    x = x(ind_x);
    ind_z = z_inf <= z & z <= z_sup;
    z = z(ind_z);
    x = x(1:L2);
    z = z(1:L1);
    fprintf("Checking difference, x: %f, z: %f\n",norm(x1-x), norm(z1-z))
    samRef = RF;
    samRef = samRef(ind_z,ind_x);
    samRef = samRef(:,1:wx*(ncol+rpt-1));
    samRef = samRef(1:wz*(nrow+rpt-1),:);
    ref{iRef} = samRef;
end

% Diffraction compesation
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
windowing = windowing*ones(1,nx);

att_ref = attenuation_phantoms_Np(f, 3, []);
%att_ref = 0.53*f/8.686; % From phantom especifications

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
        [tempSp1,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd1,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 2
        sub_block_p = ref{2}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{2}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp2,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd2,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 3
        sub_block_p = ref{3}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{3}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp3,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd3,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 4
        sub_block_p = ref{4}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{4}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp4,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd4,~] = spectra(sub_block_d,windowing,0,nw,NFFT);

        tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
        tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));
    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;
save([refDir,'\compensation.mat'],"diffraction_compensation");
