clear,clc
addpath('./functions_v7');
addpath('./AttUtils');

% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\Simulation_23_12_18'];
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\Timana\batch1'];
baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\process_simulation\23_12_31'];
targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];

if (~exist(croppedDir,"dir")), mkdir(croppedDir); end

targetFiles = dir([targetDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 10;     % Block size in wavelengths
freq_L = 3e6; freq_H = 8.5e6;
%freq_L = 4e6; freq_H = 9e6;
overlap_pc      = 0.8;
ratio_zx        = 1;
referenceAtteanuation = 0.6;

%% For looping
for iAcq = 1:length(targetFiles)
load(fullfile(targetDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);
dynRange = [-50,0];


[pxx,fpxx] = pwelch(sam1,500,400,500,fs);
figure,plot(fpxx/1e6,mean(pxx,2))

% Bmode = db(hilbert(sam1));
% Bmode = Bmode - max(Bmode(:));
% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');

% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
% end
% close,

%% Cropping and finding sample sizes
% Region for attenuation imaging
% x_inf = rect(1); x_sup = rect(1)+rect(3);
% z_inf = rect(2); z_sup = rect(2)+rect(4);
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.5; z_sup = 3.5;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% Lateral samples
wx = round(blocksize*wl*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize*wl/dz /2); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
% dynRange = [-40 -10];
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));
figure, imagesc(x,z,Bmode);
axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = referenceAtteanuation*(f.^1.05)/(20*log10(exp(1)));
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
windowing = windowing*ones(1,nx);

% For looping
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    load(fullfile(refDir,refFiles(iRef).name),"rf");
    samRef = rf;
    samRef = samRef(ind_z,ind_x); % Cropping
    %figure,imagesc(db(hilbert(samRef)))
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
            Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
        end
    end
end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref

%% Spectrum
Sp = zeros(m,n,p);
Sd = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% Saving data
% save(fullfile(croppedDir,targetFiles(iAcq).name),"Sd","Sp",...
%     "compensation","z_ACS","x_ACS","nx","nz","x0","z0p","z0d","sam1",...
%     "m","n","p","Bmode","x","z","f","L", "alpha_coeff_mean")
save(fullfile(croppedDir,targetFiles(iAcq).name),"Sd","Sp",...
    "compensation","z_ACS","x_ACS","nx","nz","x0","z0p","z0d","sam1",...
    "m","n","p","Bmode","x","z","f","L")
end
%%
diffraction_xz = mean(compensation,3);
diffraction_zf = squeeze(mean(compensation,2));
figure, tiledlayout(1,2)
nexttile,
imagesc(x_ACS,z_ACS,diffraction_xz, [-1 1]);
title('Diffraction compensation'),
xlabel('x [cm]'), ylabel('z [cm]'),
colorbar
nexttile,
imagesc(f,z_ACS,diffraction_zf, [-1 1]);
title('Diffraction compensation'),
xlabel('f [MHz]'), ylabel('z [cm]'),
colorbar
