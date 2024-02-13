% ====================================================================== %
% Script to explore homogeneous regions in thyroid tissue. 
% Created on Jan 11, 2024
% ====================================================================== %

clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');


baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Journal'];
targetDir = [baseDir,'\raw'];
figDir = [baseDir,'\fig\23-01-16'];
if (~exist("figDir","dir")), mkdir(figDir); end

targetFiles = dir([targetDir,'\*.mat']);
%disp('Patient list:')
refs     =  [3,2,2,2, 2,3, 2,3, 2,2, 3, 3, 3, 3, 2, 3, 2, 2];
for iAcq = 1:length(targetFiles)
    fprintf("Acquisition no. %i, patient %s, ref %i\n",...
        iAcq,targetFiles(iAcq).name,refs(iAcq));
end
%%

for iAcq = 1:length(targetFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
    patient = targetFiles(iAcq).name(1:end-4);

%% Setting parameters
% Selecting patient
%patient = '135418-07';
ref     = refs(iAcq);
refDir = [baseDir,'\ref',num2str(2)];

% Parameters
blocksize       = 10;     
overlap_pc      = 0.8;
ratio_zx        = 1;
c0              = 1540;
freqC           = 5.5e6;
freq_L = 1e6; freq_H = 15e6;

%% Loading and cropping
load(fullfile(targetDir,patient));
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = RF(:,:,1);
BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));

% Manual cropping
dynRange = [-50,-10];
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.1 3.5])
title(targetFiles(iAcq).name(1:end-4))

confirmation = '';
while ~strcmp(confirmation,'Yes')
    rect = getrect;
    confirmation = questdlg('Sure?');
    if strcmp(confirmation,'Cancel')
        break
    end
end
close,
disp('rect = ')
disp(rect)

%% Calculating parameters
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;

roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = RF(ind_z,ind_x,1);

% Wavelength size
wl = c0/mean(freqC);   % Wavelength (m)

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

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

% Frequency samples
NFFT = 2^(nextpow2(nz)+1);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = attenuation_phantoms_Np(f, 4, []);
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
% windowing = tukeywin(nz/2,0.25);
windowing = hamming(nz/2);
windowing = windowing*ones(1,nx);

% For looping
refFiles = dir([refDir,'\*.mat']);
Nref = length(refFiles);
swrap = saran_wrap(band); % Correction factor for phantom data

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load([refDir,'\',refFiles(iRef).name]);
    samRef = out.RF;
    samRef = samRef(ind_z,ind_x); % Cropping
    % figure,imagesc(db(hilbert(samRef)))
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,swrap,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,swrap,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
            Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
        end
    end
end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref

%% Setting up
% Spectrum
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

% System of eq
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
b = (log(Sp) - log(Sd)) - (compensation);

% Optimization constants
tol = 1e-3;
clear mask
mask = ones(m,n,p);

%% 
NptodB = log10(exp(1))*20;

sldLine = squeeze(mean(mean(b,2),1))/4/L*NptodB;
figure('Units','centimeters', 'Position',[5 5 10 10]),
plot(f,sldLine),
grid on,
xlim([0,freq_H]/1e6),
ylim([0 20]),
xlabel('Frequency [MHz]')
title('Mean SLD')

if iAcq == 1
    sldAllCases = sldLine;
else
    sldAllCases = [sldAllCases,sldLine];
end
end

%% cuec
figure,
% plot(f,median(sldAllCases,2))
errorbar(f,mean(sldAllCases,2),std(sldAllCases,[],2))
% plot(f,sldAllCases, 'k.')
grid on
ylim([0,10])
xlim([0,10])
xlabel('Frequency [MHz]')
ylabel('AC [dB/cm]')
save_all_figures_to_directory(figDir,'thyroidSLD')