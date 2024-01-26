% ====================================================================== %
% Script to create masks for heterogenenous ROIs in clinical data. 
% Created on Jan 11, 2024
% ====================================================================== %
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];

T = readtable('params.xlsx');

figDir = [baseDir,'\fig\24-01-23\uniqueBW-W1-4'];
if (~exist("figDir","dir")), mkdir(figDir); end


blocksize = 12;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 1;

%% Loading case FULL VERSION
for iAcq = 2:height(T)
patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});

rect = [];
switch patient
    case '135418'
        %rect = [1.5162    0.3772    2.2564    1.5740];
        rect = [];
        % 1.6170    0.4159    2.2796    1.5275
        x0Tumor = 1.9; z0Tumor = 1;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 3; z0Sano = 1;
        wSano = 0.5; hSano = 0.5; 

    case '190029'
        % rect = [1.3689    0.5865    1.1088    1.4267];
        % rect = [0.8183    0.4314    2.3804    1.7679];
        rect = [1.0183    0.4314    1.9804    1.7679];
        x0Tumor = 1.7; z0Tumor = 0.9;
        wTumor = 0.5; hTumor = 0.4;
        x0Sano = 1.7; z0Sano = 1.6;
        wSano = 0.5; hSano = 0.2; 

    case '203704'
        % rect = [1.1518    0.4857    2.6131    1.9927];
        rect = [1.1518    0.4857  2 2];
        x0Tumor = 1.8; z0Tumor = 1;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 1.8; z0Sano = 1.7;
        wSano = 0.5; hSano = 0.5; 

    case '254581'
        rect = [1.03; 0.49; 1.6; 1.69];
        x0Tumor = 1.4; z0Tumor = 0.7;
        wTumor = 0.7; hTumor = 0.5;
        x0Sano = 1.4; z0Sano = 1.4;
        wSano = 0.7; hSano = 0.5; 

    case '134135'
        rect = [0.0119    0.2764    1.9230    1.9695]; % 3.5-8 MHz
        % rect = [0.0817    0.2298    1.9850    2.1091]; % 3-9MHz
        x0Tumor = 1.3; z0Tumor = 0.8;
        wTumor = 0.4; hTumor = 0.6;
        x0Sano = 0.4; z0Sano = 0.8;
        wSano = 0.4; hSano = 0.6; 

    case '199031'
        rect = [0.4074    0.9199    2.5200    1.9230];
        x0Tumor = 2.1; z0Tumor = 1.6;
        wTumor = 0.6; hTumor = 0.6;
        x0Sano = 0.6; z0Sano = 1.6;
        wSano = 0.6; hSano = 0.6; 
    case '129424'
        rect = [1.669 0.837 1.625 1.654];
        x0Tumor = 1.8; z0Tumor = 1.4;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 2.6; z0Sano = 1.4;
        wSano = 0.5; hSano = 0.5; 
        
    case '189260'
        rect = [0.923 0.741 1.656 0.929];
        % rect = [0.723 0.741 2.056 1.129];
        x0Tumor = 1.8; z0Tumor = 0.9;
        wTumor = 0.5; hTumor = 0.4;
        x0Sano = 1.1; z0Sano = 0.9;
        wSano = 0.5; hSano = 0.4; 
        
    case '213712'
        rect = [1.683 0.488 1.298 0.9960];
        x0Tumor = 1.8; z0Tumor = 0.6;
        wTumor = 0.45; hTumor = 0.35;
        x0Sano = 2.4; z0Sano = 1;
        wSano = 0.45; hSano = 0.35; 
        
    case '265002'
        rect = [1.6240    0.9431    2.0236    1.4136];
        x0Tumor = 1.9; z0Tumor = 1.6;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 2.8; z0Sano = 1.6;
        wSano = 0.5; hSano = 0.5; 

    case '266844'
        % rect = [0.3531 0.6098 2.6286 2.1788];
        rect = [0.3531 0.8098 2.6286 1.9788];

    case '20192'
        rect = [0.5082 0.6951 3.0550 1.7989];
    otherwise
        rect = [];
end
%%


out =lectura_OK(samPath);
sam1 = out.RF(:,:,1);
fs = out.fs;
fc = out.fc;
x = out.x; z = out.z;
fprintf("\n Selecting acq. no. %i, patient %s\n",iAcq,patient);


% Manual cropping
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));


if isempty(rect)
    % Manual cropping
    dynRange = [-50,-10];
    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 3.5])
    title(patient)
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        rect = getrect;
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            break
        end
    end
    close,
end

%%
% Carcinoma, caso 221778-01
% rect = [0.0042    0.7270    3.1558    2.2486];
% rect = [0.0119    0.2764    1.9230    1.9695];
x_inf = rect(1); x_sup = rect(1)+rect(3); 
z_inf = rect(2); z_sup = rect(2)+rect(4);
disp(rect)

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
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

%%
% Manual cropping
dynRange = [-40,0];
mask = zeros(size(Bmode));
[X,Z] = meshgrid(x,z);

figure('Units','centimeters', 'Position',[5 5 30 15]),

tiledlayout(1,2)
t1 = nexttile;
im1 = imagesc(x,z,Bmode,dynRange); axis image;
colormap(t1,gray);
clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
title(patient)

t2 = nexttile;
imagesc(t2, x,z, mask, [0 1])
axis image
colormap(t2, 'parula')

while true
    rect = getrect(t1);
    rw = rect(3)/2;
    rh = rect(4)/2;
    xc = rect(1) + rw;
    zc = rect(2) + rh;
    mask2 = mask  |( ((X-xc)/rw).^2 +((Z-zc)/rh).^2 < 1 );  
    im2 = imagesc(t2, x,z, mask2, [0 1]);
    axis(t2,'image')
    im1.AlphaData = ~mask2 + 0.8;
    confirmation = questdlg('Is it ok?');
    if strcmp(confirmation,'Yes')
        mask = mask2;
    elseif strcmp(confirmation,'No')
        imagesc(t2, x,z, mask, [0 1]);
        im1.AlphaData = ~mask + 0.8;
        axis(t2,'image')
    elseif strcmp(confirmation,'Cancel')
        mask = mask2;
        break;
    end
end
close,

maskNodule = mask;

%%

mask = ~maskNodule;
[X,Z] = meshgrid(x,z);

figure('Units','centimeters', 'Position',[5 5 30 15]),

tiledlayout(1,2)
t1 = nexttile;
im1 = imagesc(x,z,Bmode,dynRange); axis image;
colormap(t1,gray);
clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
title(patient)

t2 = nexttile;
imagesc(t2, x,z, mask, [0 1])
axis image
colormap(t2, 'parula')

while true
    rect = getrect(t1);
    rw = rect(3)/2;
    rh = rect(4)/2;
    xc = rect(1) + rw;
    zc = rect(2) + rh;
    mask2 = mask  & ~( ((X-xc)/rw).^2 +((Z-zc)/rh).^2 < 1 );  
    im2 = imagesc(t2, x,z, mask2, [0 1]);
    axis(t2,'image')
    im1.AlphaData = ~mask2 + 0.8;
    confirmation = questdlg('Is it ok?');
    if strcmp(confirmation,'Yes')
        mask = mask2;
    elseif strcmp(confirmation,'No')
        imagesc(t2, x,z, mask, [0 1]);
        im1.AlphaData = ~mask + 0.8;
        axis(t2,'image')
    elseif strcmp(confirmation,'Cancel')
        mask = mask2;
        break;
    end
end
close,

% [Xq,Zq] = meshgrid(x_ACS,z_ACS);
% maskNodule = interp2(X,Z,mask,Xq,Zq, 'nearest');
% figure,
% imagesc(x_ACS,z_ACS,maskNodule)
% axis image
% 
maskThyroid = mask;
%%
save(fullfile('newMasks',patient), 'maskNodule', 'maskThyroid','-mat');
end