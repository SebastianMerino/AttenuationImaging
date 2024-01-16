% ====================================================================== %
% Script to explore heterogenenous ROIs in clinical data. 
% Created on Jan 11, 2024
% ====================================================================== %
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Journal'];

targetDir = [baseDir,'\raw'];
targetFiles = dir(fullfile(targetDir,'*.mat'));
figDir = [baseDir,'\fig\23-01-16'];
if (~exist("figDir","dir")), mkdir(figDir); end
refs     =  [3,2,2,2, 2,3, 2,3, 2,2, 3, 3, 3, 3, 2,3,2,2];
for iAcq = 1:length(targetFiles)
    fprintf("Acquisition no. %i, patient %s, ref %i\n",...
        iAcq,targetFiles(iAcq).name,refs(iAcq));
end

%%
% Selecting patient
iAcq = 3;
patient = targetFiles(iAcq).name(1:end-4);
ref     = refs(iAcq);
refDir = [baseDir,'\ref',num2str(2)];

% Parameters
blocksize       = 10;     
overlap_pc      = 0.8;
ratio_zx        = 1;
c0              = 1540;
freqC           = 5.5e6;
%spectrumCut     = db2mag(-30);
% freq_L = 3.5e6; freq_H = 8e6;
freq_L = 3e6; freq_H = 9e6;

%%
rect = [];
switch patient
    case '135418-07'
        rect = [1.5162    0.3772    2.2564    1.5740];
        ratioCutOff = 10;
        x0Tumor = 1.9; z0Tumor = 1;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 3; z0Sano = 1;
        wSano = 0.5; hSano = 0.5; 

    case '190029'
        rect = [1.3689    0.5865    1.1088    1.4267];
        ratioCutOff = 10;
        x0Tumor = 1.7; z0Tumor = 0.9;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 1.7; z0Sano = 1.6;
        wSano = 0.5; hSano = 0.2; 

    case '203704-03'
        rect = [1.1518    0.4857    2.6131    1.9927];
        ratioCutOff = 10;
        x0Tumor = 1.8; z0Tumor = 1;
        wTumor = 0.5; hTumor = 0.5;
        x0Sano = 1.8; z0Sano = 1.7;
        wSano = 0.5; hSano = 0.5; 

    case '254581-1'
        rect = [1.03; 0.49; 1.6; 1.69];
        x0Tumor = 1.4; z0Tumor = 0.7;
        wTumor = 0.7; hTumor = 0.5;
        x0Sano = 1.4; z0Sano = 1.4;
        wSano = 0.7; hSano = 0.5; 

    case '134135'
        rect = [0.0119    0.2764    1.9230    1.9695]; % 3.5-8 MHz
        % rect = [0.0817    0.2298    1.9850    2.1091]; % 3-9MHz
        ratioCutOff = 20;
        x0Tumor = 1.3; z0Tumor = 0.8;
        wTumor = 0.4; hTumor = 0.6;
        x0Sano = 0.4; z0Sano = 0.8;
        wSano = 0.4; hSano = 0.6; 

    case '199031'
        rect = [0.4074    0.9199    2.5200    1.9230];
        ratioCutOff = 20;
        x0Tumor = 2.1; z0Tumor = 1.6;
        wTumor = 0.6; hTumor = 0.6;
        x0Sano = 0.5; z0Sano = 1.6;
        wSano = 0.6; hSano = 0.6; 
    otherwise
        rect = [];
        ratioCutOff = 20;
        x0Tumor = 2.1; z0Tumor = 1.6;
        wTumor = 0.6; hTumor = 0.6;
        x0Sano = 0.5; z0Sano = 1.6;
        wSano = 0.6; hSano = 0.6; 
end

%% Loading and cropping
load(fullfile(targetDir,patient));
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = RF(:,:,1);
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
end

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
windowing = tukeywin(nz/2,0.25);
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

% Plotting constants
dynRange = [-40,-5];
attRange = [0.5,1.9];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% NEW WEIGHTS
% First estimation
muB0 = 1e3; muC0 = 10^0;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

% Weight function
order = 5;
% gain = 1;
reject = 0.1;
extension = 3; % 1 or 3
w = (1-reject)* (1./((bscMap/20).^(2*order) + 1)) + reject;
w = movmin(w,extension);
%w = db2mag(-abs(bscMap)*0.5);
%%
figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'New Weights',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
colormap(t2,parula)
c = colorbar;
ylabel(c,'[dB/cm]')
axis image
title('BSC change')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
colorbar;
axis image
title('Weights')

% %% Modifying weights

% w(10:14,15:17) = 0.9;
% w(10:14,16) = 0.1;
% w(15,17) = 0.1;
% w(9,17) = 0.1;
% 
% t3 = nexttile;
% imagesc(x_ACS,z_ACS,w, [0 1])
% colormap(t3,parula)
% colorbar;
% axis image
% title('Weights')
% 
% 
% 
% [X,Z] = meshgrid(x_ACS,z_ACS);
% w(X<1) = 1;
% w(Z>1.8) = 1;
% figure('Units','centimeters', 'Position',[5 5 30 8]),
% tl = tiledlayout(1,2);
% 
% t3 = nexttile;
% imagesc(x_ACS,z_ACS,w, [0 1])
% colormap(t3,parula)
% colorbar;
% axis image
% title('Weights')
% 
% % t3 = nexttile;
% % imagesc(x_ACS,z_ACS,w2, [0 1])
% % colormap(t3,parula)
% % colorbar;
% % axis image
% % title('Modified weights')
% 
%%
muBtv = 10^3.5; muCtv = 10^1.5;
muBtvl1 = 10^3.5; muCtvl1 = 10^0.5;
muBwfr2 = 10^3.5; muCwfr2 = 10^1.5;
muBwfr = 10^3.5; muCwfr = 10^0.5;
% muBwfr2 = 10^3; muCwfr2 = 10^1;
% muBwfr = 10^3; muCwfr = 10^0;

% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);


% TV + L1 (no weights)
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);

% New equations
% W = repmat(movmin(w,extension),[1 1 p]);
% w = w2;
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% WTV + WTV
[Bn,~] = optimAdmmWeightedTv(A1w,A2w,bw,muBwfr2,muCwfr2,m,n,tol,mask(:),w);
BRWFR2 = reshape(Bn*NptodB,m,n);

% WTV + WL1
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);

%%

figure('Units','centimeters', 'Position',[5 5 20 15]);
tl = tiledlayout(2,2);
title(tl,{'Comparison'})
subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t2,turbo)
axis equal tight
title('TV')
subtitle(['\mu_b=',num2str(muBtv,2),', \mu_c=',num2str(muCtv,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t2,turbo)
axis equal tight
title('TV-L1')
subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR2, attRange)
colormap(t2,turbo)
axis equal tight
title('WTV + WTV')
subtitle(['\mu_b=',num2str(muBwfr2,2),', \mu_c=',num2str(muCwfr2,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
title('WTV + WL1')
subtitle(['\mu_b=',num2str(muBwfr,2),', \mu_c=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

%% Mascaras
[X,Z] = meshgrid(x_ACS,z_ACS);
maskTumor = X>x0Tumor & X<x0Tumor+wTumor & Z>z0Tumor & Z<z0Tumor + hTumor;
maskSano = X>x0Sano & X<x0Sano+wSano & Z>z0Sano & Z<z0Sano + hSano;

fileID = fopen('casos.txt','a');
fprintf(fileID,"\nPatient %s \n",patient);
fprintf(fileID, "TV + TV\n");
fprintf(fileID, "Nodulo:\n %.2f\n", mean(BRTV(maskTumor)));
fprintf(fileID, "Tiroides:\n %.2f\n", mean(BRTV(maskSano)));

fprintf(fileID, "TV + L1\n");
fprintf(fileID, "Nodulo:\n %.2f\n", mean(BRTVL1(maskTumor)));
fprintf(fileID, "Tiroides:\n %.2f\n", mean(BRTVL1(maskSano)));

fprintf(fileID, "WTV + WTV\n");
fprintf(fileID, "Nodulo:\n %.2f\n", mean(BRWFR2(maskTumor)));
fprintf(fileID, "Tiroides:\n %.2f\n", mean(BRWFR2(maskSano)));

fprintf(fileID, "WTV + WL1\n");
fprintf(fileID, "Nodulo:\n %.2f\n", mean(BRWFR(maskTumor)));
fprintf(fileID, "Tiroides:\n %.2f\n", mean(BRWFR(maskSano)));
fclose(fileID);

%% Overlay
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure('Units','centimeters', 'Position',[5 5 10 15])
tiledlayout(2,1)
t2 = nexttile;
imagesc(xFull,zFull,BmodeFull,[-50 0])
axis image
title('B-mode')
ylim([0.1, 3.5])
c = colorbar;
c.Label.String = 'dB';
hold on
contour(xFull,zFull,roi,1,'w--')
hold off

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
contour(xFull,zFull,roi,1,'w--')
hold off

colormap(t2,gray)

%%
save_all_figures_to_directory(figDir,['pat',patient,'fig']);
close all
