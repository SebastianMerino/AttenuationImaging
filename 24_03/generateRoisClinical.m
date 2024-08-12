% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
clear,clc
close all

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];
resultsDir = fullfile(baseDir,'results','all_cases');
figDir = fullfile(resultsDir,'fig');
T = readtable('params_new.xlsx');
[~,~,~] = mkdir(resultsDir);
[~,~,~] = mkdir(figDir);

maskReady = false;

blocksize = 8;     % Block size in wavelengths
freq_L = 3e6; freq_H = 8.5e6;
freq_C = mean([freq_L freq_H]);
overlap_pc      = 0.8;
ratio_zx        = 12/8;

%%
muB0 = 1e3; muC0 = 10^0;
ratioCutOff     = 10;
order = 5;
reject = 0.3;
extension = 3; % 1 or 3

muBtv = 10^3.5; muCtv = 10^1;
muBwfr = 10^3.5; muCwfr = 10^0;

% Plotting constants
dynRange = [-50,0];
attRange = [0,2];
bsRange = [-20 20];
NptodB = log10(exp(1))*20;

%% Loading case FULL VERSION
for iAcq = 1:height(T)

maskReady = T.mask{iAcq} == 'X';
if maskReady; continue; end

patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});

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

if ~maskReady
    figure('Units','centimeters', 'Position',[3 5 35 15]),
    tiledlayout(1,2)
    nexttile,
    orig = imread(fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.png']));
    image(orig)
    axis image
    nexttile,
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 3.5])
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        h = drawfreehand('Multiclick',true);
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            break
        elseif strcmp(confirmation,'No')
            delete(h)
        end
    end
    regionMask = createMask(h);
else    
    load(fullfile(baseDir,'results','homogeneous',[patient,'.mat']))
    figure('Units','centimeters', 'Position',[3 5 20 15]),
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 min([3.5,max(zFull)])])
end

[X,Z] = meshgrid(xFull,zFull);
x_inf = min(X(regionMask));
x_sup = max(X(regionMask));
z_inf = min(Z(regionMask));
z_sup = max(Z(regionMask));


rectangle('Position', [x_inf, z_inf, x_sup-x_inf, z_sup-z_inf], ...
    'LineStyle','--', 'LineWidth',2, 'EdgeColor','w')


%%
% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
%% 
sam1 = sam1(ind_z,ind_x);
Bmode = BmodeFull(ind_z,ind_x);
Bmode = Bmode - max(Bmode(:));
% regionMask = regionMask(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/freq_C;   % Wavelength (m)

% Lateral samples
wx = round(blocksize*wl*(1-overlap_pc)/dx);  % Between windows
nx = round(blocksize*wl/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz * ratio_zx); % Between windows
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);


% Frequency samples
NFFT = 2^(nextpow2(nz)+1);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);
%%
[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;
[fL,fH] = findFreqBand(fpxx, meanSpectrum, db2mag(-30));
figure,
plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
xline([fL,fH]/1e6)
xlim([1 11])
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')

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
%% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);

% [Xq,Zq] = meshgrid(x_ACS,z_ACS);
% [X,Z] = meshgrid(xFull,zFull);
% regionMaskAcs = interp2(X,Z,double(regionMask),Xq,Zq,"linear");
maskRoi = regionMask(ind_z,ind_x);
regionMaskAcs = maskRoi(z0p,x0) & maskRoi(z0p+nz-1,x0+nx-1) & ...
    maskRoi(z0p,x0+nx-1) & maskRoi(z0p+nz-1,x0);
%% RSLD-TV
maskB = repmat(regionMaskAcs,1,1,p);
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,maskB(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

%% NEW WEIGHTS
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
BTVL1 = reshape(Bn,m,n)*NptodB;

% Weight function
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);
w(~regionMaskAcs) =0;


% figure('Units','centimeters', 'Position',[5 5 30 5]),
% tl = tiledlayout(1,3);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% colormap(t1,gray)
% colorbar
% axis equal
% xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% % axis image
% title('B-mode')
% 
% t2 = nexttile;
% imagesc(x_ACS,z_ACS,bscMap, bsRange)
% colormap(t2,parula)
% c = colorbar;
% ylabel(c,'[dB/cm]')
% axis image
% title('BSC change')
% 
% t3 = nexttile;
% imagesc(x_ACS,z_ACS,BTVL1, attRange)
% colormap(t3,turbo)
% colorbar;
% axis image
% title('ACS')

figure('Units','centimeters', 'Position',[5 5 30 5]),
tiledlayout(1,3);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,bscMap, bsRange)
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


%%
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BRWFR = (reshape(Bn*NptodB,m,n));
CRWFR = (reshape(Cn,m,n));

% figure, 
% imagesc(xFull,zFull,BmodeFull, [-60 -30])
% axis image
% colormap parula
% colorbar

%%
% [X,Z] = meshgrid(x_ACS,z_ACS);
% [Xq,Zq] = meshgrid(xFull,zFull);
% regionMaskAcsInt = interp2(X,Z,regionMaskAcs,Xq,Zq, 'nearest');
alphaImage = 0.7;
figure('Units','centimeters', 'Position',[5 5 25 10]),
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')

nexttile,
[~,~,~] = imoverlay2(BmodeFull,BR,dynRange,...
    attRange,alphaImage,x_ACS,z_ACS,regionMaskAcs,xFull,zFull);
title('TV')
ylim([0.1 min([3.5,max(zFull)])])

nexttile,
[~,~,~] = imoverlay2(BmodeFull,BRWFR,dynRange,...
    attRange,alphaImage,x_ACS,z_ACS,regionMaskAcs,xFull,zFull);
title('WFR')
ylim([0.1 min([3.5,max(zFull)])])

%% SLD fit and Weighted SLD fit
% sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*NptodB/sum(regionMaskAcs(:));
% fit1 = f\sldLine;
% fit2 = [f ones(length(f),1)]\sldLine;
% 
% figure('Units','centimeters', 'Position',[5 5 20 10]),
% tiledlayout(1,2),
% nexttile
% plot(f,sldLine),
% hold on,
% plot(f,fit1*f, '--')
% plot(f,fit2(1)*f + fit2(2), '--')
% hold off,
% grid on,
% xlim([0,freq_H]/1e6),
% ylim([0 16]),
% xlabel('Frequency [MHz]')
% ylabel('Attenuation [dB/cm]')
% title('Mean SLD')
% legend({'SLD','Fit 1', 'Fit 2'}, 'Location','northwest')

sldLine = squeeze(sum(sum(b.*w,2),1))/4/L*NptodB/sum(w(:));
% fit1 = f\sldLine;
% fit2 = [f ones(length(f),1)]\sldLine;
% nexttile,
% plot(f,sldLine),
% hold on,
% plot(f,fit1*f, '--')
% plot(f,fit2(1)*f + fit2(2), '--')
% hold off,
% grid on,
% xlim([0,freq_H]/1e6),
% ylim([0 16]),
% xlabel('Frequency [MHz]')
% ylabel('Attenuation [dB/cm]')
% title('Weighted mean SLD')
% legend({'SLD','Fit 1', 'Fit 2'}, 'Location','northwest')

%%
save_all_figures_to_directory(figDir,['pat',patient,'fig']);
close all
save(fullfile(resultsDir,[patient,'.mat']),...
    "sldLine",'regionMask','BR','BRWFR', 'regionMaskAcs', 'w')
end

