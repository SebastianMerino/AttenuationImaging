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
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];

T = readtable('params.xlsx');

figDir = [baseDir,'\fig\24-01-25\uniqueTightBW-W2-1'];
if (~exist("figDir","dir")), mkdir(figDir); end


blocksize = 12;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 1;

weightMap = 1;
ratioCutOff     = 10;
order = 5;
reject = 0.2;
extension = 3; % 1 or 3

% weightMap = 2;
% dBgain = 0.5;

muBtv = 10^3; muCtv = 10^1;
muBtvl1 = 10^3; muCtvl1 = 10^0;
muBwfr2 = 10^3; muCwfr2 = 10^1;
muBwfr = 10^3; muCwfr = 10^0;

%%
rowId = []; 
col1 = [];
col2 = [];
col3 = [];
col4 = [];
col5 = [];
col6 = [];
col7 = [];
col8 = [];
for iAcq = 1:height(T)
patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});

rect = [];
switch patient
    case '135418'
        %rect = [1.5162    0.3772    2.2564    1.5740];
        rect = [];
        % 1.6170    0.4159    2.2796    1.5275

    case '190029'
        % rect = [1.3689    0.5865    1.1088    1.4267];
        % rect = [0.8183    0.4314    2.3804    1.7679];
        rect = [1.0183    0.4314    1.9804    1.7679];

    case '203704'
        % rect = [1.1518    0.4857    2.6131    1.9927];
        rect = [1.1518    0.4857  2 2];

    case '254581'
        rect = [1.03; 0.49; 1.6; 1.69];

    case '134135'
        % rect = [0.0119    0.2764    1.9230    1.9695]; % 3.5-8 MHz
        % rect = [0.0817    0.2298    1.9850    2.1091]; % 3-9MHz
        rect = [0.0817    0.2298    1.9850    1.9091]; 

    case '199031'
        rect = [0.4074    0.9199    2.5200    1.9230];
    case '129424'
        rect = [1.669 0.837 1.625 1.654];
        
    case '189260'
        % rect = [0.923 0.741 1.656 0.929]; % BUENO
        % rect = [0.923 0.741 1.656 1.029];
        % rect = [0.723 0.741 2.056 1.129];
        rect = [1.023 0.741 1.756 0.929];
        
    case '213712'
        rect = [1.683 0.488 1.298 0.9960];
        
    case '265002'
        %rect = [1.6240    0.9431    2.0236    1.4136];
        rect = [1.6240    1.0431    2.0236    1.3136];

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
            disp(rect)
            break
        end
    end
    close,
end

%%
x_inf = rect(1); x_sup = rect(1)+rect(3); 
z_inf = rect(2); z_sup = rect(2)+rect(4);

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

%% BW from spectrogram
ratio = db2mag(-30);

[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;
[freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
% freq_L = 3e6; freq_H = 9e6;
% 
% Plotting BW
figure,plot(fpxx/1e6,meanSpectrum)
yline(max(meanSpectrum)*ratio)
xline([freq_L,freq_H]/1e6)
xlabel('Frequency [MHz]')
ylabel('Magnitude')


NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('Frequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
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

% Plotting constants
dynRange = [-40,-5];
attRange = [0.4,1.9];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% NEW WEIGHTS
% First estimation
muB0 = 1e3; muC0 = 10^0;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

if weightMap == 1
    w = (1-reject)* (1./((bscMap/ratioCutOff).^(2*order) + 1)) + reject;
    w = movmin(w,extension);
else
    w = db2mag(-dBgain*abs(bscMap));
end

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

%%
% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);


% TV + L1 (no weights)
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);

% New equations
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

figure('Units','centimeters', 'Position',[5 5 20 15]);
tl = tiledlayout(2,2);
title(tl,{'Comparison'})
subtitle(tl,{['Patient ',patient],''})

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t2,turbo)
axis equal tight
title('TV')
subtitle(['\mu_b=',num2str(muBtv,2),', \mu_c=',num2str(muCtv,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
% rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
% rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t2,turbo)
axis equal tight
title('TV-L1')
subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
% rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
% rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR2, attRange)
colormap(t2,turbo)
axis equal tight
title('WTV + WTV')
subtitle(['\mu_b=',num2str(muBwfr2,2),', \mu_c=',num2str(muCwfr2,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
% rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
% rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
title('WTV + WL1')
subtitle(['\mu_b=',num2str(muBwfr,2),', \mu_c=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
% rectangle('Position',[x0Tumor z0Tumor wTumor hTumor], 'LineStyle','--', 'EdgeColor','w')
% rectangle('Position',[x0Sano z0Sano wSano hSano], 'LineStyle','--', 'EdgeColor','k')

%% Mascaras
%[X,Z] = meshgrid(x_ACS,z_ACS);
% maskTumor = X>x0Tumor & X<x0Tumor+wTumor & Z>z0Tumor & Z<z0Tumor + hTumor;
%maskSano = X>x0Sano & X<x0Sano+wSano & Z>z0Sano & Z<z0Sano + hSano;
load(fullfile('newMasks',[patient,'.mat']));
[X,Z] = meshgrid(x,z);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskNoduleACS = interp2(X,Z,maskNodule,Xq,Zq, 'nearest');
% figure,
% imagesc(x_ACS,z_ACS,maskNoduleACS)
% axis image

maskThyroidACS = interp2(X,Z,maskThyroid,Xq,Zq, 'nearest');
rowId = [rowId;{patient}]; 
col1 = [col1;mean(BRTV(maskNoduleACS))];
col2 = [col2;mean(BRTVL1(maskNoduleACS))];
col3 = [col3;mean(BRWFR2(maskNoduleACS))];
col4 = [col4;mean(BRWFR(maskNoduleACS))];
col5 = [col5;mean(BRTV(maskThyroidACS))];
col6 = [col6;mean(BRTVL1(maskThyroidACS))];
col7 = [col7;mean(BRWFR2(maskThyroidACS))];
col8 = [col8;mean(BRWFR(maskThyroidACS))];

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
% contour(xFull,zFull,roi,1,'w--')
% contour(x,z,maskThyroid,1,'w--')
hold off

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
% contour(xFull,zFull,roi,1,'w--')
% contour(x,z,maskThyroid,1,'w--')
hold off

colormap(t2,gray)

%%
save_all_figures_to_directory(figDir,['pat',patient,'fig']);
close all

end

%%
T = table(col1,col2,col3,col4,col5,col6,col7,col8,...
          'VariableNames',{ ...
          'noduleTV','noduleTVL1','noduleWTV','noduleWTVL1', ...
          'thyroidTV','thyroidTVL1','thyroidWTV','thyroidWTVL1'},...
          'RowNames',rowId)
writetable(T,'.\24-01-25\uniqueTightBW-W2-1.xlsx'...
    , 'WriteRowNames',true);
