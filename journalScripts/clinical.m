% ====================================================================== %
% Script for clinical data. 
% Created on Jan 31, 2024
% ====================================================================== %
clear,clc
close all
    addpath('./functions_v7');
addpath('./AttUtils');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];

tableName = 'clinical.xlsx';

resultsDir = 'C:\Users\sebas\Pictures\Journal2024\24-02-13\BS_8_12';
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

T = readtable('params.xlsx');
%%
blocksize = 8;     % Block size in wavelengths
fixedBW = true;
ratio = db2mag(-30);
freq_L = 3.5e6; freq_H = 8e6;
% freq_L = 3e6; freq_H = 9e6;
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% weights
% muB0 = 1e3; muC0 = 10^0;
% ratioCutOff     = 10;
% order = 5;
% reject = 0.2;
% extension = 3; % 1 or 3

% reg
muBtv = 10^3; muCtv = 10^1;
muBswtv = 10^2.5; muCswtv = 10^0.5;
muBtvl1 = 10^3; muCtvl1 = 10^0;
muBwfr = 10^3; muCwfr = 10^0;

% ================= NEW THINGS  BS 6X12 ==============
% muBtv = 10^3.5; muCtv = 10^1.5;
% muBswtv = 10^3; muCswtv = 10^1;
% muBtvl1 = 10^3.5; muCtvl1 = 10^0.5;
% muBwfr = 10^3.5; muCwfr = 10^0.5;

muB0 = 10^3; muC0 = 10^0;
ratioCutOff     = 10;
order = 5;
reject = 0.2;
extension = 3; % 1 or 3
% =============================================

% swtv weights
aSNR = 1; bSNR = 0.1;
desvMin = 15;

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
        rect = [0.0119    0.2764    1.9230    1.9695]; % 3.5-8 MHz
        % rect = [0.0817    0.2298    1.9850    2.1091]; % 3-9MHz

    case '199031'
        rect = [0.4074    0.9199    2.5200    1.9230];
    case '129424'
        rect = [1.669 0.837 1.625 1.654];
        
    case '189260'
        % rect = [0.923 0.741 1.656 0.929];
        % rect = [0.723 0.741 2.056 1.129];
        % rect = [1.023 0.741 1.756 0.929]; % EL BUENO 31-01
        rect = [0.923 0.741 2.056 0.929]; 


    case '213712'
        rect = [1.683 0.488 1.298 0.9960];
        
    case '265002'
        % rect = [1.6240    0.9431    2.0236    1.4136];
        % rect = [1.6240    1.1431    2.0236    1.2136];
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
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% BW from spectrogram
[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;

if ~fixedBW
    [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
end

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


%% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);

%% SWTV
% Calculating SNR
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% Calculating weights
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% Method
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

%% TV + L1 (no weights)
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);

%% WFR
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
w = (1-reject)* (1./((bscMap/ratioCutOff).^(2*order) + 1)) + reject;
w = movmin(w,extension);

% New equations
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;
% WTV + WL1
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);

%%
figure('Units','centimeters', 'Position',[5 5 30 8]);
tl = tiledlayout(1,3, 'TileSpacing','tight');

t2 = nexttile; 
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
colormap(t2,parula)
axis equal tight
title('TV')
% subtitle(['\mu_b=',num2str(muBtv,2),', \mu_c=',num2str(muCtv,2)])
c = colorbar;
c.Label.String = '\Delta BSC [db/cm]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t2,parula)
axis equal tight
title('SWTV')
% subtitle(['\mu_b=',num2str(muBwfr,2),', \mu_c=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = '[a.u.]';


t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
title('WFR')
% subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';


%%
figure('Units','centimeters', 'Position',[5 5 12 12]);
tl = tiledlayout(2,2);
title(tl,{'Comparison'})
subtitle(tl,{['Patient ',patient],''})

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t2,turbo)
axis equal tight
title('TV')
% subtitle(['\mu_b=',num2str(muBtv,2),', \mu_c=',num2str(muCtv,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t2,turbo)
axis equal tight
title('SWTV')
% subtitle(['\mu_b=',num2str(muBwfr,2),', \mu_c=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';


t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t2,turbo)
axis equal tight
title('TV-L1')
% subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
title('WFR')
% subtitle(['\mu_b=',num2str(muBwfr,2),', \mu_c=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';


%% Mascaras
load(fullfile('newMasks',[patient,'.mat']));
[X,Z] = meshgrid(x,z);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskNoduleACS = interp2(X,Z,maskNodule,Xq,Zq, 'nearest');
maskThyroidACS = interp2(X,Z,maskThyroid,Xq,Zq, 'nearest');

se = strel('diamond',1);
maskThyroidACS = imerode(maskThyroidACS,se);
maskNoduleACS = imerode(maskNoduleACS,se);
%figure, imagesc(x_ACS,z_ACS,maskThyroidACS|maskNoduleACS)

rowId = [rowId;{patient}]; 
col1 = [col1;mean(BRTV(maskNoduleACS))];
col2 = [col2;mean(BRSWTV(maskNoduleACS))];
col3 = [col3;mean(BRTVL1(maskNoduleACS))];
col4 = [col4;mean(BRWFR(maskNoduleACS))];
col5 = [col5;mean(BRTV(maskThyroidACS))];
col6 = [col6;mean(BRSWTV(maskThyroidACS))];
col7 = [col7;mean(BRTVL1(maskThyroidACS))];
col8 = [col8;mean(BRWFR(maskThyroidACS))];

%% Overlay
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure('Units','centimeters', 'Position',[5 5 8 12])
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
contour(x,z,maskThyroid,1,'w--')
hold off
xlabel('x [cm]')
ylabel('z [cm]')

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('TV')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
contour(xFull,zFull,roi,1,'w--')
contour(x,z,maskThyroid,1,'w--')
hold off
colormap(t2,gray)
xlabel('x [cm]')
ylabel('z [cm]')

% nexttile,
% [~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('B-mode and attenuation map')
% hColor.Label.String = 'dB/cm/MHz';
% ylim([0.1, 3.5])
% hold on
% contour(xFull,zFull,roi,1,'w--')
% contour(x,z,maskThyroid,1,'w--')
% hold off


%%
save_all_figures_to_directory(resultsDir,['pat',patient,'fig']);
close all

end

%%
T = table(col1,col2,col3,col4,col5,col6,col7,col8,...
          'VariableNames',{ ...
          'noduleTV','noduleSWTV','noduleTVL1','noduleWFR', ...
          'thyroidTV','thyroidSWTV','thyroidTVL1','thyroidWFR'},...
          'RowNames',rowId)
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
