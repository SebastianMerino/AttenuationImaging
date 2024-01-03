% ====================================================================== %
% Script to explore homogeneous ROIs in clinical data. 
% Created on Dec 26, 2023
% ====================================================================== %
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\carcinoma'];


targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
figDir = [baseDir,'\fig\12-26'];
if (~exist("figDir","dir")), mkdir(figDir); end

targetFiles = dir([targetDir,'\*.mat']);
disp('Patient list:')
for iAcq = 1:length(targetFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
end
blocksize = 10;     % Block size in wavelengths
overlap_pc      = 0.8;
ratio_zx        = 1;

%% Loading case FULL VERSION
iAcq = 1;
load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("\n Selecting acq. no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
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
c0 = 1540; freqC = 5.5e6;
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

%% Choosing bandwidth
ratio = db2mag(-30);

% BW from spectrogram
[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
[freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
freq_L = 2e6; freq_H = 10e6;
% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plotting BW
figure,plot(fpxx/1e6,meanSpectrum)
yline(max(meanSpectrum)*ratio)
xline([freq_L,freq_H]/1e6)
xlabel('Frequency [MHz]')
ylabel('Magnitude')

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = attenuation_phantoms_Np(f, 3, []);
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
attRange = [0.3,1.7];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

figure,
plot(f,squeeze(mean(mean(b,2),1))'),
grid on,
xlim([0,freq_H]/1e6),ylim([0 1]),
xlabel('Frequency [MHz]')
title('Mean SLD')
%% RSLD

%muB = 10.^(2.5:0.5:3.5);
muC = 10.^(0:1:2)*10;
muB = 10^4;

for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        disp(mean(BR(:)))
        CR = reshape(Cn*NptodB,m,n);
        % figure('Units','centimeters', 'Position',[5 5 20 8]);
        % tl = tiledlayout(1,3);
        figure('Units','centimeters', 'Position',[5 5 15 15]);
        tl = tiledlayout(3,1);

        title(tl,'RSLD-TV')
        subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
        t1 = nexttile;
        imagesc(x,z,Bmode,dynRange)
        axis equal
        xlim([x_ACS(1) x_ACS(end)]),
        ylim([z_ACS(1) z_ACS(end)]),
        colormap(t1,gray)
        colorbar(t1,'westoutside')
        title('Bmode')

        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR, attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD, \mu=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';

        t3 = nexttile; 
        imagesc(x_ACS,z_ACS,CR, bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD, \mu=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio (a.u.)';
    end
end


%% Minimizing BS log ratio

%muB = 10.^(2.5:0.5:3.5)*10;
muB = 10^4;
muC = 10.^(0:0.5:1)*10;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);
        disp(mean(BR(:)))
        %figure('Units','centimeters', 'Position',[5 5 20 8]);
        %tl = tiledlayout(1,3);
        figure('Units','centimeters', 'Position',[5 5 15 15]);
        tl = tiledlayout(3,1);
        title(tl,'TV-L1')
        subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
        t1 = nexttile;
        imagesc(x,z,Bmode,dynRange)
        axis equal
        xlim([x_ACS(1) x_ACS(end)]),
        ylim([z_ACS(1) z_ACS(end)]),
        colormap(t1,gray)
        colorbar(t1,'westoutside')
        title('Bmode')

        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR, attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD, \mu=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';

        t3 = nexttile; 
        imagesc(x_ACS,z_ACS,CR, bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD, \mu=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio (a.u.)';
    end
end


%% NEW WEIGHTS
% First estimation
muB0 = 1e3; muC0 = 10^0;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

% Weight function
ratioCutOff = 20;
order = 5;
reject = 0.1;
extension = 3; % 1 or 3
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% New equations
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

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
ylabel(c,'dB')
axis image
title('BS ratio')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')

%% Weighting equation and regularizations
muB = 10.^(3:0.5:4);
muC = 10.^(0:0.5:1);

for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,...
          muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);
        figure('Units','centimeters', 'Position',[5 5 30 8]);
        tl = tiledlayout(1,3);
        title(tl,'WFR')
        subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
        t1 = nexttile;
        imagesc(x,z,Bmode,dynRange)
        axis equal
        xlim([x_ACS(1) x_ACS(end)]),
        ylim([z_ACS(1) z_ACS(end)]),
        colormap(t1,gray)
        colorbar(t1,'westoutside')
        title('Bmode')

        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR, attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD, \mu=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';

        t3 = nexttile; 
        imagesc(x_ACS,z_ACS,CR, bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD, \mu=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio (a.u.)';
    end
end

% ======================================================================
% ======================================================================
% ======================================================================
%% ACUMULADO

muBtv = 10^4; muCtv = 10^3;
muBtvl1 = 10^4; muCtvl1 = 10^2;

% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);

% British Columbia
% [Bn,~] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
% m,n,tol,mask(:),wBC);
% BRSWTV = reshape(Bn*NptodB,m,n);

% TV + L1 (no weights)
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);

% WFR   
% [Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
% BRWFR = reshape(Bn*NptodB,m,n);

%%
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
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRTV,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
colormap(t2,gray)


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
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRTVL1,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
colormap(t2,gray)


%%
save_all_figures_to_directory(figDir,[targetFiles(iAcq).name(1:end-4),'fig']);
close all