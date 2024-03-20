% ====================================================================== %
% Script to explore homogeneous ROIs in clinical data. 
% Created on Mar 5, 2024
% ====================================================================== %
clear,clc
close all

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];

T = readtable('params.xlsx');

resultsDir = 'C:\Users\sebas\Pictures\Journal2024\24-03-05\';
if (~exist("figDir","dir")), mkdir(resultsDir); end


blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
% freq_L = 3e6; freq_H = 9e6;
freq_C = 5.75e6;
overlap_pc      = 0.8;
ratio_zx        = 12/8;

muB0 = 1e3; muC0 = 10^0;
ratioCutOff     = 15;
order = 5;
reject = 0.2;
extension = 3; % 1 or 3

muBtv = 10^3.5; muCtv = 10^1;
muBwfr = 10^3.5; muCwfr = 10^0;

% Plotting constants
dynRange = [-50,0];
attRange = [0.4,1.9];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;
%% Loading case FULL VERSION
% for iAcq = 1:height(T)
% 
% iAcq = 7;
% rect1 = [2.7, 1, 1.1, 1.3];
% rect2 = [1.6, 1.5, 1.1, 0.8];
% rect1 = [1.6, 0.7, 2.5,2];

% iAcq = 6;
% rect2 = [1.4, 0.4, 0.7, 0.8];
% rect1 = [1.4, 1.3, 0.7, 0.8];
% 
% iAcq = 4;
% rect1 = [0.7    1.1    0.9    1.2];
% rect2 = [2.1   1.1    0.9    1.2];
% ratioCutOff     = 10;

iAcq = 1;
rect1 = [0.4    0.75    0.4    0.8];
rect2 = [1.6    0.75    0.4    0.8];

patient = num2str(T.patient(iAcq));
samPath = fullfile(baseDir,patient,[patient,'-',T.sample{iAcq},'.rf']);
refDir = fullfile(refsDir,T.reference{iAcq});


out =lectura_OK(samPath);
sam1 = out.RF(:,:,1);
fs = out.fs;
fc = out.fc;
x = out.x; z = out.z;

dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
colormap gray; % clim([-70 -20]);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.1 3.5])
rectangle('Position',rect1, 'EdgeColor','w', 'LineStyle','--')
rectangle('Position',rect2, 'EdgeColor','w', 'LineStyle','--')


%% ===================================================================== %

for iRoi = 1:2
out =lectura_OK(samPath);
sam1 = out.RF(:,:,1);
fs = out.fs;
fc = out.fc;
x = out.x; z = out.z;

dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]


BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));

if iRoi == 1
    rect = rect1;
else
    rect = rect2;
end
%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = sam1(ind_z,ind_x);
Bmode = BmodeFull(ind_z,ind_x);

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
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
% Bmode = db(hilbert(sam1));
% Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%%
%[pxx,fpxx] = pwelch(sam1,nz,nz-wz,nz,fs);
%plot(fpxx,db(meanSpectrum))

%%
pxx = fft(sam1);
fpxx = (0:size(pxx,1)-1)/size(pxx,1)*fs;
meanSpectrum = mean(pxx,2);
spec{iRoi} =  db(meanSpectrum);

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

b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);

%% SLD fit
sldLine = squeeze(mean(mean(b,2),1))/4/L*NptodB;
fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;
figure('Units','centimeters', 'Position',[5 5 10 10]),
plot(f,sldLine),
hold on,
plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), '--')
hold off,
grid on,
xlim([0,freq_H]/1e6),
ylim([0 16]),
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Mean SLD')
legend({'SLD','Fit 1', 'Fit 2'}, 'Location','northwest')

fprintf("\nACS when forzing zero cross: %f \n",fit1);
fprintf("ACS when not forzing zero cross: %f \n\n",fit2(1));
%% RSLD-TV
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

% figure('Units','centimeters', 'Position',[5 5 10 8]);
% iRoi = 1;
% imageData.x = x_ACS;
% imageData.z = z_ACS;
% imageData.roi = roi;
% imageData.TV = BR;
% % imageData.WFR = BRWTik;
% dataRoi{iRoi} = imageData;
% [~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.TV,dynRange,...
%     attRange,0.7,...
%     dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);

%% NEW WEIGHTS
tol = 1e-3;
mask = ones(m,n,p);
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

% Weight function
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

figure('Units','centimeters', 'Position',[5 5 15 5]),
tl = tiledlayout(1,3);
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
% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
% W = speye(m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*NptodB,m,n));
CRWTik = (reshape(Cn,m,n));

imageData.x = x_ACS;
imageData.z = z_ACS;
imageData.roi = roi;
imageData.TV = BR;
imageData.WFR = BRWTik;
dataRoi{iRoi} = imageData;

%% Weighted SLD fit
sldLine = squeeze(sum(sum(b.*w,2),1))/4/L*NptodB/sum(w(:));
fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;
figure('Units','centimeters', 'Position',[5 5 10 10]),
plot(f,sldLine),
hold on,
plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), '--')
hold off,
grid on,
xlim([0,freq_H]/1e6),
ylim([0 16]),
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Weighted mean SLD')
legend({'SLD','Fit 1', 'Fit 2'}, 'Location','northwest')

fprintf("\nACS when forzing zero cross: %f \n",fit1);
fprintf("ACS when not forzing zero cross: %f \n\n",fit2(1));

end

%% Mean spectra
figure('Units', 'centimeters', 'Position',[5 5 10 10]),
plot(fpxx/1e6,spec{1})
hold on
plot(fpxx/1e6,spec{2})
hold off
grid on
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
xline(freq_L/1e6,'k--')
xline(freq_H/1e6,'k--')
legend('Thyroid', 'Nodule')
xlim([0 20])
%% ACS images

alpha = 0.7;
% dynRange = [-40 -10];

figure('Units','centimeters', 'Position',[5 5 10 8]);
imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('B-mode')
fontsize(gcf,8,'points')


figure('Units','centimeters', 'Position',[5 5 10 8]);
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.TV,dynRange,...
    attRange,alpha,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.TV,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
iRoi = 2;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('RSLD-TV')
hColor.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')


figure('Units','centimeters', 'Position',[5 5 10 8]);
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.WFR,dynRange,...
    attRange,alpha,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.WFR,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
iRoi = 2;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('RSLD-WFR')
hColor.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

%%
fprintf("Thyroid: \n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{1}.TV(:)),...
    std(dataRoi{1}.TV(:)))
fprintf("Nodule: \n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{2}.TV(:)),...
    std(dataRoi{2}.TV(:)))

fprintf("Thyroid: \n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{1}.WFR(:)),...
    std(dataRoi{1}.WFR(:)))
fprintf("Nodule: \n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{2}.WFR(:)),...
    std(dataRoi{2}.WFR(:)))

%%

% save_all_figures_to_directory(resultsDir,['pat',patient,'fig']);
% close all

