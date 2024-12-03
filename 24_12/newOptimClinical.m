% ====================================================================== %
% Script for clinical data. 
% Created on Jan 31, 2024
% ====================================================================== %
setup,
baseDir = 'C:\Users\sebas\Documents\Data\Attenuation\Thyroid_Data_PUCP_UTD';
refsDir = 'C:\Users\sebas\Documents\Data\Attenuation\REFERENCES';
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\UFFC2024Results\newOptim';

tableName = 'clinical.xlsx';
T = readtable('params.xlsx');
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

%%
blocksize = 8;     % Block size in wavelengths
fixedBW = true;
ratio = db2mag(-30);
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% Weight parameters
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% reg FINAL VERSION
muBtv = 10^2.5; muCtv = 10^2.5;
muBswtv = 10^2.5; muCswtv = 10^-0.5;
muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
muBswift = 10^4; muCswift = 10^1.5;

% Plotting constants
dynRange = [-50,0];
attRange = [0.2,2];
bsRange = [-15 15];

dataCols = zeros(3,12);

iAcq = 1;
%%
for iAcq = 1:3
patient = num2str(T.patient(iAcq));
class = T.class(iAcq);
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

if T.xSup(iAcq) == 0
    % Manual cropping
    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 3.0])
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

    x_inf = rect(1); x_sup = rect(1)+rect(3);
    z_inf = rect(2); z_sup = rect(2)+rect(4);

else
    x_inf = T.xInf(iAcq); x_sup = T.xSup(iAcq);
    z_inf = T.zInf(iAcq); z_sup = T.zSup(iAcq);

end
% x_sup = 2.5; z_sup = 2.3; 
%%

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
% figure,
% plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
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
[Bn,~] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);

%% SWIFT
muBswift = 10^3.5; muCswift = 10^1;
% First iteration
[Bn,Cn,ite] = optimSwift(A1,A2,b(:),muBswift,muCswift,m,n,tol,ones(m,n));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
% w = (1-reject)*(1 ./ (1 + (bscMap./ratioCutOff).^10))+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
[Bn,Cn,ite] = optimSwift(A1w,A2w,bw,muBswift,muCswift,m,n,tol,w);
BSWIFT = reshape(Bn*NptodB,m,n);

% Weight map
figure('Units','centimeters', 'Position',[5 5 18 4]);
tl = tiledlayout(1,3, 'TileSpacing','tight', 'Padding','compact');

t2 = nexttile; 
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
colormap(t2,parula)
axis equal tight
title('BSC map')
c = colorbar;
c.Label.String = '\Delta BSC [db/cm]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t2,parula)
axis equal tight
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';


t2 = nexttile; 
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(t2,turbo)
axis equal tight
title('SWIFT')
% subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';



%% Mascaras
load(fullfile('newMasks',[patient,'.mat']));

[X,Z] = meshgrid(xFull,zFull);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskNoduleACS = interp2(X,Z,maskNodule,Xq,Zq, 'nearest');
maskThyroidACS = interp2(X,Z,maskThyroid,Xq,Zq, 'nearest');

se = strel('diamond',2);
maskThyroidACS = imerode(maskThyroidACS,se);
maskNoduleACS = imerode(maskNoduleACS,se);
% figure, imagesc(x_ACS,z_ACS,maskThyroidACS|maskNoduleACS)

patCol(iAcq) = {patient}; 
classCol(iAcq) = {class};
dataCols(iAcq,:) = [mean(BRTV(maskNoduleACS)), std(BRTV(maskNoduleACS)),...
    mean(BRSWTV(maskNoduleACS)), std(BRSWTV(maskNoduleACS)),...
    mean(BSWIFT(maskNoduleACS)), std(BSWIFT(maskNoduleACS)), ...
    mean(BRTV(maskThyroidACS)), std(BRTV(maskThyroidACS)),...
    mean(BRSWTV(maskThyroidACS)), std(BRSWTV(maskThyroidACS)),...
    mean(BSWIFT(maskThyroidACS)), std(BSWIFT(maskThyroidACS))];


fprintf("D ACS, TV: %.2f\n",...
    mean(BRTV(maskThyroidACS)) - mean(BRTV(maskNoduleACS)));
fprintf("D ACS, SWIFT: %.2f\n",...
    mean(BSWIFT(maskThyroidACS)) - mean(BSWIFT(maskNoduleACS)));
%% Overlay
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

[X,Z] = meshgrid(x,z);
roiSub = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units','centimeters', 'Position',[5 5 24 4.6])
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile();
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
title('B-mode')
ylim([0.1, 3])
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
hBm = colorbar;
hBm.Label.String = 'dB';
hBm.Location = 'westoutside';

nexttile,
[~,~,hColor] = imOverlayInterp(BmodeFull,BRTV,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('RSLD')
colorbar off
ylim([0.1, 3])
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskThyroid & roi,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [cm]')

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRSWTV,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWTV-ACE')
colorbar off
ylim([0.1, 3])
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskThyroid & roi,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [cm]')


nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BSWIFT,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWIFT')
% colorbar off
ylim([0.1, 3])
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskThyroid & roi,1,'w--')
hold off
xlabel('Lateral [cm]')
% hColor.Location = 'northoutside';
% hColor.Layout.Tile = 'northoutside';
hColor.Label.String = 'ACS [dB/cm/MHz]';
colormap(t1,'gray')
fontsize(gcf,9,'points')

%%
save_all_figures_to_directory(resultsDir,['pat',patient,'fig'],'svg');
close all

end

%%
infoTable = table(patCol',classCol',...
          'VariableNames',{'patient','type'});
dataTable = array2table(dataCols,...
    'VariableNames',{ ...
    'nod-TV-mean','nod-TV-std', ...
    'nod-SWTV-mean','nod-SWTV-std', ...
    'nod-SWIFT-mean','nod-SWIFT-std', ...
    'thy-TV-mean','thy-TV-std', ...
    'thy-SWTV-mean','thy-SWTV-std', ...
    'thy-SWIFT-mean','thy-SWIFT-std', ...
    });

writetable([infoTable,dataTable],fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
