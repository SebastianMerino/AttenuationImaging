% ====================================================================== %
% Script for clinical data. 
% Created on June 27th, 2024
% ====================================================================== %
clear,clc
close all

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Liver_24_06_28\set1'];
targetDir = fullfile(baseDir,'sample');
refsDir = fullfile(baseDir,'ref');
resultsDir = fullfile(baseDir,'results','24-07-02-intento3');

% convertRfFilesToMat(targetDir)
% convertRfFilesToMat(refDir)

if (~exist(resultsDir,"dir")), mkdir(resultsDir); end
targetFiles = dir([targetDir,'\*.mat']);

% CHANGE ACCORDINGLY
maskReady = false;

%%
blocksize = 12;   % Axial block size in wavelengths
blocklines = 8;   % Num of lines, lateral block size
overlap_pc      = 0.8;

rect = [];

% Bandwidth
fixedBW = true;
ratio = db2mag(-30);
freq_L = 1.5e6; freq_H = 4.5e6;
% freq_L = 1.8e6; freq_H = 3.5e6;

% Weight parameters
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Optimal hetero phantom
% muBtv = 10^3; muCtv = 10^3;
% muBswtv = 10^2.5; muCswtv = 10^-0.5;
% muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
% muBwfr = 10^3; muCwfr = 10^0.5;

% Optimal thyroid
% muBtv = 10^2.5; muCtv = 10^2.5;
% muBswtv = 10^2.5; muCswtv = 10^-0.5;
% muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
% muBwfr = 10^3; muCwfr = 10^0.5;

muBtv = 10^3.5; muCtv = 10^3.5;
muBswtv = 10^3; muCswtv = 10^2.5;
muBtvl1 = 10^3.5; muCtvl1 = 10^2;
muBwfr = 10^3.5; muCwfr = 10^2;

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.5];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%%
dataCols = zeros(length(targetFiles),16);
for iAcq = 4:4 %1:length(targetFiles)

fprintf("Patient no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(targetDir,targetFiles(iAcq).name));

BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));

figure();
pcolor(xPolar*1e2,zPolar*1e2,BmodeFull)
colorbar;
clim(dynRange);
colormap gray
title('Bmode image')
ylabel('[cm]')
shading interp
axis equal ij tight



if maskReady
    load(fullfile(baseDir,'masks',targetFiles(iAcq).name)),
else
    maskLiver = true(size(rf));
end
xFull = th; % [deg]
zFull = (r-r(1))*1e2; % [cm]

rect = [];
% switch iAcq
%     case 1
%         rect = [-15.9636    4.0296   50.1058    3.9175];
%         % rect = [-12    4.0296   44    3.9175];
%     case 2
%         % rect = [-26.0459    4.0504   55.6052    3.8]; % cajita 4
%         % rect = [-26.0459    4.0504   55.6052    4.2]; % cajita 4
%         rect = [ -7.7145    3.9046   39.2597    4.6260];
%         % rect = [];
%     otherwise
%         rect = [];
% end

if isempty(rect)
    % Manual cropping
    figure('Units','centimeters', 'Position',[5 5 15 15]),
    imagesc(xFull,zFull,BmodeFull);
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('Lateral [deg]'); ylabel('Axial [cm]');
    % ylim([0.1 8])
    title('Liver')
    
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

disp(rect)
%%
x_inf = rect(1); x_sup = rect(1)+rect(3); 
z_inf = rect(2); z_sup = rect(2)+rect(4);
dz = (zFull(2) - zFull(1))/100;

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = rf(ind_z,ind_x);
Bmode = BmodeFull(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% Lateral samples
wx = round(blocklines*(1-overlap_pc));  % Between windows
nx = blocklines;                 % Window size
x0 = 1:wx:length(x)-nx;
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blocksize*wl*(1-overlap_pc)/dz ); % Between windows
nz = 2*round(blocksize*wl/dz /2); % Window size
% nz = 2*round(blocksize*wl/dz /2); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% BW from spectrogram
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = Bmode - max(Bmode(:));

fprintf('Frequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f lines, z: %.2f mm\n',blocklines,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = 0.54*f/NptodB; % From 20960001 _ID203458544
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
refFiles = dir([refsDir,'\*.mat']);
Nref = length(refFiles);
swrap = 0; % Correction factor for phantom data

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load([refsDir,'\',refFiles(iRef).name]);
    samRef = out.rf;
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

% compensation = zeros(m,n,p);

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref
% %% 
% figure,
% tiledlayout(1,2)
% t1 = nexttile; t2 = nexttile;
% for kk = 1:p
%     imagesc(t1,x_ACS,z_ACS, db(Sp(:,:,kk)), [100 200])
%     colorbar(t1)
% 
%     imagesc(t2,x_ACS,z_ACS, compensation(:,:,kk), [-2 2])
%     colorbar(t2)
%     pause(0.1)
% end
% figure,
% plot(f,squeeze(mean(mean(compensation,1),2))*NptodB)

%% Setting up

spectrumEnd = zeros(NFFT,1);
spectrumMid = zeros(NFFT,1);

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

        if ii == m
            spectrumEnd = spectrumEnd + tempSd/m;
        end
        if ii == round(m/2)
            spectrumMid = spectrumMid + tempSd/m;
        end

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


%% Power spectrum
normFactor = max(spectrumMid);
spectrumEnd = spectrumEnd/normFactor;
spectrumMid = spectrumMid/normFactor;

figure, hold on
plot(band/1e6,db(spectrumMid))
plot(band/1e6,db(spectrumEnd))
xline(freq_L/1e6, 'k--')
xline(freq_H/1e6, 'k--')
xlim([0, fs/2e6])
hold off
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
legend('Mid','End')
grid on

%% Mascaras
[X,Z] = meshgrid(xFull,zFull);
[Xq,Zq] = meshgrid(x_ACS,z_ACS);
maskLiverACS = interp2(X,Z,maskLiver,Xq,Zq, 'nearest');
maskMuscleACS = interp2(X,Z,~maskLiver,Xq,Zq, 'nearest');


se = strel('diamond',5); % 1 initially
maskLiverACS = imerode(maskLiverACS,se);
maskMuscleACS = imerode(maskMuscleACS,se);
%% Spectrum fit
% Heterogeneous
region1 = maskLiverACS;
sld1 = squeeze(sum(sum(b.*region1,1),2))/sum(region1(:)) * NptodB /4/L;
acs1 = f\sld1;

acs2 = [f ones(size(f))]\sld1;
figure, plot(f,sld1)
hold on
% plot(f,sld2)
plot(f,acs1*f, 'k--')
plot(f,[f ones(size(f))]*acs2, 'r--')
hold off
grid on,
xlim([0,max(f)]), ylim([0 5]),
xlabel('Frequency [MHz]')
ylabel('Att. [dB/cm]')
title('Mean SLD')
legend('Liver',sprintf('ACS ZC = %.2f\n',acs1), ...
    sprintf('ACS NZC = %.2f\n',acs2(1)), 'Location','northwest')


% region1 = maskMuscleACS;
% sld1 = squeeze(sum(sum(b.*region1,1),2))/sum(region1(:)) * NptodB /4/L;
% acs1 = f\sld1;
% acs2 = [f ones(size(f))]\sld1;
% figure, plot(f,sld1)
% hold on
% % plot(f,sld2)
% % plot(f,acs1*f, 'k--')
% % plot(f,[f ones(size(f))]*acs2, 'r--')
% hold off
% grid on,
% xlim([0,max(f)]), ylim([0 15]),
% xlabel('Frequency [MHz]')
% ylabel('Att. [dB/cm]')
% title('Mean SLD')
% legend('Muscle',sprintf('ACS ZC = %.2f\n',acs1), ...
%     sprintf('ACS NZC = %.2f\n',acs2(1)), 'Location','northwest')

%% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);
disp(mean(BRTV(:)))
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

% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBwfr,muCwfr,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);

% Weight map
figure('Units','centimeters', 'Position',[5 5 18 4]);
tl = tiledlayout(1,3, 'TileSpacing','tight', 'Padding','compact');

t2 = nexttile; 
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
colormap(t2,parula)
% axis equal tight
title('BSC map')
c = colorbar;
c.Label.String = '\Delta BSC [db/cm]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t2,parula)
% axis equal tight
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';


t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
% axis equal tight
title('SWIFT')
% subtitle(['\mu_b=',num2str(muBtvl1,2),', \mu_c=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';



%% Mascaras
[~,Z] = meshgrid(x_ACS,z_ACS);
[~,Zbm] = meshgrid(xFull,zFull);
maskLiver = Zbm>4.7; 
maskLiverACS = Z>4.7;
dataCols(iAcq,:) = ...
    [mean(BRTV(maskLiverACS),'all'), std(BRTV(maskLiverACS),[],'all'),...
    mean(BRSWTV(maskLiverACS),'all'), std(BRSWTV(maskLiverACS),[],'all'),...
    mean(BRTVL1(maskLiverACS),'all'), std(BRTVL1(maskLiverACS),[],'all'),...
    mean(BRWFR(maskLiverACS),'all'), std(BRWFR(maskLiverACS),[],'all'),...
    mean(BRTV(maskMuscleACS),'all'), std(BRTV(maskMuscleACS),[],'all'),...
    mean(BRSWTV(maskMuscleACS),'all'), std(BRSWTV(maskMuscleACS),[],'all'),...
    mean(BRTVL1(maskMuscleACS),'all'), std(BRTVL1(maskMuscleACS),[],'all'),...
    mean(BRWFR(maskMuscleACS),'all'), std(BRWFR(maskMuscleACS),[],'all')];


%% Overlay
yLimits = [0,10];

[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units','centimeters', 'Position',[5 5 24 8])
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile();
imagesc(xFull,zFull,BmodeFull,dynRange); % axis image; 
title('B-mode')
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
% contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
xlabel('Lateral [deg]')
ylabel('Axial [cm]')
hBm = colorbar;
hBm.Label.String = 'dB';
hBm.Location = 'westoutside';

nexttile,
[~,~,hColor] = imOverlayInterp(BmodeFull,BRTV,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('RSLD')
axis normal
colorbar off
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [deg]')

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRSWTV,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWTV-ACE')
colorbar off
axis normal
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [deg]')


nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWIFT')
axis normal
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
xlabel('Lateral [deg]')
% hColor.Location = 'northoutside';
% hColor.Layout.Tile = 'northoutside';
hColor.Label.String = 'ACS [dB/cm/MHz]';
colormap(t1,'gray')
fontsize(gcf,9,'points')


%%
figure('Units','centimeters', 'Position',[5 5 14 8])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','tight')

t1 = nexttile;
imagesc(xFull,zFull,BmodeFull,dynRange); % axis image; 
title('B-mode')
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
% contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
xlabel('Lateral [deg]')
ylabel('Axial [cm]')
hBm = colorbar;
hBm.Label.String = 'dB';
hBm.Location = 'westoutside';

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRTVL1,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('TVL1')
axis normal
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
contour(xFull,zFull,maskLiver,1,'w--')
% contour(x_ACS,z_ACS,maskMuscleACS,1,'w--')
% contour(x_ACS,z_ACS,maskLiverACS,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [deg]')


% hColor.Layout.Tile = 'east';
% hColor.Label.String = 'ACS [dB/cm/MHz]';
colormap(t1,'gray')
fontsize(gcf,9,'points')


% ylabel('z [cm]')



%%
save_all_figures_to_directory(resultsDir,['liverAcq',num2str(iAcq),'fig'],'fig');
close all
save(fullfile(resultsDir,['rectAcq',num2str(iAcq),'.mat']),'rect')

end
%%
% infoTable = table(patCol',classCol',...
%           'VariableNames',{'patient','type'});
dataTable = array2table(dataCols,...
    'VariableNames',{ ...
    'liver-TV-mean','liver-TV-std', ...
    'liver-SWTV-mean','liver-SWTV-std', ...
    'liver-TVL1-mean','liver-TVL1-std', ...
    'liver-WFR-mean','liver-WFR-std', ...
    'muscle-TV-mean','muscle-TV-std', ...
    'muscle-SWTV-mean','muscle-SWTV-std', ...
    'muscle-TVL1-mean','muscle-TVL1-std', ...
    'muscle-WFR-mean','muscle-WFR-std', ...
    });
disp(dataTable)
tableName = 'clinicalLiver.xlsx';
writetable(dataTable,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
