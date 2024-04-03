% ====================================================================== %
% Script for plotting figures and results for simulation on circular
% inclusions, one resembling thyroid tissue
% Created on Jan 31, 2024
% ====================================================================== %

clear,clc

% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\simulations_processed\inc_journal'];
% resultsDir = 'C:\Users\smerino.C084288\Pictures\JOURNAL\24-02-20\BS_8_12';

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\inc_journal'];
resultsDir = 'C:\Users\sebas\Pictures\Journal2024\24-04-02';

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];


tableName = 'simuInc.xlsx';
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

targetFiles = dir([targetDir,'\rf*.mat']);
% targetFiles = targetFiles(2:3);
refFiles = dir([refDir,'\rf*.mat']);

if (~exist(resultsDir,"dir")), mkdir(resultsDir); end

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 8;     % Block size in wavelengths
ratio_zx        = 12/8;

freq_L = 3e6; freq_H = 8e6; % GOOD
% freq_L = 3e6; freq_H = 9e6; % Also GOOD

overlap_pc      = 0.8;
referenceAtt    = 0.6;

% Weights SWTV
aSNR = 1; bSNR = 0.1;
desvMin = 15;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% groundTruthThyroid = [0.8,1.5];
% groundTruthNodule = [1.5,0.8];
groundTruthThyroid = [0.8,1.5,0.8,1.5];
groundTruthNodule = [1.5,0.8,1.5,0.8];


%attRange = [0.6 1.7];
attRange = [0.4,1.7];

% Plotting
dynRange = [-60,0];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% For looping
% figure('Units','centimeters', 'Position',[5 5 25 8]);
% tl = tiledlayout(2,5, "Padding","tight");

for iAcq = 1:4

load(fullfile(targetDir,targetFiles(iAcq).name));

switch iAcq
    % OPTIMAL WEIGHTS FOR BS 8x12
    case 1
        muBtv = 10^3; muCtv = 10^1;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3; muCtvl1 = 10^1;
        muBwfr = 10^4; muCwfr = 10^1.5;
        % muBwfr = 10^3.5; muCwfr = 10^1.5; % WRONG WINDOW
    case 2
        muBtv = 10^3; muCtv = 10^1;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3; muCtvl1 = 10^0;
        muBwfr = 10^3.5; muCwfr = 10^-0.5;
        % muBwfr = 10^3; muCwfr = 10^-0.5; % WRONG WINDOW
    case 3
        muBtv = 10^3.5; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBwfr = 10^3.5; muCwfr = 10^2;
    case 4
        muBtv = 10^3.5; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBwfr = 10^3.5; muCwfr = 10^2;
end

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);

%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.3; z_sup = 3.9;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
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
% nz = 2*round(blocksize*wl/dz /2); % Window size
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% BW and spectrogram
% ratio = db2mag(-50);
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% figure,plot(fpxx/1e6,meanSpectrum)
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
% xline([freq_L,freq_H]/1e6)
% xlim([0 15])
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% grid on

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
% dynRange = [-40 -10];
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));
% figure, imagesc(x,z,Bmode);
% axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest rows: %i, col: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = referenceAtt*(f.^1.05)/(20*log10(exp(1)));
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
% windowing = hamming(nz/2);
windowing = windowing*ones(1,nx);

% For looping
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    load(fullfile(refDir,refFiles(iRef).name),"rf","medium");
    % disp(mean(medium.alpha_coeff(:)))
    samRef = rf;
    samRef = samRef(ind_z,ind_x); % Cropping
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
%cuec = zeros([NFFT,1]);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = sam1(zd:zd+nz/2-1,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);
%         if (ii==floor(m/2))
%             cuec = cuec + tempSp;
%         end
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% Setting Up

% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

%% TV
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
cz = 2; cx = 0;
rInc = 0.7;
noduleMask = (Zq-cz).^2 + (Xq-cx).^2 <= rInc^2;
thyroidMask = ~noduleMask;

c1x = 0; c1z = 2;
roiL = 0.7; roiD = 0.5;
roiLz = 1.2;
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
% [thyroidMask,noduleMask] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanTop = mean(AttInterp(thyroidMask),"omitnan");
r.stdTop = std(AttInterp(thyroidMask),"omitnan");
r.meanBottom = mean(AttInterp(noduleMask),"omitnan");
r.stdBottom = std(AttInterp(noduleMask),"omitnan");
r.biasTop = mean( AttInterp(thyroidMask) - groundTruthThyroid(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(noduleMask) - groundTruthNodule(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(thyroidMask) - groundTruthThyroid(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(noduleMask) - groundTruthNodule(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;


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

AttInterp = interp2(X,Z,BRSWTV,Xq,Zq);
r.meanTop = mean(AttInterp(thyroidMask),"omitnan");
r.stdTop = std(AttInterp(thyroidMask),"omitnan");
r.meanBottom = mean(AttInterp(noduleMask),"omitnan");
r.stdBottom = std(AttInterp(noduleMask),"omitnan");
r.biasTop = mean( AttInterp(thyroidMask) - groundTruthThyroid(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(noduleMask) - groundTruthNodule(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(thyroidMask) - groundTruthThyroid(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(noduleMask) - groundTruthNodule(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% TVL1
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanTop = mean(AttInterp(thyroidMask),"omitnan");
r.stdTop = std(AttInterp(thyroidMask),"omitnan");
r.meanBottom = mean(AttInterp(noduleMask),"omitnan");
r.stdBottom = std(AttInterp(noduleMask),"omitnan");
r.biasTop = mean( AttInterp(thyroidMask) - groundTruthThyroid(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(noduleMask) - groundTruthNodule(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(thyroidMask) - groundTruthThyroid(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(noduleMask) - groundTruthNodule(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% WFR
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Computing weights
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% %%
% w = ones(size(BRTVL1));
% borderMask = (Z-cz).^2 + (X-cx).^2 <= (rInc-0.15)^2 | ...
%     (Z-cz).^2 + (X-cx).^2 >= (rInc+0.15)^2 ;
% w(borderMask) = 1;
% w(~borderMask) = 0.1;
% w(18:24,:) = 1;
% 
% figure, imagesc(x_ACS,z_ACS,w)
%% -------
% Setting up new system
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Method
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);
CRWFR = reshape(Cn*NptodB,m,n);


AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanTop = mean(AttInterp(thyroidMask),"omitnan");
r.stdTop = std(AttInterp(thyroidMask),"omitnan");
r.meanBottom = mean(AttInterp(noduleMask),"omitnan");
r.stdBottom = std(AttInterp(noduleMask),"omitnan");
r.biasTop = mean( AttInterp(thyroidMask) - groundTruthThyroid(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(noduleMask) - groundTruthNodule(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(thyroidMask) - groundTruthThyroid(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(noduleMask) - groundTruthNodule(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;

% TESTING
%%
attIdeal{iAcq} = ones(m,n)*groundTruthThyroid(iAcq);
incACS = (Z-cz).^2 + (X-cx).^2 <= (rInc)^2;
attIdeal{iAcq}(incACS) = groundTruthNodule(iAcq);

%% Plotting
figure('Units','centimeters', 'Position',[5 5 20 4]);
tiledlayout(1,6, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
% c = colorbar(t1, 'westoutside');
% c.Label.String = '[db]';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(x_ACS,z_ACS,attIdeal{iAcq},attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('TV')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% hold off

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% hold off

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% hold off

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t4,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% hold off

fontsize(gcf,8,'points')

%%
figure('Units','centimeters', 'Position',[5 5 10 8])
tl = tiledlayout(2,2, "Padding","tight", 'TileSpacing','tight');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'eastoutside');
c.Label.String = '[dB]';
title('Bmode')
xlabel('Lateral [cm]')
ylabel('Axial [cm]')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('\DeltaBSC')
c = colorbar(t2, 'eastoutside');
c.Label.String = '[dB/cm]';
xlabel('Lateral [cm]')
%ylabel('Axial [cm]')

t3 = nexttile([1,2]); 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar(t3, 'eastoutside');
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
c.Label.String = '[a.u.]';

fontsize(gcf,8,'points')

% figure,
% [~,hB,hColor] = imOverlayInterp(Bmode,w,[-50 0],[0 1],0.2,...
%     x_ACS,z_ACS,ones(size(Bmode)),x,z);
% colormap(hot)
%%
axialTV{iAcq} = mean(BRTV(:,41:49),2);
axialSWTV{iAcq} = mean(BRSWTV(:,41:49),2);
axialTVL1{iAcq} = mean(BRTVL1(:,41:49),2);
axialWFR{iAcq} = mean(BRWFR(:,41:49),2);



%%
lateralTV{iAcq} = mean(BRTV(24:26,:),1);
lateralSWTV{iAcq} = mean(BRSWTV(24:26,:),1);
lateralTVL1{iAcq} = mean(BRTVL1(24:26,:),1);
lateralWFR{iAcq} = mean(BRWFR(24:26,:),1);



end

save_all_figures_to_directory(resultsDir,'thyroidSimFig');
close all
%%
layered = load('simuLayered.mat');
figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile([1,2]),
plot(layered.z_ACS, layered.axialTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(layered.z_ACS, layered.axialSWTV{1}, 'r', 'LineWidth',1),
plot(layered.z_ACS, layered.axialTVL1{1}, 'b:', 'LineWidth',1.5),
plot(layered.z_ACS, layered.axialWFR{1}, 'b', 'LineWidth',1),
plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([layered.z_ACS(1) layered.z_ACS(end)])
%title('Axial profile')
xlabel('Axial [cm]')
ylabel('ACS [dB/cm/MHz]')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northwest') 

figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile([1,2]),
plot(layered.z_ACS, layered.axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(layered.z_ACS, layered.axialSWTV{2}, 'r', 'LineWidth',1),
plot(layered.z_ACS, layered.axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(layered.z_ACS, layered.axialWFR{2}, 'b', 'LineWidth',1),
plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([layered.z_ACS(1) layered.z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')
%legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
% figure('Units','centimeters', 'Position',[5 5 12 12])
nexttile,
plot(z_ACS, axialTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{1}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{1}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{1}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal{1}(:,41:49),2), 'k--')
hold off
grid on
ylim([0.6 1.7])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 


nexttile,
plot(x_ACS, lateralTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(x_ACS, lateralSWTV{1}, 'r', 'LineWidth',1),
plot(x_ACS, lateralTVL1{1}, 'b:', 'LineWidth',1.5),
plot(x_ACS, lateralWFR{1}, 'b', 'LineWidth',1),
plot(x_ACS,mean(attIdeal{1}(24:26,:),1), 'k--')
hold off
grid on
ylim([0.6 1.7])
xlim([x_ACS(1) x_ACS(end)])
%title('Lateral profile')
xlabel('Lateral [cm]')
%legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile,
plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal{2}(:,41:49),2), 'k--')
hold off
grid on
ylim([0.6 1.7])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 


nexttile,
plot(x_ACS, lateralTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(x_ACS, lateralSWTV{2}, 'r', 'LineWidth',1),
plot(x_ACS, lateralTVL1{2}, 'b:', 'LineWidth',1.5),
plot(x_ACS, lateralWFR{2}, 'b', 'LineWidth',1),
plot(x_ACS,mean(attIdeal{2}(24:26,:),1), 'k--')
hold off
grid on
ylim([0.6 1.7])
xlim([x_ACS(1) x_ACS(end)])
% title('Lateral profile')
xlabel('Lateral [cm]')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

fontsize(gcf,8,'points')

save_all_figures_to_directory(resultsDir,'SimAxialLat');

%% ALL IN ONE FIGURE
% layered = load('simuLayered.mat');
% figure('Units','centimeters', 'Position',[5 5 10 15])
% tiledlayout(4,2, 'TileSpacing','compact', 'Padding','compact')
% nexttile([1,2]),
% plot(layered.z_ACS, layered.axialTV{1}, 'r:', 'LineWidth',1.5),
% hold on
% plot(layered.z_ACS, layered.axialSWTV{1}, 'r', 'LineWidth',1),
% plot(layered.z_ACS, layered.axialTVL1{1}, 'b:', 'LineWidth',1.5),
% plot(layered.z_ACS, layered.axialWFR{1}, 'b', 'LineWidth',1),
% plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
% hold off
% grid on
% ylim([0.4 1.4])
% xlim([layered.z_ACS(1) layered.z_ACS(end)])
% %title('Axial profile')
% xlabel('Axial [cm]')
% ylabel('ACS [dB/cm/MHz]')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northwest') 
% 
% nexttile([1,2]),
% plot(layered.z_ACS, layered.axialTV{2}, 'r:', 'LineWidth',1.5),
% hold on
% plot(layered.z_ACS, layered.axialSWTV{2}, 'r', 'LineWidth',1),
% plot(layered.z_ACS, layered.axialTVL1{2}, 'b:', 'LineWidth',1.5),
% plot(layered.z_ACS, layered.axialWFR{2}, 'b', 'LineWidth',1),
% plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
% hold off
% grid on
% ylim([0.4 1.4])
% xlim([layered.z_ACS(1) layered.z_ACS(end)])
% %title('Axial profile')
% ylabel('ACS [dB/cm/MHz]')
% xlabel('Axial [cm]')
% %legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% 
% % figure('Units','centimeters', 'Position',[5 5 12 12])
% nexttile,
% plot(z_ACS, axialTV{1}, 'r:', 'LineWidth',1.5),
% hold on
% plot(z_ACS, axialSWTV{1}, 'r', 'LineWidth',1),
% plot(z_ACS, axialTVL1{1}, 'b:', 'LineWidth',1.5),
% plot(z_ACS, axialWFR{1}, 'b', 'LineWidth',1),
% plot(z_ACS,mean(attIdeal{1}(:,41:49),2), 'k--')
% hold off
% grid on
% ylim([0.6 1.7])
% xlim([z_ACS(1) z_ACS(end)])
% %title('Axial profile')
% ylabel('ACS [dB/cm/MHz]')
% xlabel('Axial [cm]')
% % legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% 
% % figure('Units','centimeters', 'Position',[5 5 12 12])
% % tiledlayout(2,1)
% nexttile,
% plot(x_ACS, lateralTV{1}, 'r:', 'LineWidth',1.5),
% hold on
% plot(x_ACS, lateralSWTV{1}, 'r', 'LineWidth',1),
% plot(x_ACS, lateralTVL1{1}, 'b:', 'LineWidth',1.5),
% plot(x_ACS, lateralWFR{1}, 'b', 'LineWidth',1),
% plot(x_ACS,mean(attIdeal{1}(24:26,:),1), 'k--')
% hold off
% grid on
% ylim([0.6 1.7])
% xlim([x_ACS(1) x_ACS(end)])
% %title('Lateral profile')
% xlabel('Lateral [cm]')
% %legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% 
% nexttile,
% plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
% hold on
% plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
% plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
% plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
% plot(z_ACS,mean(attIdeal{2}(:,41:49),2), 'k--')
% hold off
% grid on
% ylim([0.6 1.7])
% xlim([z_ACS(1) z_ACS(end)])
% %title('Axial profile')
% ylabel('ACS [dB/cm/MHz]')
% xlabel('Axial [cm]')
% % legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% 
% nexttile,
% plot(x_ACS, lateralTV{2}, 'r:', 'LineWidth',1.5),
% hold on
% plot(x_ACS, lateralSWTV{2}, 'r', 'LineWidth',1),
% plot(x_ACS, lateralTVL1{2}, 'b:', 'LineWidth',1.5),
% plot(x_ACS, lateralWFR{2}, 'b', 'LineWidth',1),
% plot(x_ACS,mean(attIdeal{2}(24:26,:),1), 'k--')
% hold off
% grid on
% ylim([0.6 1.7])
% xlim([x_ACS(1) x_ACS(end)])
% % title('Lateral profile')
% xlabel('Lateral [cm]')
% % legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% fontsize(gcf,8,'points')
% save_all_figures_to_directory(resultsDir,'SimAxialLat');


%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasTop)
disp(results2.biasTop)
disp(results3.biasTop)
disp(results4.biasTop)

disp('Bias Bottom')
disp(results1.biasBottom)
disp(results2.biasBottom)
disp(results3.biasBottom)
disp(results4.biasBottom)

disp('RMSE Top')
disp(results1.rmseTop)
disp(results2.rmseTop)
disp(results3.rmseTop)
disp(results4.rmseTop)

disp('RMSE Bottom')
disp(results1.rmseBottom)
disp(results2.rmseBottom)
disp(results3.rmseBottom)
disp(results4.rmseBottom)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)


T = [results1;results2;results3;results4];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);


%%
