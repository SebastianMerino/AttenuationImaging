% ======================================================================
% Script to explore weights in frequency
% Jan 12, 2024
% ======================================================================
%% PHANTOMSSS
clear,
% clc
close all

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];

targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\phantoms\ID316V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\phantoms\ID544V2\06-08-2023-Generic'];

rawFiles = dir([targetDir,'\*.rf']);
targetFiles = dir([targetDir,'\*.mat']);

resultsDir = fullfile(targetDir,'results','24-02-15','wFreq');
if ~exist("resultsDir","dir"); mkdir(resultsDir); end
tableName = 'phantoms.xlsx';

%% Constants
blocksize = 8;     % Block size in wavelengths
ratio_zx        = 12/8;

% freq_L = 2.5e6; freq_H = 7.5e6;
% freq_L = 1e6; freq_H = 10e6;
freq_L = 2e6; freq_H = 9e6;
fixedBW = true;
freq_C = 4.5e6;

overlap_pc      = 0.8;
x_inf = 0.1; x_sup = 3.8;
z_inf = 0.2; z_sup = 3.5;
NptodB = log10(exp(1))*20;
freqCutOff = db2mag(-15); % 15 is VERY GOOD

% Weights SWTV
aSNR = 1; bSNR = 0.1;
desvMin = 15;

% Weight map
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

groundTruthTargets = [0.52,0.55,0.74,0.61,0.75,0.97,0.95,0.95,0.55];

% Plotting constants
dynRange = [-50,0];
attRange = [0.4,1.1];
bsRange = [-20,20];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%% For looping each phantom
% iAcq = 8;
for iAcq = 1:8
switch mod(iAcq,3)
    % Optimal reg for BS 8x12
    case 0
        muBwfr = 10^2.5; muCwfr = 10^2; % Frequency wieght v1
        muBwfr2 = 10^3; muCwfr2 = 10^2.5; % Frequency wieght v2
    case 1
        muBwfr = 10^2.5; muCwfr = 10^0; % Frequency wieght v1
        muBwfr2 = 10^3; muCwfr2 = 10^0.5; % Frequency wieght v2

    case 2
        muBwfr = 10^2.5; muCwfr = 10^0; % Frequency wieght v1
        muBwfr2 = 10^3; muCwfr2 = 10^0.5; % Frequency wieght v2
end

%% Setting up
% Loading data
fprintf("Phantom no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(targetDir,targetFiles(iAcq).name));
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
sam1 = RF(:,:,1);

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z';
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

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

% Selecting BW 
windowing = tukeywin(nz/2,0.25);
NFFT = 2^(nextpow2(nz/2)+2);
[pxx,fpxx] = pwelch(sam1-mean(sam1),windowing,round(nz/4),NFFT,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum = meanSpectrum./max(meanSpectrum);
if ~fixedBW
    [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, freqCutOff);
end
figure,plot(fpxx/1e6,db(meanSpectrum))
xline([freq_L,freq_H]/1e6)
xline(freq_C/1e6, 'k--')
xlim([0 20])
ylim([-90 0])
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
grid on

% Frequency samples
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);


%% Generating Diffraction compensation
if true
    % Generating references
    att_ref = 0.53*f/8.686; % From phantom especifications
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
end

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

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
% rI = 0.6; rB = 1.2; % Both
% inc = (Xq - c1x).^2 + (Zq - c1z).^2 < rI^2;
% back = (Xq - c1x).^2 + (Zq - c1z).^2 > rB^2;
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

%% Checking spectrum
rInc = .9;
c1x = 1.93;
c1z = 1.97;
% incBmode = (Xq - c1x).^2 + (Zq - c1z).^2 < rInc^2;
% incAcs = (X - c1x).^2 + (Z - c1z).^2 < (rInc-0.1)^2;
% backAcs = (X - c1x).^2 + (Z - c1z).^2 >= (rInc+0.1)^2;
[backAcs, incAcs] = getRegionMasks(x_ACS,z_ACS,c1x,c1z,roiL,roiD,roiLz);
Spd = db(( abs(Sp) + abs(Sd) )/ 2);
Spd = Spd - max(Spd(:));

% %% Spectrum by depth
% spectUp = squeeze(mean(Spd(10,:,:),2));
% spectDown = squeeze(mean(Spd(35,:,:),2));
% 
% figure('Units','centimeters', 'Position',[5 5 10 10])
% plot(f,spectUp)
% hold on 
% plot(f,spectDown)
% hold off
% grid on
% legend('up','down')

%% SLD
test = reshape(b(repmat(incAcs,1,1,p)),sum(incAcs(:)),p);
sld1 = mean(test)' * NptodB /4/L;

test = reshape(b(repmat(backAcs,1,1,p)),sum(backAcs(:)),p);
sld2 = mean(test)' * NptodB /4/L;

acs1 = [f,ones(size(f))]\sld1;
fprintf('Attenuation is %.2f\n',acs1(1))
acs2 = [f,ones(size(f))]\sld2;
fprintf('Attenuation is %.2f\n',acs2(1))
figure('Units','centimeters', 'Position',[5 5 10 10])
plot(f,sld1)
hold on
plot(f,sld2)
plot(f,acs1(1)*f + acs1(2), 'k--')
plot(f,acs2(1)*f + acs2(2), 'k--')
hold off
grid on,
xlim([0,max(f)]), ylim([-3 10]),
xlabel('Frequency [MHz]')
ylabel('Att. [dB/cm]')
title('Mean SLD')
legend('Inc','Back')


 
%% NEW WEIGHTS
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

meanSpect = squeeze(mean(mean(Spd,2),1));
meanSpect = meanSpect - max(meanSpect);
wFreq = db2mag(meanSpect);

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2);

t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Spatial weights')
c = colorbar;
c.Label.String = '[a.u.]';
xlabel('x [cm]'), ylabel('z [cm]')

nexttile
plot(f,wFreq)
xlabel('Frequency [MHz]')
ylim([0 1])
xlim([f(1) f(end)])
grid on
title('Frequency weights')

W = repmat(w,[1 1 p]);
W = W.*reshape(wFreq,1,1,p);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BRWFR = (reshape(Bn*NptodB,m,n));
CRWFR = (reshape(Cn*NptodB,m,n));

AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'WFR1';
MetricsWFR(iAcq) = r;

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis image
title(['WFR, \mu=',num2str(muBwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRWFR, bsRange)
colormap(t3,parula)
axis image
title(['WFR, \mu=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

%% Frequency dependant weights
spectZF = squeeze(mean(Spd,2));
spectZF = spectZF - max(spectZF(:));
% figure,
% imagesc(f,z_ACS,spectZF)
% colorbar
% xlabel('f [MHz]')
% ylabel('z [mm]')
% 
spectT1 = -20;
spectT2 = -25;
wFreq = (spectZF>spectT1) + (spectZF<spectT1 & spectZF>spectT2).*...
    (spectZF - spectT2)./(spectT1 - spectT2);

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2);

t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Spatial weights')
c = colorbar;
c.Label.String = '[a.u.]';
xlabel('x [cm]'), ylabel('z [cm]')

nexttile,
imagesc(f,z_ACS,wFreq)
colorbar
xlabel('f [MHz]')
ylabel('z [mm]')
title('Frequency weights')
%%
W = repmat(w,[1 1 p]);
W = W.*reshape(wFreq,m,1,p);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr2,muCwfr2,m,n,tol,mask(:),w);
toc
BRWFR2 = (reshape(Bn*NptodB,m,n));
CRWFR2 = (reshape(Cn*NptodB,m,n));

AttInterp = interp2(X,Z,BRWFR2,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'WFR2';
MetricsWFR2(iAcq) = r;

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR2, attRange)
colormap(t2,turbo)
axis image
title(['WFR, \mu=',num2str(muBwfr2,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRWFR2, bsRange)
colormap(t3,parula)
axis image
title(['WFR, \mu=',num2str(muCwfr2,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


%%
figure('Units','centimeters', 'Position',[5 5 15 4]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';
fontsize(gcf,8,'points')

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWFR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('WFR')
% c = colorbar;
% c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWFR2, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off
%%
save_all_figures_to_directory(resultsDir,['phantom',num2str(iAcq),'fig']);
close all
end

results1 = struct2table(MetricsWFR);
results2 = struct2table(MetricsWFR2);

T = [results1;results2];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);

%%
