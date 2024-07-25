% ====================================================================== %
% Script for clinical data.
% Created on June 19th, 2024
% ====================================================================== %

setup,

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'LiverAcquisition_24_09_07'];

sampleDir = fullfile(baseDir,'samples');
refsDir = fullfile(baseDir,'refs','joined');
resultsDir = 'C:\Users\sebas\Pictures\Journal2024\liver';
figsDir = 'C:\Users\sebas\Pictures\Journal2024\liver';

[~,~,~] = mkdir(resultsDir);
[~,~,~] = mkdir(figsDir);

sampleFiles = dir(fullfile(sampleDir,'015*.mat'));
sampleFiles = [sampleFiles; dir(fullfile(sampleDir,'016*.mat'))];

for iFile = 1:length(sampleFiles)

samName = sampleFiles(iFile).name(1:end-4);

%% ACS estimation hyperparameters
blocksize = 12;   % Axial block size in wavelengths
blocklines = 8;   % Num of lines, lateral block size
overlap_pc      = 0.8;

% Bandwidth
fixedBW = true;
ratio = db2mag(-30);
freq_L = 1.5e6; freq_H = 4e6;

% Weight parameters
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Reg weights (same as in thyroid cases)
muBtv = 10^3; muCtv = 10^3;
muBswtv = 10^2.5; muCswtv = 10^-0.5;
muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
muBwfr = 10^3; muCwfr = 10^0.5;

% Plotting constants
dynRange = [-70,0];
attRange = [0,1.5];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% Loading file and variables
load(fullfile(sampleDir,samName+".mat"));
fprintf("Loading sample %s \n", samName)
RcvData = cell2mat(RcvData);
n_frame = size(RcvData,3); % Frame selector
RcvData = RcvData(:, :, n_frame); % Select frame of RF Data

% Additional variables
central_freq = Receive(1).demodFrequency*1e6; % Central frequency of pulse
fs = Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
n_pulses = P.numRays; % number of pulses
n_elements = Trans.numelements; % number of transducer elements
num_samples = Receive(1).endSample - Receive(1).startSample +1; % samples per channel
sound_speed = Resource.Parameters.speedOfSound; % [m/s]
wvl = sound_speed/central_freq; % [m] wavelength in meters
scalemm2wvl = 1/wvl;

% Initialize variables
rf_channel = zeros(num_samples , n_elements, n_pulses);
rx_apods = zeros(1, n_elements, n_pulses);
rf = zeros(num_samples, n_pulses);

%% Organize data
for n = 1:n_pulses % Iterate through pulses
    rf_channel(:, :, n) = RcvData(Receive(n).startSample:Receive(n).endSample, :);
end


% Acommodate to time delays and rf signals
focus = 20/1000;
t = (0:(num_samples-1))/fs; % [sec.] time domain 0:T:(N_sample-1)*T
[rx_delays] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl);


% Dynamic Aperture
f_num = 3;
z = sound_speed*t/2;
elem_pitch = Trans.spacingMm*1e-3;
maxAprSz = 32;
dyn_aperture = zeros(length(z), n_elements, n_pulses);
for n = 1:n_pulses
    for z_i = 1:length(z)
        a = z(z_i)/(2*f_num);
        hlfAprSz = floor(a / elem_pitch);
        if (hlfAprSz > maxAprSz/2)
            hlfAprSz = floor(maxAprSz / 2);
        end
        a_i = -hlfAprSz: hlfAprSz;    % aperture indices
        fulAprSz = 2*hlfAprSz + 1;
        aper_center = n;
        aper = aper_center + a_i;
        aper = aper(aper>=1);
        aper = aper(aper<=128);
        dyn_aperture(z_i, aper, n) = 1;
    end
end


%% Beamforming

% Delay-and-sum
for n = 1:n_pulses
    % Delaying
    for e = 1:n_elements
        % rx_apods(e, n);
        rf_channel(:, e, n) = ...
            interp1(t, rf_channel(:, e, n), rx_delays(:,e, n), 'linear', 0);
    end
    rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
    % Summing
    rf(:, n) = sum(rf_channel(:, :, n), 2);
end

%% B-Mode and coordinates
b_mode = 20*log10(abs(hilbert(rf)));
b_mode = b_mode-max(b_mode(:));

param = getparam('C5-2v');
siz = size(rf);
z = sound_speed*t/2;
zmax = z(end);
R = param.radius;
p = param.pitch;
N = param.Nelements;
L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
d = sqrt(R^2-L^2/4); % apothem
z0 = -d;

th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
r = linspace(R+p,-z0+zmax,siz(1));

% To Polar Coordinates
[xPolar,zPolar, z0Polar] = impolgrid(size(b_mode), z(end),param);


%% Selecting ROI
BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));

xFull = th; % [deg]
r0 = r(1);
zFull = (r-r0)*1e2; % [cm]

% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(xFull,zFull,BmodeFull,dynRange);
% ylim([0 10])
% colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
% title('Liver')
% 
% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
%     if strcmp(confirmation,'Cancel')
%         disp(rect)
%         break
%     end
% end
% close,

load(fullfile(resultsDir,samName+".mat"),'rect');


%% Setting up
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);
dz = (zFull(2) - zFull(1))/100;

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = rf(ind_z,ind_x);
Bmode = BmodeFull(ind_z,ind_x);
Bmode = Bmode - max(Bmode(:));

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

% Plot
[pxx,fpxx] = pwelch(sam1-mean(sam1),nz,nz-wz,nz,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum(1) = 0;
figure,
plot(fpxx/1e6,db(meanSpectrum/max(meanSpectrum))),grid on
xlim([0, fs/2e6])
hold on
xline(freq_L/1e6, 'k--')
xline(freq_H/1e6, 'k--')
hold off
xlabel('Frequency [MHz]')
ylabel('Magnitude [dB]')
saveas(gcf,fullfile(figsDir,"sample"+samName+"_Spectrum.png"))
%close

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

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref

%% Spectral log difference

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

%% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BR = reshape(Bn*NptodB,m,n);

%% SWTV-ACE
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
BSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

%% SWIFT

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
BSWIFT = reshape(Bn*NptodB,m,n);

% %% Weight plot
% figure('Units','centimeters', 'Position',[5 5 18 4]);
% tl = tiledlayout(1,3, 'TileSpacing','tight', 'Padding','compact');
% t2 = nexttile;
% imagesc(x_ACS,z_ACS,bscMap, [-20 20])
% colormap(t2,parula)
% title('BSC map')
% c = colorbar;
% c.Label.String = '\Delta BSC [db/cm]';
% t2 = nexttile;
% imagesc(x_ACS,z_ACS,w, [0 1])
% colormap(t2,parula)
% title('Weights')
% c = colorbar;
% c.Label.String = '[a.u.]';
% t2 = nexttile;
% imagesc(x_ACS,z_ACS,BSWIFT, attRange)
% colormap(t2,turbo)
% title('SWIFT')
% c = colorbar;
% c.Label.String = 'ACS [db/cm/MHz]';

%% SLD and SWIFT fit
regionMaskAcs = ones(m,n);
sldLine = squeeze(sum(sum(b.*regionMaskAcs,2),1))/4/L*NptodB/sum(regionMaskAcs(:));
fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2),
nexttile
plot(f,sldLine),
hold on,
plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), '--')
hold off,
grid on,
xlim([0,freq_H*1.1]/1e6),
ylim([-1 5]),
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Mean SLD')
legend('Liver',sprintf('ACS ZC = %.2f\n',fit1), ...
    sprintf('ACS NZC = %.2f\n',fit2(1)), 'Location','northwest')

sldLine = squeeze(sum(sum(b.*wExt,2),1))/4/L*NptodB/sum(wExt(:));
fit1 = f\sldLine;
fit2 = [f ones(length(f),1)]\sldLine;
nexttile,
plot(f,sldLine),
hold on,
plot(f,fit1*f, '--')
plot(f,fit2(1)*f + fit2(2), '--')
hold off,
grid on,
xlim([0,freq_H*1.1]/1e6),
ylim([-1 5]),
xlabel('Frequency [MHz]')
ylabel('Attenuation [dB/cm]')
title('Weighted mean SLD')
legend('Liver',sprintf('ACS ZC = %.2f\n',fit1), ...
    sprintf('ACS NZC = %.2f\n',fit2(1)), 'Location','northwest')

saveas(gcf,fullfile(figsDir,"sample"+samName+"_sldFit.png"))
close
%% Overlay
yLimits = [0, 12];
meanRsld = mean(BR,'all');
meanSwtv = mean(BSWTV,'all');
meanSwift = mean(BSWIFT,'al');

[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);

figure('Units','centimeters', 'Position',[5 5 20 7])
tiledlayout(1,4, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile();
imagesc(xFull,zFull,BmodeFull,dynRange); % axis image;
title('B-mode')
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
xlabel('Lateral [deg]')
ylabel('Axial [cm]')
hBm = colorbar('Ticks',-60:20:0);
hBm.Label.String = 'dB';
hBm.Location = 'westoutside';

nexttile,
[~,~,hColor] = imOverlayInterp(BmodeFull,BR,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('RSLD')
subtitle(['ACS = ',num2str(meanRsld,2),' dB/cm/MHz'])
axis normal
colorbar off
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [deg]')

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BSWTV,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWTV-ACE')
subtitle(['ACS = ',num2str(meanSwtv,2),' dB/cm/MHz'])
colorbar off
axis normal
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
% axis off
%xlabel('x [cm]')
xlabel('Lateral [deg]')


nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BSWIFT,dynRange,attRange,0.7,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('SWIFT')
subtitle(['ACS = ',num2str(meanSwift,2),' dB/cm/MHz'])
axis normal
ylim(yLimits)
hold on
contour(xFull,zFull,roi,1,'w--')
hold off
xlabel('Lateral [deg]')
% hColor.Location = 'northoutside';
% hColor.Layout.Tile = 'northoutside';
hColor.Label.String = 'ACS [dB/cm/MHz]';
colormap(t1,'gray')
% fontsize(gcf,9,'points')

exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_Polar.png"), ...
    'Resolution','300')

%% Plot in cartesian cords
[TH_acs,R_acs] = meshgrid(-x_ACS*pi/180 + pi/2,z_ACS/100 + r0);
[xPolarACS,zPolarACS] = pol2cart(TH_acs,R_acs);
zPolarACS = zPolarACS + z0Polar;


figure('Units','centimeters', 'Position',[5 5 6 6]);
[ax1,~] = imOverlayPolar(BmodeFull,ones(m,n),dynRange,attRange,0, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'B-mode')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-9 9]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_Bm.png"), ...
    'Resolution','300')

figure('Units','centimeters', 'Position',[5 5 6 6]);
[ax1,~] = imOverlayPolar(BmodeFull,BR,dynRange,attRange,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'RSLD')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-9 9]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_rsld.png"), ...
    'Resolution','300')


figure('Units','centimeters', 'Position',[5 5 6 6]);
[ax1,~] = imOverlayPolar(BmodeFull,BSWTV,dynRange,attRange,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'SWTV-ACE')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-9 9]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_swtv.png"), ...
    'Resolution','300')

figure('Units','centimeters', 'Position',[5 5 6 6]);
[ax1,~] = imOverlayPolar(BmodeFull,BSWIFT,dynRange,attRange,0.7, ...
    xPolar,zPolar,xPolarACS,zPolarACS);
yticks(ax1,[4 8 12 16])
title(ax1,'SWIFT')
xlabel(ax1,'Lateral [cm]'), ylabel(ax1,'Axial [cm]')
xlim([-9 9]), ylim([0 18])
hold on
contour(xPolar*1e2, zPolar*1e2, roi,1,'w--')
hold off
exportgraphics(gcf,fullfile(figsDir,"sample"+samName+"_swift.png"), ...
    'Resolution','300')
pause(0.5)
close all


%% Saving ACS maps and ROI
save(fullfile(resultsDir,samName+".mat"), ...
    'BR','BSWTV','BSWIFT','rect')

end

%%
resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles = length(resultFiles);


varNames = {'acs','method','patient','sample'};
varTypes = {'double','string','int16','int16'};
T = table('Size',[nFiles*3,4], 'VariableTypes',varTypes,'VariableNames',varNames);


resultFiles = dir(fullfile(resultsDir,'*.mat'));
nFiles = length(resultFiles);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample};
    T(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample};
    T(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample};
end
writetable(T,fullfile(resultsDir,'meanACS.xlsx')) 

%% Auxiliary functions

% Get delays
function [t_delay] = getRXDelays(Trans, t, n_elements, n_pulses, sound_speed, wvl)

t_delay = zeros(length(t), n_elements, n_pulses);
% (x, z) [m] Obtain positions of center of every element
element_pos_x = Trans.ElementPos(:, 1)*wvl;
element_pos_z = Trans.ElementPos(:, 3)*wvl;
phi = Trans.ElementPos(:, 4);

for n = 1:n_pulses
    for e = 1:n_elements
        focus = sound_speed*t(:)/2;
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        t_delay(:,e,n) = (focus + sqrt((zfocus- element_pos_z(e)).^2 + ...
            (xfocus - element_pos_x(e)).^2))/sound_speed;
    end
end

end