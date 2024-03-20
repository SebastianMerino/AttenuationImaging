% ====================================================================== %
% Script to explore homogeneous ROIs in clinical data. 
% Created on Mar 14, 2024
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
if (~exist(resultsDir,"dir")), mkdir(resultsDir); end


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
dynRange = [-60,0];
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
% rect1 = [0.4    0.75    0.4    0.8];
% rect2 = [1.6    0.75    0.4    0.8];
rect1 = [0.4    1.1    0.4    0.5];
rect2 = [1.6    1.1    0.4    0.5];

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
%sam1-mean(sam1)
pxx = fft(sam1(:,end),1024);
fpxx = (0:1024-1)/size(pxx,1)*fs;
meanSpectrum = mean(pxx,2);
spec{iRoi} =  db(meanSpectrum);

%%
figure('Units','centimeters', 'Position',[5 5 10 8]);
imagesc(x,z,Bmode); axis image; colormap gray; %clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
%ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('B-mode')
fontsize(gcf,8,'points')



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