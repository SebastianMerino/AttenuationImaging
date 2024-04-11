% ====================================================================== %
% Script to fin the acoustic enhancement posterior to a colloid nodule.
% Used for ISBI 2024.
% ====================================================================== %
clear,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\CUELLO#3'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\references'];
croppedDir = [baseDir,'\cropped'];
figDir = 'C:\Users\sebas\Pictures\ISBI2024\v2';

targetFiles = dir([targetDir,'\*.mat']);

blocksize = 15;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 1;
NptodB = log10(exp(1))*20;

% CASOS INTERESANTES: 2,4,6,8,9,13

%% Loading case
iAcq = 8;
load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = RF(:,:,1);

% dynRange = [-50,0];
% BmodeFull = db(hilbert(sam1));
% BmodeFull = BmodeFull - max(BmodeFull(:));
% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
% ylim([0.05 3])

% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
% end
rect = [0.8726    0.1720    2.8224    2.3804];
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

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

dynRange = [-40,-5];
attRange = [0.3,1.7];
%attRange = [0,1]; % Just for 13 acq
bsRange = [-2 2];



%% ACOUSTIC ENHANCEMENT

fs = 40000000;
[pxx,fpxx] = pwelch(sam1,300,250,512,fs);

figure('Units','centimeters', 'Position',[5 5 18 6]),
tiledlayout(1,2)
nexttile, imagesc(x,z,Bmode,dynRange)
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')

nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on


%%
fc = 5E6;
[bFilt,aFilt] = butter(1,[fc-0.5E6 fc+0.5E6]/fs*2, "bandpass");
samFilt = filtfilt(bFilt,aFilt,sam1);
[pxx,fpxx] = pwelch(samFilt,300,250,512,fs);

BmodeFilt = db(hilbert(samFilt));
BmodeFilt = BmodeFilt - max(BmodeFilt(:));

figure('Units','centimeters', 'Position',[5 5 18 6]),
tiledlayout(1,2)
nexttile, imagesc(x,z,BmodeFilt,[-50 0])
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')
xline(1.4, 'g--', 'LineWidth',2)
xline(2.1, 'g--', 'LineWidth',2)
xline(2.3, 'r--', 'LineWidth',2)


nexttile,plot(fpxx/1e6,mean(pxx,2))
title('Spectrum')
xlabel('f [MHz]')
grid on

%%
[~,Z] = meshgrid(x,z);
mask = Z>1.5 & Z<2.5;
%latProfile = sum(BmodeFilt.*mask)./sum(mask);
BmodeFilt(~mask) = NaN;
latProfile = median(BmodeFilt,"omitmissing");
figure('Units','centimeters', 'Position',[5 5 9 6]),
plot(x,latProfile)
grid on
axis tight
xlabel('Lateral [cm]')
title('Acoustic enhancement')
%%
underInclusion = mean(latProfile(x>1.4 & x <2.1))
ousideInclusion = mean(latProfile(x>2.3 & x <3.2))
acEnhancement = underInclusion - ousideInclusion
% Inclusion in 0.85 cm long
attDiff = acEnhancement/1.7/fc*1E6
