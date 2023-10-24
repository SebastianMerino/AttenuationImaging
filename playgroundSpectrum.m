clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\data'];
refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\ref'];

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];

% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']); 
figDir = [targetDir,'\fig\24-10-Playground'];
if (~exist(figDir,"dir")), mkdir(figDir); end
%%
for iAcq = 1:length(croppedFiles)
    %iAcq = 2;
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load([refDir,'\compensation.mat']);
attRange = [0.4,1.1];
bsRange = [-2 2];
%% Spectrum
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);
Sp = zeros(m,n,length(f));
Sd = zeros(m,n,length(f));
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% 
% %% Examining spectrum (vertical slices)
% figure,
% for ix=1:n
%     SLDslice = squeeze(b(:,ix,:));
%     imagesc(f,z_ACS,SLDslice, [-3,3])
%     colorbar
%     xlabel('Frequency [Hz]')
%     ylabel('Depth [cm]')
%     title(num2str(ix))
%     pause(0.1)
% end
% 
% %% Examining spectrum (frequency slices)
% figure,
% for iFreq=1:p
%     SLDslice = squeeze(b(:,:,iFreq));
%     im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
%     im.AlphaData = (maskBack | maskInc) + 0.5;
% 
%     % SLDslice = squeeze(log(Sd(:,:,iFreq)));
%     % imagesc(x_ACS,z_ACS,SLDslice, [11,20])
% 
%     colorbar
%     xlabel('x [cm]')
%     ylabel('Depth [cm]')
%     title(['Frequency=',num2str(f(iFreq),2),'MHz'])
%     pause(0.1)
% end

% %% Center of image
% ix = floor(n/2);
% figure,
% for iz = 1:m
%     plot(f,squeeze(log(Sp(iz,ix,:))))
%     hold on
%     plot(f,squeeze(log(Sd(iz,ix,:))))
%     hold off
%     ylim([11 20])
%     legend('Proximal','Distal')
%     xlim([f(1) f(end)])
%     grid on
%     pause(0.1)
% end
% 
% %% Spectral log difference
% figure,
% for iz = 1:m
%     plot(f,squeeze(b(iz,ix,:)))   
%     ylim([-4 4])
%     %legend('Proximal','Distal')
%     xlim([f(1) f(end)])
%     grid on
%     pause(0.1)
% end
%% Masks
%cx = 1.85; cz = 1.9; r = 0.6; rB = 1.2; % Both
%cx = 1.85; cz = 1.9; r = 0.6; rB = 1.1; % T8
%cx = 1.85; cz = 1.9; r = 0.7; rB = 1.2; % T7

% Simulation
cx = 0; cz = 2.75; r = 0.8; rB = 1.2; 
[X,Z] = meshgrid(x_ACS,z_ACS);
maskInc = (X - cx).^2 + (Z - cz).^2 < r^2;
maskBack = (X - cx).^2 + (Z - cz).^2 > rB^2;
im = imagesc(x_ACS,z_ACS,mean(b,3), [-2,2]);
im.AlphaData = (maskBack | maskInc) + 0.5;
colorbar
xlabel('x [cm]')
ylabel('Depth [cm]')
axis image
title('SLD, frequency mean')


%% Spectrums in Inclusion and background
spectrumBack = zeros(size(f));
spectrumInc = zeros(size(f));
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    spectrumInc(iFreq) = mean(SLDslice(maskInc));
    spectrumBack(iFreq) = mean(SLDslice(maskBack));
end

% Fitting lines
slopeInc = f\spectrumInc;
slopeBack = f\spectrumBack;
attIncBsc0(iAcq) = slopeInc*8.686;
attBackBsc0(iAcq) = slopeBack*8.686;

% Fitting lines
fitInc = [ones(size(f)) f]\spectrumInc;
fitBack = [ones(size(f)) f]\spectrumBack;
attInc(iAcq) = fitInc(2)*8.686;
attBack(iAcq) = fitBack(2)*8.686;

figure('Units','centimeters', 'Position',[5 5 25 10]) 
tl = tiledlayout(1,2);
title(tl,'Mean SLD at homogeneous regions')

nexttile,
plot(f,spectrumBack)
hold on
plot(f,spectrumInc)
plot(f,slopeBack*f, 'k--')
plot(f,slopeInc*f, 'k--')
xlim([0 f(end)])
ylim([0 max(spectrumInc)])
grid on
legend('Back','Inc', 'Location','northwest')

nexttile,
plot(f,spectrumBack)
hold on
plot(f,spectrumInc)
plot(f,fitBack(1) + fitBack(2)*f, 'k--')
plot(f,fitInc(1) + fitInc(2)*f, 'k--')
xlim([0 f(end)])
ylim([0 max(spectrumInc)])
grid on
legend('Back','Inc', 'Location','northwest')

%% Masks
%cx = 1.85; cz = 1.9; r = 0.6; rB = 1.2; % Both
% Simulation
cx = 0; cz = 2.75; r = 0.8; rB = 1.2; 
[X,Z] = meshgrid(x_ACS,z_ACS);
maskBorder = (X - cx).^2 + (Z - cz).^2 > r^2 &...
    (X - cx).^2 + (Z - cz).^2 < rB^2; 
maskUp = maskBorder & Z > cz + 0.3;
maskDown = maskBorder & Z < cz - 0.3;

% Examining spectrum (frequency slices)
figure,
im = imagesc(x_ACS,z_ACS,mean(b,3), [-2,2]);
im.AlphaData = (maskUp | maskDown) + 0.5;
colorbar
xlabel('x [cm]')
ylabel('Depth [cm]')
title('SLD, frequency mean')
axis image

%% Spectrums in borders
spectrumUp = zeros(size(f));
spectrumDown = zeros(size(f));
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    spectrumUp(iFreq) = mean(SLDslice(maskUp));
    spectrumDown(iFreq) = mean(SLDslice(maskDown));
end

% Fitting lines
slopeUp = f\spectrumUp;
slopeDown = f\spectrumDown;
% attUp = slopeUp*8.686
% attDown = slopeDown*8.686

figure('Units','centimeters', 'Position',[5 5 25 10]) 
tl = tiledlayout(1,2);
title(tl,'Mean SLD at borders')
nexttile,
plot(f,spectrumDown)
hold on
plot(f,spectrumUp)
plot(f,slopeDown*f, 'k--')
plot(f,slopeUp*f, 'k--')

xlim([0 f(end)])
%ylim([0 max(spectrumUp)])
grid on
legend('Down','Up', 'Location','northwest')
title('Considering same BSC')

% Fitting lines
fitUp = [ones(size(f)) f]\spectrumUp;
fitDown = [ones(size(f)) f]\spectrumDown;
% attUp = fitUp(2)*8.686
% attDown = fitDown(2)*8.686

nexttile,
plot(f,spectrumDown)
hold on
plot(f,spectrumUp)
plot(f,fitDown(1) + fitDown(2)*f, 'k--')
plot(f,fitUp(1) + fitUp(2)*f, 'k--')

xlim([0 f(end)])
%ylim([0 max(spectrumUp)])
grid on
legend('Down','Up', 'Location','northwest')
title('Considering different BSC')

%%
newDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(newDir,"dir")), mkdir(newDir); end
save_all_figures_to_directory(newDir);
close all

end

%% Error
close all
%groundTruthT = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];
groundTruthBack = [0.5 0.5 0.5 0.6 0.6 0.6];
groundTruthInc = [1 1 1 1.4 1.4 1.4];

% PEInc = ( attInc - groundTruthT(1:8) )./groundTruthT(1:8)*100;
% PEBack = ( attBack - groundTruthT(9) )./groundTruthT(9)*100;
% PEIncBsc0 = ( attIncBsc0 - groundTruthT(1:8) )./groundTruthT(1:8)*100;
% PEBackBsc0 = ( attBackBsc0 - groundTruthT(9) )./groundTruthT(9)*100;
PEInc = ( attInc - groundTruthInc )./groundTruthInc*100;
PEBack = ( attBack - groundTruthBack )./groundTruthBack*100;
PEIncBsc0 = ( attIncBsc0 - groundTruthInc )./groundTruthInc*100;
PEBackBsc0 = ( attBackBsc0 - groundTruthBack )./groundTruthBack*100;
figure('Units','centimeters', 'Position',[5 5 12 9])
plot(PEInc,'o', 'LineWidth',2)
hold on
plot(PEIncBsc0,'x', 'LineWidth',2)
hold off
yline(0, 'k--', 'LineWidth',2)
legend('\DeltaBSC\neq0','\DeltaBSC=0')
grid on
title('Percentage error in Inclusion')

figure('Units','centimeters', 'Position',[5 5 12 9])
plot(PEBack,'o', 'LineWidth',2)
hold on
plot(PEBackBsc0,'x', 'LineWidth',2)
hold off
yline(0, 'k--', 'LineWidth',2)
legend('\DeltaBSC\neq0','\DeltaBSC=0')
grid on
title('Percentage error in Background')

figure('Units','centimeters', 'Position',[5 5 12 9])
plot(attInc - attBack,'o', 'LineWidth',2)
hold on
plot(attIncBsc0 - attBackBsc0,'x', 'LineWidth',2)
%plot(groundTruthT(1:8) - groundTruthT(9),'k_', 'LineWidth',2)
plot(groundTruthInc - groundTruthBack,'k_', 'LineWidth',2)
hold off
grid on
title('Difference in attenuation')

save_all_figures_to_directory(figDir);
