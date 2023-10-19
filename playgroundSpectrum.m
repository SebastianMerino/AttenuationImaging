clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID316V2\06-08-2023-Generic'];

% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig\18-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end
%%
iAcq = 8;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
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


%% Examining spectrum (vertical slices)
figure,
for ix=1:n
    SLDslice = squeeze(b(:,ix,:));
    imagesc(f,z_ACS,SLDslice, [-3,3])
    colorbar
    xlabel('Frequency [Hz]')
    ylabel('Depth [cm]')
    title(num2str(ix))
    pause(0.1)
end

%% Examining spectrum (frequency slices)
figure,
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
    %im.AlphaData = (maskBack | maskInc) + 0.5;

    %SLDslice = squeeze(log(Sd(:,:,iFreq)));
    %imagesc(x_ACS,z_ACS,SLDslice, [11,20])
    colorbar
    xlabel('x [cm]')
    ylabel('Depth [cm]')
    title(['Frequency=',num2str(f(iFreq),2),'MHz'])
    pause(0.1)
end

%% Masks
cx = 1.85; cz = 1.9; r = 0.6; rB = 1.2; % Both
%cx = 1.85; cz = 1.9; r = 0.6; rB = 1.1; % T8
%cx = 1.85; cz = 1.9; r = 0.7; rB = 1.2; % T7
[X,Z] = meshgrid(x_ACS,z_ACS);
maskInc = (X - cx).^2 + (Z - cz).^2 < r^2;
maskBack = (X - cx).^2 + (Z - cz).^2 > rB^2;
im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
im.AlphaData = (maskBack | maskInc) + 0.5;

%% Spectrums in Inclusion and background
spectrumBack = zeros(size(f));
spectrumInc = zeros(size(f));
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    spectrumInc(iFreq) = mean(SLDslice(maskInc));
    spectrumBack(iFreq) = mean(SLDslice(maskBack));
end

%% Fitting lines
disp('Forcing no BSC difference')
slopeInc = f\spectrumInc;
slopeBack = f\spectrumBack;
attInc = slopeInc*8.686
attBack = slopeBack*8.686

figure('Units','centimeters', 'Position',[5 5 25 10]) 
tiledlayout(1,2)
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

% Fitting lines
disp('Considering BSC difference')
fitInc = [ones(size(f)) f]\spectrumInc;
fitBack = [ones(size(f)) f]\spectrumBack;
attInc = fitInc(2)*8.686
attBack = fitBack(2)*8.686

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
cx = 1.85; cz = 1.9; r = 0.6; rB = 1.2; % Both
[X,Z] = meshgrid(x_ACS,z_ACS);
maskBorder = (X - cx).^2 + (Z - cz).^2 > 0.8^2 &...
    (X - cx).^2 + (Z - cz).^2 < 1^2; 
maskUp = maskBorder & Z > cz + 0.2;
maskDown = maskBorder & Z < cz - 0.2;

% Examining spectrum (frequency slices)
figure,
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
    im.AlphaData = (maskUp | maskDown) + 0.5;

    %SLDslice = squeeze(log(Sd(:,:,iFreq)));
    %imagesc(x_ACS,z_ACS,SLDslice, [11,20])
    colorbar
    xlabel('x [cm]')
    ylabel('Depth [cm]')
    title(['Frequency=',num2str(f(iFreq),2),'MHz'])
    pause(0.1)
end

%% Spectrums in borders
spectrumUp = zeros(size(f));
spectrumDown = zeros(size(f));
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    spectrumUp(iFreq) = mean(SLDslice(maskUp));
    spectrumDown(iFreq) = mean(SLDslice(maskDown));
end

% Fitting lines
disp('Forcing no BSC difference')
slopeUp = f\spectrumUp;
slopeDown = f\spectrumDown;
attUp = slopeUp*8.686
attDown = slopeDown*8.686

figure('Units','centimeters', 'Position',[5 5 25 10]) 
tiledlayout(1,2)
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

% Fitting lines
disp('Considering BSC difference')
fitUp = [ones(size(f)) f]\spectrumUp;
fitDown = [ones(size(f)) f]\spectrumDown;
attUp = fitUp(2)*8.686
attDown = fitDown(2)*8.686

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
