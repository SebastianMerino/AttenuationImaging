clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig\18-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end
%%
iAcq = 7;
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
    pause(0.2)
end

%% Examining spectrum (frequency slices)
figure,
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
    im.AlphaData = (maskBack | maskInc) + 0.5;

    %SLDslice = squeeze(log(Sd(:,:,iFreq)));
    %imagesc(x_ACS,z_ACS,SLDslice, [11,20])
    colorbar
    xlabel('x [cm]')
    ylabel('Depth [cm]')
    title(['Frequency=',num2str(f(iFreq),2),'MHz'])
    pause(0.1)
end

%% Masks
%cx = 1.85; cz = 1.9; r = 0.6; rB = 1.1; % T8
cx = 1.85; cz = 1.9; r = 0.7; rB = 1.2; % T7
[X,Z] = meshgrid(x_ACS,z_ACS);
maskInc = (X - cx).^2 + (Z - cz).^2 < r^2;
maskBack = (X - cx).^2 + (Z - cz).^2 > rB^2;
im = imagesc(x_ACS,z_ACS,SLDslice, [-3,3]);
im.AlphaData = (maskBack | maskInc) + 0.5;

%%
spectrumBack = zeros(size(f));
spectrumInc = zeros(size(f));
for iFreq=1:p
    SLDslice = squeeze(b(:,:,iFreq));
    spectrumInc(iFreq) = mean(SLDslice(maskInc));
    spectrumBack(iFreq) = mean(SLDslice(maskBack));
end


figure,
plot(f,spectrumBack)
hold on
plot(f,spectrumInc)
xlim([0 f(end)])
ylim([0 0.7])
grid on
legend('Back','Inc', 'Location','northwest')

%% Fitting lines
slopeInc = f\spectrumInc;
slopeBack = f\spectrumBack;
attInc = slopeInc*8.686
attBack = slopeBack*8.686

figure,
plot(f,spectrumBack)
hold on
plot(f,spectrumInc)
plot(f,slopeBack*f, 'k--')
plot(f,slopeInc*f, 'k--')

xlim([0 f(end)])
ylim([0 0.7])
grid on
legend('Back','Inc', 'Location','northwest')

