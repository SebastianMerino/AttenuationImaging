clear,clc
addpath('./functions_v7');
addpath('./AttUtils');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered_14_11_23'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];
figDir = [baseDir,'\fig\17-11'];

if ~exist(figDir,"dir"), mkdir(figDir); end
%% Generating cropped data
% SETTING PARAMETERS
%for blocksize = 10:5:20     % Block size in wavelengths
blocksize = 15;
freq_L = 4e6; freq_H = 9e6;

%freq_L = 3e6; freq_H = 9e6;
overlap_pc      = 0.8;
ratio_zx        = 1;
NptodB = log10(exp(1))*20;
targetFiles = dir([targetDir,'\*.mat']);

%% For looping
iAcq = 5;
load(fullfile(targetDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = rf(:,:,1);
dynRange = [-50,0];

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));
% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');

% confirmation = '';
% while ~strcmp(confirmation,'Yes')
%     rect = getrect;
%     confirmation = questdlg('Sure?');
% end
% close,

%% Cropping and finding sample sizes
% Region for attenuation imaging
% x_inf = rect(1); x_sup = rect(1)+rect(3);
% z_inf = rect(2); z_sup = rect(2)+rect(4);
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.5; z_sup = 3.5;

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
nz = 2*round(blocksize*wl/dz /2); % Window size
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
dynRange = [-40 -10];
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
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = 0.7*(f.^1.05)/8.6858;
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
refFiles = dir([refDir,'\rf*.mat']);
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load([refDir,'\',refFiles(iRef).name]);
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


% ====================================================================== %
%% RSLD and plotting
dynRange = [-50,0];
attRange = [0.4,1.4];
bsRange = [-15 15];

groundTruthTop = [0.6,0.6,0.6,1.2,1.2,1.2];
groundTruthBottom = [1.2,1.2,1.2,0.6,0.6,0.6];
[~,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
attIdeal(Z<=2) = groundTruthTop(iAcq);
attIdeal(Z>2) = groundTruthBottom(iAcq);

%% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );

% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
% muB = 10.^(3:0.5:4);
% muC = 10.^(1:0.5:2);
muB = 10^3.5; muC = 10^1.5;
BRTV = zeros(m,n,length(muC));
CRTV = zeros(m,n,length(muC));
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BRTV(:,:,mmC) = (reshape(Bn*NptodB,m,n));
        CRTV(:,:,mmC) = (reshape(Cn*NptodB,m,n));
    end
        
    figure('Units','centimeters', 'Position',[5 5 8 12]);
    tl = tiledlayout(2,1);
    title(tl,['RSLD-TV, BS = ',num2str(blocksize),'\lambda'])    
    
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BRTV(:,:,1), attRange)
    colormap(t2,turbo)
    axis image
    title(['RSLD, \mu_B=',num2str(muB(mmB),2)])
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    t3 = nexttile; 
    imagesc(x_ACS,z_ACS,CRTV(:,:,1), bsRange)
    colormap(t3,parula)
    axis image
    title(['RSLD, \mu_C=',num2str(muC(mmC),2)])
    c = colorbar;
    c.Label.String = 'BS log ratio [dB]';
end

% figure('Units','centimeters', 'Position',[5 5 6 8]),
% imOverlayInterp(BmodeFull,BRTV(:,:,2),dynRange,attRange,0.5,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('RSLD-TV')
% xlabel('x [cm]'), ylabel('z [cm]'), 
% 
% fprintf('ACS = %.2f\n',mean(BRTV(:,:,2),'all'))

%% TV on Att an L1 norm on BSC
% muB = 10.^(2.5:0.5:3.5);
% muC = 10.^(-0.5:0.5:0.5);
muB = 10^3.5; muC = 10^0.5;
BWFR = zeros(m,n,length(muC));
CWFR = zeros(m,n,length(muC));
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BWFR(:,:,mmB) = reshape(Bn*NptodB,m,n);
        CWFR(:,:,mmC) = reshape(Cn*NptodB,m,n);
        % Plotting
        figure('Units','centimeters', 'Position',[5 5 8 12]);
        tl = tiledlayout(2,1);
        title(tl,['RSLD-TVL1, BS = ',num2str(blocksize),'\lambda'])

        t2 = nexttile;
        imagesc(x_ACS,z_ACS,BWFR(:,:,mmB), attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD, \mu_B=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';

        t3 = nexttile;
        imagesc(x_ACS,z_ACS,CWFR(:,:,mmC), bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD, \mu_C=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio [dB]';

    end

end

% figure('Units','centimeters', 'Position',[5 5 6 8]),
% imOverlayInterp(BmodeFull,BWFR(:,:,2),dynRange,attRange,0.5,...
%     x_ACS,z_ACS,roi,xFull,zFull);
% title('RSLD-TV')
% xlabel('x [cm]'), ylabel('z [cm]'), 
% 
% fprintf('ACS = %.2f\n',mean(BWFR(:,:,2),'all'))


%% Weights
% order = 1;
% 
% [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),10^3.5,10^0.5,m,n,tol,mask(:));
% bscMap = (reshape(Cn,m,n));
% 
% logBscRatio = bscMap*log10(exp(1))*20;
% w = 1./((logBscRatio/6).^(2*order) + 1);
% %w = ones(size(BRTV(:,:,1)));
% figure('Units','centimeters','Position',[5 5 10 10]),
% imagesc(x_ACS,z_ACS,w,[0 1])
% axis image
% colormap(parula)
% colorbar
% title(['Weights, order =',num2str(order)])
% 
% w = movmin(w,5);
% w = medfilt2(w,[5 5],'symmetric');
% figure('Units','centimeters','Position',[5 5 10 10]),
% imagesc(x_ACS,z_ACS,w,[0 1])
% axis image
% colormap(parula)
% colorbar
% title('New weights')
% 
% W = repmat(w,[1 1 p]);
% W = spdiags(W(:),0,m*n*p,m*n*p);
% bw = W*b(:);        
% A1w = W*A1;
% A2w = W*A2;
ratioCutOff = 10;
xWeights = -20:0.05:20;
legends = {};
figure('Units','centimeters', 'Position',[5 5 15 8])
for order = 1:5
    yWeights = 1./((xWeights/ratioCutOff).^(2*order) + 1);
    plot(xWeights,yWeights)
    legends{order} = ['order=',num2str(order)];
    hold on
end
xline(ratioCutOff,'--')
xline(-ratioCutOff,'--')
hold off
legend(legends)
grid on
xlabel('Log BSC ratio [dB]')
ylabel('Weights')
title(['Weights, cut-off = ',num2str(ratioCutOff),'dB'])

%% Weighting equation and regularizations
% muB = 10.^(3:0.5:4);
% muC = 10.^(1:0.5:2);
ratioCutOff = 10;
order = 2;
extension = 1;
for extension = 1:2:7
muB = 10^4; muC = 10^0.5;
BWFR = zeros(m,n,length(muC));
CWFR = zeros(m,n,length(muC));

for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));        
        logBscRatio = reshape(Cn*NptodB,m,n);
        w = (1./((logBscRatio/ratioCutOff).^(2*order) + 1));
        w = movmin(w,extension);

        W = repmat(w,[1 1 p]);
        W = spdiags(W(:),0,m*n*p,m*n*p);
        bw = W*b(:);        
        A1w = W*A1;
        A2w = W*A2;

        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        toc

        BWFR(:,:,mmC) = reshape(Bn*NptodB,m,n);
        CWFR(:,:,mmC) = reshape(Cn*NptodB,m,n);
    end
    fprintf('ACS Bottom: %.2f\n',mean(BWFR(22:end,:),'all'))

    % Plotting
    figure('Units','centimeters', 'Position',[5 5 8 12]);
    tl = tiledlayout(2,1);
    title(tl,['RSLD-WFR, BS = ',num2str(blocksize),'\lambda'])    
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BWFR(:,:,1), attRange)
    colormap(t2,turbo)
    axis image
    title(['RSLD, \mu_B=',num2str(muB(mmB),2)])
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    t3 = nexttile; 
    imagesc(x_ACS,z_ACS,w,[0 1])
    colormap(t3,parula)
    colorbar
    %title(['Weights, CutOff=',num2str(ratioCutOff),'dB'])
    %title(['Weights, order=',num2str(order)])
    title(['Weights, extension=',num2str(extension)])
    axis image
end

end
% save_all_figures_to_directory(figDir,'extension');
% close all;

%% MANUAL WEIGHTS
[~,border] = gradient(attIdeal);
manualW = double(border==0);

manualW(manualW == 0) = 10^-2;
extension = 5;

dxACS = x_ACS(2)-x_ACS(1);
BScm = blocksize*wl*100;
manualW = movmin(manualW,extension);
fprintf("\nBS = %.2f cm, Border width = %.2f cm\n",BScm,dxACS*(extension+1))

w = manualW;
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;


muB = 10^4;
muC = 10.^(0.5:1:2.5);
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        toc

        BWFR(:,:,mmC) = reshape(Bn*NptodB,m,n);
        CWFR(:,:,mmC) = reshape(Cn*NptodB,m,n);
        fprintf('ACS Bottom: %.2f\n',mean(BWFR(22:end,:,mmC),'all'))
    
        % Plotting
        figure('Units','centimeters', 'Position',[5 5 20 6]);
        tl = tiledlayout(1,3);
        t3 = nexttile; 
        imagesc(x_ACS,z_ACS,w,[0 1])
        colormap(t3,parula)
        colorbar
        title('Manual weights')
        axis image
    
        title(tl,['RSLD-WFR, BS = ',num2str(blocksize),'\lambda'])    
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BWFR(:,:,mmC), attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD, \mu_B=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';
        
        t3 = nexttile;
        imagesc(x_ACS,z_ACS,CWFR(:,:,mmC), bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD, \mu_C=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio [dB]';
    end


end

save_all_figures_to_directory(figDir,'manualR100');
close all;

%% Saving data
% save_all_figures_to_directory(figDir,['L',num2str(blocksize),'figure']);
% close all;


% end


