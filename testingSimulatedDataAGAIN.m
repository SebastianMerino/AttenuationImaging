clear,clc
addpath('./functions_v7');
addpath('./AttUtils');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered_14_11_23'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 15;     % Block size in wavelengths
freq_L = 4e6; freq_H = 9e6;
%freq_L = 3e6; freq_H = 9e6;

overlap_pc      = 0.8;
ratio_zx        = 1;

targetFiles = dir([targetDir,'\*.mat']);


%% For looping
iAcq = 2;
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
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');

confirmation = '';
while ~strcmp(confirmation,'Yes')
    rect = getrect;
    confirmation = questdlg('Sure?');
end
close,

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



%% RSLD and plotting
dynRange = [-50,0];
attRange = [0.4,1.4];
bsRange = [-2 2];

%% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );

% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3)*10;
muC = logspace(-0.5,1.5,3)*10;
BRTV = zeros(m,n,length(muC));
CRTV = zeros(m,n,length(muC));
for mm = 1:length(mu)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),muC(mmC),m,n,tol,mask(:));
        toc
        BRTV(:,:,mmC) = (reshape(Bn*8.686,m,n));
        CRTV(:,:,mmC) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BRTV,3)+1);
    title(tl,{'RSLD - TV',''})
    %subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BRTV,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BRTV(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    nexttile;
    axis off
    
    for ii = 1:size(BRTV,3)
        t4 = nexttile; 
        imagesc(x_ACS,z_ACS,CRTV(:,:,ii), bsRange)
        colormap(t4,parula)
        axis image
        title(['RSLD, \mu=',num2str(muC(ii),2)])
    end
    c = colorbar(t4);
    c.Label.String = 'BS log ratio (a.u.)';
end
%%

figure('Units','centimeters', 'Position',[5 5 6 8]),
imOverlayInterp(BmodeFull,BRTV(:,:,2),dynRange,attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('RSLD-TV')
xlabel('x [cm]'), ylabel('z [cm]'), 

fprintf('ACS = %.2f\n',mean(BRTV(:,:,2),'all'))
%% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );

mu = logspace(2.5,3.5,3)*10;
muC = logspace(-0.5,1.5,3)*10;
BWFR = zeros(m,n,length(muC));
CWFR = zeros(m,n,length(muC));
for mm = 1:length(mu)
    for mmC = 1:length(muC)
        tic
        [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu(mm),muC(mmC),m,n,tol,mask(:));
        bscMap = (reshape(Cn,m,n));

        logBscRatio = bscMap*log10(exp(1))*20;
        w = 1./((logBscRatio/6).^4 + 1);

        W = repmat(w,[1 1 p]);
        W = spdiags(W(:),0,m*n*p,m*n*p);
        bw = W*b(:);        
        A1w = W*A1;
        A2w = W*A2;

        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu(mm),muC(mmC),m,n,tol,mask(:),w);
        toc
        BWFR(:,:,mmC) = (reshape(Bn*8.686,m,n));
        CWFR(:,:,mmC) = (reshape(Cn,m,n));
    end

    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BWFR,3)+1);
    title(tl,{'TV, Tikhonov reg and weights',''})
    %subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')

    for ii = 1:size(BWFR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BWFR(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';

    t3 = nexttile;
    imagesc(x_ACS,z_ACS,w,[0 1])
    axis image
    colormap(t3,parula)
    colorbar(t3,'westoutside')
    title('Weights')

    for ii = 1:size(BWFR,3)
        t4 = nexttile; 
        imagesc(x_ACS,z_ACS,CWFR(:,:,ii), bsRange)
        colormap(t4,parula)
        axis image
        title(['RSLD, \mu=',num2str(muC(ii),2)])
    end
    c = colorbar(t4);
    c.Label.String = 'BS log ratio (a.u.)';
    pause()
end
%%
figure('Units','centimeters', 'Position',[5 5 6 8]),
imOverlayInterp(BmodeFull,BWFR(:,:,3),dynRange,attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('RSLD-WFR')
xlabel('x [cm]'), ylabel('z [cm]'), 

fprintf('ACS = %.2f\n',mean(BWFR(:,:,3),'all'))

%% Saving data
% save(fullfile(croppedDir,targetFiles(iAcq).name),"Sd","Sp",...
%     "compensation","z_ACS","x_ACS","nx","nz","x0","z0p","z0d","sam1",...
%     "m","n","p","Bmode","x","z","f","L","xFull","zFull","BmodeFull")

% end



