clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered'];
% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\ref'];

% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets' ...
%     '\Attenuation\Timana\data'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\Timana\ref'];

croppedDir = [targetDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']); 
figDir = [targetDir,'\fig\24-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% Loading data
for iAcq = 1:length(croppedFiles)
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load([targetDir,'\compensation.mat']);
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

%% Standard SLD
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

[u,~] = cgs(A'*A,A'*b(:),1e-6,20);

% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);
CS = reshape(CS,m,n);

figure('Units','centimeters', 'Position',[5 5 30 8]);
tl = tiledlayout(1,3);
title(tl,'Standard RSLD')
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BS, attRange)
colormap(t2,turbo)
axis image
title('SLD')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CS, bsRange)
colormap(t3,parula)
axis image
title('SLD')
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';
%% Diffraction compensation
diffraction_zf = squeeze(mean(diffraction_compensation,2));
SLD = squeeze(mean((log(Sp) - log(Sd)),2));
SLDcomp = squeeze(mean(b,2));

figure('Units','centimeters', 'Position',[5 5 25 12]);
tiledlayout(1,3)
nexttile,
imagesc(f,z_ACS,SLD, [-2 2]);
title('SLD'),
xlabel('f [cm]'), ylabel('z [cm]'),
colorbar
nexttile,
imagesc(f,z_ACS,diffraction_zf, [-2 2]);
title('Diffraction compensation'),
xlabel('f [MHz]'), ylabel('z [cm]'),
colorbar

nexttile,
imagesc(f,z_ACS,SLDcomp, [-2 2]);
title('SLD compensated'),
xlabel('f [MHz]'), ylabel('z [cm]'),
colorbar


%% RSLD
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0,1,3);
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'Isotropic RSLD')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    nexttile;
    axis off
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end
%% JUST FOR LAYERED MEDIA
[X,Z] = meshgrid(x,z);
layerMask = Z>1.5 & Z < 2.95;
samEnhanced = sam1;
factor = std(sam1(~layerMask))/std(sam1(layerMask));
samEnhanced(layerMask) = samEnhanced(layerMask)*factor/2;
BmodeEnhanced = db(hilbert(samEnhanced));
BmodeEnhanced = BmodeEnhanced - max(BmodeEnhanced(:));

% Spectrum and RSLD
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);
Sp2 = zeros(m,n,length(f));
Sd2 = zeros(m,n,length(f));
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = samEnhanced(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = samEnhanced(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        Sp2(ii,jj,:) = (tempSp(rang));
        Sd2(ii,jj,:) = (tempSd(rang));
    end
end

% RSLD
b = (log(Sp2) - log(Sd2)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0,1,3)*10;
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end

    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'Isotropic RSLD with equalized BSC')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')

    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';

    t3 = nexttile;
    imagesc(x,z,BmodeEnhanced,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t3,gray)
    colorbar(t3,'westoutside')
    title('Equalized Bmode')

    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end


%% British Columbia Approach
envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        
        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

figure,
imagesc(x_ACS,z_ACS,SNR)
colorbar

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
a = 1; b = 0.1;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2,3,3);
mu2 = logspace(0,1,3);
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));

for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),mu(mm),mu2(mm2),...
        m,n,tol,mask(:),w);
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'British Columbia Approach')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
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
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end


%% Minimizing BS log ratio
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(-0.5,0.5,3)*10;
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'RSLD with isotropic TV and Tikhonov reg.')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t3,gray)
    colorbar(t3,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end


%% PROPOSED WEIGHTS
envelope = abs(hilbert(sam1));

BSratio = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        BSratio(ii,jj) = mean(sub_block_p(:))/mean(sub_block_d(:));
    end
end
w = 1./((db(BSratio)/3).^2 + 1);

% Weights
figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,'Proposed Weights');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,db(BSratio))
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('BS ratio')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')



%% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tol = 1e-3;

mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(-0.5,0.5,3);
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu(mm),mu2(mm2),m,n,tol,mask(:),w);
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'TV, Tikhonov reg and weights')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
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
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end
%% NEW WEIGHTS
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );


% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
mu = 1e3;
mu2 = 1;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));
bscMap = (reshape(Cn,m,n));

logBscRatio = bscMap*log10(exp(1))*20;
w = 1./((logBscRatio/10).^2 + 1);


figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,'New Weights');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,logBscRatio)
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('BS ratio')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')

%% Weighting equation and regularizations
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
% W = speye(m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tol = 1e-3;

mask = ones(m,n,p);
mu = logspace(3,4,3);
mu2 = logspace(0,1,3);
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu(mm),mu2(mm2),m,n,tol,mask(:),w);
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'TV, Tikhonov reg and weights')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
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
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis image
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar(t2);
    c.Label.String = 'BS log ratio (a.u.)';
end

%%
newDir = fullfile(figDir,croppedFiles(iAcq).name(1:end-4));
if(~exist(newDir,"dir")), mkdir(newDir); end
save_all_figures_to_directory(newDir);
close all
end