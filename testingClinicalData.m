%clear,clc
close all
addpath('./functions_v7');

% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\ThyroidSelected\CUELLO#3'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\CUELLO#2'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

figDir = [baseDir,'\fig\28-11'];
if (~exist(figDir,"dir")), mkdir(figDir); end

% CASOS INTERESANTES CUELLO 3:  2,4,6,8,9,13

% CASOS INTERESANTES CUELLO 2:  1,3,9

%% Loading data
for iAcq = 1:length(croppedFiles)

iAcq =1;
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
dynRange = [-40,-5];
attRange = [0.3,1.7];
%attRange = [0,1]; % Just for 13 acq
bsRange = [-2 2];

%% System
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
b = (log(Sp) - log(Sd)) - (compensation);

tol = 1e-3;
clear mask
mask = ones(m,n,p);

%% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
[u,~] = cgs(A'*A,A'*b(:),1e-6,20);
BS = u(1:end/2); CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);
CS = reshape(CS,m,n);

figure('Units','centimeters', 'Position',[5 5 30 8]);
tl = tiledlayout(1,3);
title(tl,'Standard RSLD')
subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
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


%% RSLD

mu = logspace(2,3,3);
mu2 = logspace(-1,1,3);
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
    subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
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

%%
tic
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),1E3,10,m,n,tol,mask(:));
toc
BRTV = (reshape(Bn*8.686,m,n));



%% British Columbia Approach
envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
        
        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% figure,
% imagesc(x_ACS,z_ACS,SNR)
% colorbar

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2,3,3);
mu2 = logspace(-1,1,3);
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
    subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
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
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2,3,3);
mu2 = logspace(-1.5,0.5,3);
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
    title(tl,{'RSLD with isotropic TV and Tikhonov reg.',''})
    subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
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
    
    nexttile; axis off
    
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

%%
tic
% CASO 1
% [Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),1E3,10^0.5,m,n,tol,mask(:));

% CASO 3
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),1E3,10^-0.5,m,n,tol,mask(:));
toc
BRTVL1 = (reshape(Bn*8.686,m,n));


%% NEW WEIGHTS
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
NptodB = 20*log10(exp(1));

% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);

mu = 1e3;
mu2 = 10;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
% Weight function
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;


figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'New Weights',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
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
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
% W = speye(m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;
%%
% Regularization: Au = b
tol = 1e-3;

mask = ones(m,n,p);
mu = 10.^(3:0.5:4);
mu2 = 10.^(-0.5:1:1.5);
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
    title(tl,{'TV, Tikhonov reg and weights',''})
    subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
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
tic
% CASO 1
% [Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,10^3.5,10^0.5,m,n,tol,mask(:),w);
% CASO 3
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,2*10^3,10^-0.5,m,n,tol,mask(:),w);
toc
BRWFR = (reshape(Bn*8.686,m,n));


figure('Units','centimeters', 'Position',[5 5 30 8]);
tl = tiledlayout(1,4);
title(tl,{'Comparison'})
subtitle(tl,['Patient ',croppedFiles(iAcq).name(1:end-4)])
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t2,turbo)
axis equal tight
% title(['TV, \mu=',num2str(mu(mm),2)])
title('TV')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t2,turbo)
axis equal tight
%title(['SWTV, \mu=',num2str(mu(mm),2)])
title('TV-L1')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
%title(['WFR, \mu=',num2str(mu(mm),2)])
title('WFR')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%%
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure,
[~,~,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';

ylim([0.1, 3])
%%

[X,Z] = meshgrid(x_ACS,z_ACS);
% inc = X>1.4 & Z<1.6 & Z>0.8;
inc = X>1.8 & Z<2.2 & Z>1.3;

mean(BRTV(inc),'all')
mean(BRTVL1(inc),'all')
mean(BRWFR(inc),'all')
im = imagesc(BRWFR);
im.AlphaData = (0.5 + inc);
%%
newDir = fullfile(figDir,croppedFiles(iAcq).name(1:end-4));
if(~exist(newDir,"dir")), mkdir(newDir); end
save_all_figures_to_directory(newDir,'figure');
close all


end