clear,clc
close all
addpath('./functions_v7');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\', ...
    'Attenuation\DataQUS_4_Merino'];
% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets' ...
%     '\Attenuation\DataQUS_4_Merino'];
targetDir = [baseDir,'\Hashimoto'];
refDir = [baseDir,'\References\P4-CUELLO-3'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig\18-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end
%% Loading data
for iAcq = 1:10
iAcq = 9;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);
attRange = [0.4,1.6];
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
BS = u(1:end/2); %CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);

figure('Units','centimeters', 'Position',[5 5 20 8]);
tl = tiledlayout(1,2);
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
axis equal tight
title('SLD')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%% RSLD
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = mu/100;
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),mu2(mm),m,n,tol,mask(:));
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 12]);
tl = tiledlayout(2,size(BR,3)+1);
title(tl,'RSLD')
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
    axis image
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile;
imagesc(x,z,BmodeEqual,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t3,gray)
colorbar(t3,'westoutside')
title('Equalized Bmode')

for ii = 1:size(CR,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
    colormap(t2,parula)
    axis image
    title(['RSLD, \mu=',num2str(mu2(ii),2)])
end
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';

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

% Weights
figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'Weights proposed by BC',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,db(SNR))
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')


SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
a = 1; b = 0.1;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')
%title(['Weights, order=',num2str(gamma)])

%% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
mu2 = mu/10;
for mm = 1:length(mu)
    tic
    [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),mu(mm),mu2(mm),...
        m,n,tol,mask(:),w);
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
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
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar(t3,'westoutside')
axis image
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

%% ---------------------- MANUAL COMPENSATION ----------------------
% iAcq = 1;
% load([croppedDir,'\T',num2str(iAcq),'.mat'])
h = fspecial("average",[50 5]);
blurred = imfilter(Bmode,h,"symmetric");
% [BW,~] = segmentImage(blurred);
save(['./masks/T',num2str(iAcq),'.mat'],'BW','BW1');

%% Refining mask 
% Launch the SEGMENTATION TOOL
fig = figure(10);
fig.Units = 'centimeters'; fig.Position = [5 5 30 8]; 

tiledlayout(1,3)
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis image
title('B-mode')

% Equalization
t2 = nexttile;
load(['./masks/T',num2str(iAcq),'.mat']);
mask = BW;
mask1 = BW1;
factor = ones(size(mask));
if sum(~mask(:)) ~=0
    factor(mask) = std(sam1(mask1))/std(sam1(mask));
    h = fspecial("average",[50 5]);
    factor = imfilter(factor,h,"symmetric");
    imagesc(x,z,factor)
else
    imagesc(x,z,factor,[0.5 1.5])
end
colormap(t2,parula)
colorbar
axis image
title('Factor')

samEnhanced = sam1.*factor;
BmodeEqual = db(hilbert(samEnhanced));
BmodeEqual = BmodeEqual - max(BmodeEqual(:));
t3 = nexttile;
imagesc(x,z,BmodeEqual,dynRange)
colormap(t3,gray)
colorbar
axis image
title('Equalized B-mode')

%% Calculating spectra
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

        sub_block_p = samEnhanced(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = samEnhanced(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% Au = b
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];


% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = mu/100;
BRE = zeros(m,n,length(mu));
CRE = zeros(m,n,length(mu));
for mm = 1:length(mu)
    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),mu2(mm),m,n,tol,mask(:));
    toc
    BRE(:,:,mm) = (reshape(Bn*8.686,m,n));
    CRE(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 12]);
tl = tiledlayout(2,size(BR,3)+1);
title(tl,'Equalized RSLD')
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

for ii = 1:size(BRE,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BRE(:,:,ii), attRange)
    colormap(t2,turbo)
    axis image
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile;
imagesc(x,z,BmodeEqual,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t3,gray)
colorbar(t3,'westoutside')
title('Equalized Bmode')

for ii = 1:size(CRE,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,CRE(:,:,ii), bsRange)
    colormap(t2,parula)
    axis image
    title(['RSLD, \mu=',num2str(mu2(ii),2)])
end
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';

%%
targetDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);
close all
end