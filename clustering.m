clear,clc
close all
addpath('./functions_att');

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% Loading sample
iAcq = 2;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);


%% Segmentation by clustering
figure('Units','centimeters', 'Position',[5 5 30 7]), 
tiledlayout(1,4)
nexttile,
imagesc(x,z,Bmode, dynRange)
colormap gray
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Bmode')

h = fspecial("average",[70 7]);
blurred = imfilter(Bmode,h,"symmetric"); 
nexttile, imagesc(x,z,blurred,dynRange);
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Blurred image')

[X,Z] = meshgrid(x,z);
Data = normalize([X(:) Z(:) blurred(:)]);
Data = Data.*[1.2 1 1.5];
[mBm,nBm] = size(Bmode);
nClusters = 3; % number of clusters
idx = kmeans(Data,nClusters);
ID = reshape(idx,[mBm,nBm]);

nexttile;
segmented = labeloverlay(reshape(normalize(Bmode(:),'range'),[mBm,nBm]),ID);
imagesc(x,z,segmented)
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Segmented image')

% Equalizing image
mask = ID==1;
energy1 = std(sam1(mask));
factor = ones(size(sam1));

for iCluster = 1:nClusters
    mask = ID==iCluster;
    factor(mask) = energy1/std(sam1(mask));
end
h = fspecial("average",[50 5]);
factor = imfilter(factor,h,"symmetric");

% figure('Units','centimeters', 'Position',[5 5 30 8]), 
% tiledlayout(1,2)
% t2 = nexttile;
% imagesc(x,z,factor)
% colormap(t2,parula)
% colorbar
% axis image
% title('Factor')

samEnhanced = sam1.*factor;
Bmode2 = db(hilbert(samEnhanced));
Bmode2 = Bmode2 - max(Bmode2(:));
t3 = nexttile;
imagesc(x,z,Bmode2,dynRange)
colormap(t3,gray)
colorbar
axis image
title('Equalized B-mode')

saveas(gcf,[figDir,'\segmentedT',num2str(iAcq),'.png'])

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

% jj = round(n/2);
% ii = round(m/4);
% figure('Units','centimeters', 'Position',[5 5 10 10])
% plot(f,10*log10(squeeze(Sp(ii,jj,:)/max(Sd(ii,jj,:)))),'k');
% hold on
% plot(f,10*log10(squeeze(Sd(ii,jj,:)/max(Sd(ii,jj,:)))),'r');
% hold off
% title('SLD at interface');
% xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
% axis([f(1) f(end) -30 10]);

%% Au = b
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = 1e4*[0.4,1.2,3.6];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

%% Plotting
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
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
saveas(gcf,[figDir,'\equalizedT',num2str(iAcq),'.png'])
