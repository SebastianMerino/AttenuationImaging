clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% NOTA: LAS MASCARAS SE OBTUVIERON CON VENTANAS DE 20X20 WAVELENGTHS

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

%% FOR LOOPING
% Attenuation values from T1 to T8 and background, in that order
groundTruthTargets = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];
c1x = 1.9; c1z = 1.93; roiL = 0.9; roiD = 0.6;
attRange = [0.4,1.1];
bsRange = [-2 2];
%%
for iAcq = 1:8
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);

%% ----------------------- CURRENT APPROACH -----------------------

% Spectrum
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

%% SOLVING SYSTEM
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];


% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = mu/31;
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
    title(['RSLD, \mu=',num2str(mu(ii),2)])
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


%% ---------------------- British Columbia Approach ----------------------
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
title(tl,'Weights proposed by BC');
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
a = 5/4; b = 0.09;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')

%% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
BRW = zeros(m,n,length(mu));
CRW = zeros(m,n,length(mu));
mu2 = mu/31;
for mm = 1:length(mu)
    tic
    [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),mu(mm),mu2(mm),...
        m,n,tol,mask(:),w);
    toc
    BRW(:,:,mm) = (reshape(Bn*8.686,m,n));
    CRW(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 12]);
tl = tiledlayout(2,size(BRW,3)+1);
title(tl,'British Columbia Approach')
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')
for ii = 1:size(BRW,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BRW(:,:,ii), attRange)
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

for ii = 1:size(BRW,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,CRW(:,:,ii), bsRange)
    colormap(t2,parula)
    axis equal tight
    title(['RSLD, \mu=',num2str(mu2(ii),2)])
end
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';


%% ---------------------- MANUAL COMPENSATION ----------------------
% iAcq = 1;
% load([croppedDir,'\T',num2str(iAcq),'.mat'])
% h = fspecial("average",[50 5]);
% blurred = imfilter(Bmode,h,"symmetric");
% [BW,~] = segmentImage(blurred);
% save(['./masks/T',num2str(iAcq),'.mat'],'BW');

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
factor = ones(size(mask));
if sum(~mask(:)) ~=0
    factor(mask) = std(sam1(~mask))/std(sam1(mask));
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
mu2 = mu/10;
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

%% Getting metrics
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR(:,:,3),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
RSLD(iAcq) = r;

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRE(:,:,2),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
ERSLD(iAcq) = r;

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRW(:,:,3),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
WRSLD(iAcq) = r;

%% ROI
fig = figure(13);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
[~,~,hC] = imoverlay2(Bmode,AttInterp,dynRange,attRange,0.5,x,z,back|inc,x,z);
xlabel('x [cm]')
ylabel('y [cm]')
hC.Label.String = 'Attenuation [db/cm/MHz]';

%%
saveDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(saveDir,"dir")), mkdir(saveDir); end
save_all_figures_to_directory(saveDir);
close all

end


%% Results
% CAMBIAR NUMERO A LAS FIGURAS
resultsRSLD = struct2table(RSLD);
results1 = struct2table(ERSLD);
results2 = struct2table(WRSLD);

fig = figure(15);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
errorbar(resultsRSLD.meanInc, resultsRSLD.stdInc,'o', 'LineWidth',2)
hold on,
errorbar(results1.meanInc, results1.stdInc, 'x', 'LineWidth',2)
errorbar(results2.meanInc, results2.stdInc, 's', 'LineWidth',2)
plot(groundTruthTargets(1:8), "_k"	, 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 1.5])
grid on
legend({'RSLD','ERSLD','WRSLD'})
title('Results inclusion')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

fig = figure(16);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
errorbar(resultsRSLD.meanBack, resultsRSLD.stdBack,'o', 'LineWidth',2)
hold on,
errorbar(results1.meanBack, results1.stdBack, 'x', 'LineWidth',2)
errorbar(results2.meanBack, results2.stdBack, 's', 'LineWidth',2)
yline(groundTruthTargets(9), 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 1.5])
grid on
legend({'RSLD','ERSLD','WRSLD'})
title('Results background')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

%% CV
fig = figure(17);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
plot(resultsRSLD.stdInc./resultsRSLD.meanInc*100,'o', 'LineWidth',2)
hold on,
plot(results1.stdInc./results1.meanInc*100, 'x', 'LineWidth',2)
plot(results2.stdInc./results2.meanInc*100, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 30])
grid on
legend({'RSLD','ERSLD','WRSLD'})
title('Results inclusion')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

fig = figure(18);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
plot(resultsRSLD.stdBack./resultsRSLD.meanBack*100,'o', 'LineWidth',2)
hold on,
plot(results1.stdBack./results1.meanBack*100, 'x', 'LineWidth',2)
plot(results2.stdBack./results2.meanBack*100, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 30])
grid on
legend({'RSLD','ERSLD','WRSLD'})
title('Results background')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

%% Mean Percentage error
fig = figure(19);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
plot(resultsRSLD.MPEInc, 'o', 'LineWidth',2)
hold on,
plot(results1.MPEInc, 'x', 'LineWidth',2)
plot(results2.MPEInc, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([-30 100])
grid on
legend({'RSLD','ERSLD','WRSLD'}, 'Location','northwest')
title('Results inclusion')
xlabel('Target')
ylabel('MPE %')

fig = figure(20);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
plot(resultsRSLD.MPEBack, 'o', 'LineWidth',2)
hold on,
plot(results1.MPEBack, 'x', 'LineWidth',2)
plot(results2.MPEBack, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([-30 100])
grid on
legend({'RSLD','ERSLD','WRSLD'}, 'Location','northwest')
title('Results background')
xlabel('Target')
ylabel('MPE %')

%%
targetDir = fullfile(figDir,'Overall');
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);



