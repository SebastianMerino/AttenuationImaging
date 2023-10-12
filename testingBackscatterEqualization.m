clear,clc
close all
addpath('./functions_v7');

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID316V2\06-08-2023-Generic'];

% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% FOR LOOPING
% Attenuation values from T1 to T8 and background, in that order
groundTruthTargets = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];
c1x = 1.9; c1z = 1.93; roiL = 0.9; roiD = 0.6;
%%
for iAcq = 1:8
% iAcq = 8;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);

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
% 
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

% save('T4.mat')

%% ----------------------- SOLVING SYSTEM -----------------------
% load('T4.mat')

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];


% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = 1e4*[0.4,3.6];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

%% ---------------------- MANUAL COMPENSATION ----------------------

% Inclusion
% cx = 1.92; cz = 1.98; r = 0.98; % T1
% cx = 1.87; cz = 1.92; r = 0.98; % T2
% cx = 1.90; cz = 1.87; r = 0.96; % T4
% cx = 1.95; cz = 1.93; r = 0.95; % T7
% cx = 1.85; cz = 1.93; r = 0.93; % T8
% figure('Units','centimeters', 'Position',[5 5 30 8]), 
% tiledlayout(1,3)
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% colormap(t1,gray)
% colorbar
% axis image
% title('B-mode')
% hold on
% rectangle('Position',[cx-r,cz-r,r+r,r+r], 'Curvature',[1 1], ...
%     'LineStyle','--', 'LineWidth',2)
% hold off
% 
% % Equalization
% [X,Z] = meshgrid(x,z);
% mask = ( (X-cx).^2 + (Z-cz).^2 ) < r*r;
% samEnhanced = sam1;
% factor = ones(size(mask));
% factor(mask) = std(samEnhanced(~mask))/std(samEnhanced(mask));
% h = fspecial("average",[50 5]);
% factor = imfilter(factor,h,"symmetric");
% t2 = nexttile;
% imagesc(x,z,factor)
% colormap(t2,parula)
% colorbar
% axis image
% title('Mask')
% 
% samEnhanced = sam1.*factor;
% Bmode2 = db(hilbert(samEnhanced));
% Bmode2 = Bmode2 - max(Bmode2(:));
% t3 = nexttile;
% imagesc(x,z,Bmode2,dynRange)
% colormap(t3,gray)
% colorbar
% axis image
% title('Equalized B-mode')
% 
% %BW = ones(size(Bmode));
% h = fspecial("average",[50 5]);
% blurred = imfilter(Bmode,h,"symmetric");

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
load(['maskT',num2str(iAcq),'.mat']);
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
Bmode2 = db(hilbert(samEnhanced));
Bmode2 = Bmode2 - max(Bmode2(:));
t3 = nexttile;
imagesc(x,z,Bmode2,dynRange)
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
BRE = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BRE(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

%% Plotting
% figure('Units','centimeters', 'Position',[5 5 20 15]);
% tiledlayout(2,2);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
% colormap(t1,gray)
% title('Bmode')
% 
% t3 = nexttile;
% imagesc(x,z,Bmode2,dynRange)
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
% colormap(t3,gray)
% colorbar(t3)
% title('Equalized Bmode')
% 
% 
% t2 = nexttile; 
% imagesc(x_ACS,z_ACS,BR(:,:,1), attRange)
% colormap(t2,turbo)
% axis equal tight
% title(['RSLD, \mu=',num2str(mu(1),2)])
% 
% t2 = nexttile; 
% imagesc(x_ACS,z_ACS,BRE(:,:,1), attRange)
% colormap(t2,turbo)
% axis equal tight
% title(['RSLD, \mu=',num2str(mu(1),2)])
% 
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';

%% Getting local weights
%----------
gamma = 6;
%-------------

envelope = abs(hilbert(samEnhanced));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        
        % SNR(ii,jj) = min([mean(sub_block_p)/std(sub_block_p),...
        %     mean(sub_block_d)/std(sub_block_d)]);

        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% Weights
% fig = figure(11);
% fig.Units = 'centimeters'; fig.Position = [5 5 30 8]; 
% tiledlayout(1,3)
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% colormap(t1,gray)
% colorbar
% axis equal
% xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% % axis image
% title('B-mode')
% 
% t2 = nexttile;
% imagesc(x_ACS,z_ACS,db(SNR))
% colormap(t2,parula)
% c = colorbar;
% ylabel(c,'dB')
% axis image
% title('SNR')
% 
% 
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = (SNRopt./SNR);
w = desvSNR.^gamma.*exp(1-desvSNR.^gamma);
% 
% t3 = nexttile;
% imagesc(x_ACS,z_ACS,w, [0 1])
% colormap(t3,parula)
% colorbar;
% axis image
% title(['Weights, order=',num2str(gamma)])

%% Au = b
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;
Aw = [A1w A2w];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
BREW = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1w,A2w,bw,mu1,mu2,m,n,tol,mask(:));

    BREW(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

fig = figure(12);
fig.Units = 'centimeters'; fig.Position = [5 5 30 15]; 

tiledlayout(2,3);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
title('Bmode')

t3 = nexttile;
imagesc(x,z,Bmode2,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t3,gray)
colorbar(t3)
title('Equalized Bmode')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
colorbar;
axis image
title(['Weights, order=',num2str(gamma)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BR(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['RSLD, \mu=',num2str(mu(1),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRE(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['RSLD, \mu=',num2str(mu(1),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BREW(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['RSLD, \mu=',num2str(mu(1),2)])

c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%% Getting metrics
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR(:,:,1),Xq,Zq);

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
AttInterp = interp2(X,Z,BRE(:,:,1),Xq,Zq);

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
AttInterp = interp2(X,Z,BREW(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
WERSLD(iAcq) = r;

%% ROI
fig = figure(13);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
[~,~,hC] = imoverlay2(Bmode,AttInterp,dynRange,attRange,0.5,x,z,back|inc,x,z);
xlabel('x [cm]')
ylabel('y [cm]')
hC.Label.String = 'Attenuation [db/cm/MHz]';

%%
targetDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);
close all

end


%% Results
% CAMBIAR NUMERO A LAS FIGURAS
resultsRSLD = struct2table(RSLD);
results1 = struct2table(ERSLD);
results2 = struct2table(WERSLD);

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
legend({'RSLD','ERSLD','WERSLD'})
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
legend({'RSLD','ERSLD','WERSLD'})
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
plot(groundTruthTargets(1:8), "_k"	, 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 60])
grid on
legend({'RSLD','ERSLD','WERSLD'})
title('Results inclusion')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

fig = figure(18);
fig.Units = 'centimeters'; fig.Position = [5 5 12 8]; 
plot(resultsRSLD.stdBack./resultsRSLD.meanBack*100,'o', 'LineWidth',2)
hold on,
plot(results1.stdBack./results1.meanBack*100, 'x', 'LineWidth',2)
plot(results2.stdBack./results2.meanBack*100, 's', 'LineWidth',2)
yline(groundTruthTargets(9), 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 60])
grid on
legend({'RSLD','ERSLD','WERSLD'})
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
legend({'RSLD','ERSLD','WERSLD'}, 'Location','northwest')
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
legend({'RSLD','ERSLD','WERSLD'}, 'Location','northwest')
title('Results background')
xlabel('Target')
ylabel('MPE %')

%%
targetDir = fullfile(figDir,'Overall');
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);



