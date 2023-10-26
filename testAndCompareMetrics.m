clear,clc
close all
addpath('./functions_v7');
%%
targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig\26-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% FOR LOOPING
% Attenuation values from T1 to T8 and background, in that order
%groundTruthTargets = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];
groundTruthTargets = [0.97,0.95,0.95,0.55];

muTV = [3.2e3,1e3,1e3];
mu2TV = [1e3,10,3.2];
muTik = [3.2e3,1e3,1e3];
mu2Tik = [1e2,3.2,1];
muWTik = [3.2e3,3.2e3,3.2e3];
mu2WTik = [1e2,10,3.2];

c1x = 1.9; c1z = 1.93;
%roiL = 0.9; roiD = 0.6;
roiL = 1;roiD = 0.6;

%%
for iAcq = 1:3
load([croppedDir,'\T',num2str(iAcq+5),'.mat'])
iAcq = iAcq - 5;
load([refDir,'\compensation.mat']);

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

        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

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


%% Regularized SLD
% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);

mu = muTV(iAcq);
mu2 = mu2TV(iAcq);

BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR(:,:,1),Xq,Zq);

rI = 0.6; rB = 1.2; % Both
inc = (Xq - c1x).^2 + (Zq - c1z).^2 < rI^2;
back = (Xq - c1x).^2 + (Zq - c1z).^2 > rB^2;

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
RSLD(iAcq) = r;


% Plotting
% figure('Units','centimeters', 'Position',[5 5 30 8]);
% tiledlayout(1,size(BR,3)+1);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% for ii = 1:size(BR,3)
%     t2 = nexttile; 
%     imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
%     colormap(t2,turbo)
%     axis equal tight
%     title(['RSLD, \mu=',num2str(mu(ii),2)])
% end
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';

%saveas(gcf,[figDir,'\T',num2str(iAcq),'.png'])


%% Au = b
% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = muTik(iAcq);
mu2 = mu2Tik(iAcq);
tic
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));
toc
BRT(:,:,1) = (reshape(Bn*8.686,m,n));
CR(:,:,1) = (reshape(Cn,m,n));
% mu = logspace(2.5,3.5,3);
% mu2 = logspace(-0.5,0.5,3);
% BRT = zeros(m,n,length(mu2));
% CR = zeros(m,n,length(mu2));
% for mm = 1:length(mu)
%     for mm2 = 1:length(mu2)
%         tic
%         [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
%         toc
%         BRT(:,:,mm2) = (reshape(Bn*8.686,m,n));
%         CR(:,:,mm2) = (reshape(Cn,m,n));
%     end
% 
%     % Plotting
%     figure('Units','centimeters', 'Position',[5 5 30 12]);
%     tl = tiledlayout(2,size(BRT,3)+1);
%     title(tl,'RSLD with isotropic TV and Tikhonov reg.')
%     t1 = nexttile;
%     imagesc(x,z,Bmode,dynRange)
%     axis equal
%     xlim([x_ACS(1) x_ACS(end)]),
%     ylim([z_ACS(1) z_ACS(end)]),
%     colormap(t1,gray)
%     colorbar(t1,'westoutside')
%     title('Bmode')
% 
%     for ii = 1:size(BRT,3)
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,BRT(:,:,ii), attRange)
%         colormap(t2,turbo)
%         axis equal tight
%         title(['RSLD, \mu=',num2str(mu(mm),2)])
%     end
%     c = colorbar;
%     c.Label.String = 'Att. [db/cm/MHz]';
% 
%     t3 = nexttile;
%     imagesc(x,z,Bmode,dynRange)
%     axis equal
%     xlim([x_ACS(1) x_ACS(end)]),
%     ylim([z_ACS(1) z_ACS(end)]),
%     colormap(t3,gray)
%     colorbar(t3,'westoutside')
%     title('Bmode')
% 
%     for ii = 1:size(BRT,3)
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
%         colormap(t2,parula)
%         axis equal tight
%         title(['RSLD, \mu=',num2str(mu2(ii),2)])
%     end
%     c = colorbar;
%     c.Label.String = 'BS log ratio (a.u.)';
% end

%%

% Getting metrics
%[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

AttInterp = interp2(X,Z,BRT(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
RSLDT(iAcq) = r;


%% NEW WEIGHTS
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );


% Regularization: Au = b
mu = 1e3; mu2 = 1;
tol = 1e-3;
mask = ones(m,n,p);
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));
bscMap = (reshape(Cn,m,n));

logBscRatio = bscMap*log10(exp(1))*20;
w = 1./((logBscRatio/6).^2 + 1);


% Weights
figure('Units','centimeters', 'Position',[5 5 30 8]),
tiledlayout(1,3)
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
title('Log BSC Ratio')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
colorbar;
axis image
title('Weights')


%% Au = b
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;
Aw = [A1w A2w];

% Regularization: Au = b
tol = 1e-3;

mask = ones(m,n,p);
mu = muWTik(iAcq);
mu2 = mu2WTik(iAcq);
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu,mu2,m,n,tol,mask(:),w);
BRW2(:,:,1) = (reshape(Bn*8.686,m,n));
CR(:,:,1) = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRW2(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
WRSLD(iAcq) = r;

%% Plotting three results
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BR(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['RSLD, \mu=',num2str(muTV(iAcq),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRT(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['Tikhonov and TV, \mu=',num2str(muTik(iAcq),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRW2(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['Weighted SLD, \mu=',num2str(muWTik(iAcq),2)])

c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%% ROI
figure('Units','centimeters', 'Position',[5 5 12 8])
[~,~,hC] = imoverlay2(Bmode,AttInterp,dynRange,attRange,0.5,x,z,back|inc,x,z);
xlabel('x [cm]')
ylabel('y [cm]')
hC.Label.String = 'Attenuation [db/cm/MHz]';


end

%% Results
resultsRSLD = struct2table(RSLD);
resultsWRSLD1 = struct2table(RSLDT);
resultsWRSLD2 = struct2table(WRSLD);

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(resultsRSLD.meanInc, resultsRSLD.stdInc,'o', 'LineWidth',2)
hold on,
errorbar(resultsWRSLD1.meanInc, resultsWRSLD1.stdInc, 'x', 'LineWidth',2)
errorbar(resultsWRSLD2.meanInc, resultsWRSLD2.stdInc, 's', 'LineWidth',2)
plot(groundTruthTargets(1:end-1), "_k"	, 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0.3 1.3])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'})
title('Results inclusion')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(resultsRSLD.meanBack, resultsRSLD.stdBack,'o', 'LineWidth',2)
hold on,
errorbar(resultsWRSLD1.meanBack, resultsWRSLD1.stdBack, 'x', 'LineWidth',2)
errorbar(resultsWRSLD2.meanBack, resultsWRSLD2.stdBack, 's', 'LineWidth',2)
yline(groundTruthTargets(end), 'k--', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0.3 1.3])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'})
title('Results background')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

%% CV
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.stdInc./resultsRSLD.meanInc*100,'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.stdInc./resultsWRSLD1.meanInc*100, 'x', 'LineWidth',2)
plot(resultsWRSLD2.stdInc./resultsWRSLD2.meanInc*100, 's', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0 40])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'})
title('Results inclusion')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.stdBack./resultsRSLD.meanBack*100,'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.stdBack./resultsWRSLD1.meanBack*100, 'x', 'LineWidth',2)
plot(resultsWRSLD2.stdBack./resultsWRSLD2.meanBack*100, 's', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0 40])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'})
title('Results background')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

%% Mean Percentage error
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.MPEInc, 'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.MPEInc, 'x', 'LineWidth',2)
plot(resultsWRSLD2.MPEInc, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 3.5])
%ylim([-30 80])
ylim([-30 30])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'}, 'Location','northeast')
title('Results inclusion')
xlabel('Target')
ylabel('MPE %')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.MPEBack, 'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.MPEBack, 'x', 'LineWidth',2)
plot(resultsWRSLD2.MPEBack, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 3.5])
%ylim([-30 80])
ylim([-30 30])
grid on
legend({'RSLD','RSLDv2','WRSLDv2'}, 'Location','northeast')
title('Results background')
xlabel('Target')
ylabel('MPE %')
%%
save_all_figures_to_directory(figDir);


