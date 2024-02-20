% Script that generates images for ISBI 2024
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% ======================================================================
% ======================================================================
% ======================================================================
%% SIMULATION
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\layeredNew'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered_14_11_23'];
croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);

figDir = 'C:\Users\sebas\Pictures\ISBI2024\v5';
mkdir(figDir)
%%
for iAcq = [1,2,3,6]
switch iAcq
    case 1                
        muTV = 10^3; mu2TV = 10^0.5;
        muWTV = 10^2.5; mu2WTV = 10^0;
        muWTik = 10^4.5; mu2WTik = 10^2;
    case 2
        muTV = 10^4; mu2TV = 10^2;
        muWTV = 10^3.5; mu2WTV = 10^3;
        muWTik = 10^4; mu2WTik = 10^2;
    case 3
        muTV = 10^4; mu2TV = 10^5;
        muWTV = 10^3; mu2WTV = 10^0;
        muWTik = 10^4; mu2WTik = 10^2;
    case 6
        muTV = 10^3.5; mu2TV = 10^1;
        muWTV = 10^3; mu2WTV = 10^0;
        muWTik = 10^4.5; mu2WTik = 10^1.5;
end
muB = 10^3.5; muC = 10^1;

NptodB = log10(exp(1))*20;
fprintf("Simulation no. %i, %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");

% Setting up equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);

% GT
% groundTruthTop = [0.5,1,1,0.5,1,1];
% groundTruthBottom = [1,0.5,1,1,0.5,1];
groundTruthTop = [0.6,0.6,0.6,1.2,1.2,1.2];
groundTruthBottom = [1.2,1.2,1.2,0.6,0.6,0.6];

% Creating reference
[~,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
attIdeal(Z<=2) = groundTruthTop(iAcq);
attIdeal(Z>2) = groundTruthBottom(iAcq);
%% RSLD
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muTV,mu2TV,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

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

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% RSLD ANISOTROPIC AND BS WEIGHTED
tic
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muWTV,mu2WTV,m,n,tol,mask(:),w);
toc
BRBC = (reshape(Bn*NptodB,m,n));
CRBC = (reshape(Cn,m,n));

%% WEIGHTS FROM TVL1
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight function
% w = 1./((bscMap/10).^2 + 1);
ratioCutOff = 6;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);
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

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muWTik,mu2WTik,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*NptodB,m,n));
CRWTik = (reshape(Cn,m,n));

%% Metrics
fprintf('\nMetrics in simulation %i\n',iAcq);


[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR,Xq,Zq);

% top = Zq < 2.3;
% bottom = Zq > 2.5;
top = Zq < 1.9;
bottom = Zq > 2.1;

r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
    groundTruthTop(iAcq),"omitnan") * 100;
r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
    groundTruthBottom(end),"omitnan") * 100;
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;
fprintf('RSLD-TV\n');
fprintf('%.2f\n',r.MPETop);
fprintf('%.2f\n',r.MPEBottom);
fprintf('%.2f\n',r.cnr);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRBC,Xq,Zq);

r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
    groundTruthTop(iAcq),"omitnan") * 100;
r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
    groundTruthBottom(end),"omitnan") * 100;
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;
fprintf('RSLD-SWTV\n');
fprintf('%.2f\n',r.MPETop);
fprintf('%.2f\n',r.MPEBottom);
fprintf('%.2f\n',r.cnr);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRWTik,Xq,Zq);

r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
    groundTruthTop(iAcq),"omitnan") * 100;
r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
    groundTruthBottom(end),"omitnan") * 100;
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;
fprintf('RSLD-WFR\n');
fprintf('%.2f\n',r.MPETop);
fprintf('%.2f\n',r.MPEBottom);
fprintf('%.2f\n',r.cnr);


%% Plotting three results
rmseTV = sqrt(mean((BR-attIdeal).^2,'all'));
rmseSWTV = sqrt(mean((BRBC-attIdeal).^2,'all'));
rmseWFR = sqrt(mean((BRWTik-attIdeal).^2,'all'));

dynRange = [-50,0];
attRange = [0.4,1.4];
bsRange = [-2 2];

figure('Units','centimeters', 'Position',[5 5 5 4.3]);
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(gray)
title('B-mode')
subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 4.3]);
imagesc(x_ACS,z_ACS,attIdeal,attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('Ideal')
subtitle(' ')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 4.3]);
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-TV')
subtitle(['RMSE:',num2str(rmseTV,2),', CNR:',num2str(MetricsTV(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 4.3]);
imagesc(x_ACS,z_ACS,BRBC, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-SWTV')
subtitle(['RMSE:',num2str(rmseSWTV,2),', CNR:',num2str(MetricsSWTV(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 4.3]);
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-WFR')
subtitle(['RMSE:',num2str(rmseWFR,2),', CNR:',num2str(MetricsWFR(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

%% Profile

% attTV = mean(BR,2);
% attSWTV = mean(BRBC,2);
% attSWTVTik = mean(BRWTik,2);
% 
% bscTV = mean(CR,2);
% bscSWTV = mean(CRBC,2);
% bscSWTVTik = mean(CRWTik,2);
% 
% attGT = zeros(size(z_ACS));
% attGT(z_ACS < 2) = groundTruthTop(iAcq);
% attGT(z_ACS > 2) = groundTruthBottom(iAcq);
% 
% figure('Units','centimeters', 'Position',[5 5 8 4]);
% 
% plot(z_ACS,[attTV,attSWTV,attSWTVTik], 'LineWidth',1.5)
% hold on
% plot(z_ACS,attGT, 'k--')
% hold off
% title('Axial profile')
% grid on
% xlim([z_ACS(1),z_ACS(end)])
% ylim([0.5 1.3])
% xlabel('Depth [cm]'), ylabel('ACS [dB/cm/MHz]')
% legend({'TV','SWTV','WFR'}, 'Location','southeast')
% fontsize(gcf,8,'points')

%% Save
save_all_figures_to_directory(figDir,['simulation',num2str(iAcq),'fig']);
close all

end

%%


% ======================================================================
% ======================================================================
% ======================================================================
%% PHANTOMSSS
clear,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
     '\ID316V2\06-08-2023-Generic'];
croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
NptodB = log10(exp(1))*20;

figDir = 'C:\Users\sebas\Pictures\ISBI2024\v5';

%% For looping each phantom

for iAcq = 1:3
fprintf("Phantom no. %i, %s\n",iAcq,croppedFiles(iAcq+5).name);
load(fullfile(croppedDir,croppedFiles(iAcq+5).name));

switch iAcq
    case 1
        muTV = 10^3; mu2TV = 10^2.5;
        muWTV = 10^3; mu2WTV = 10^3;
        muWTik = 10^3.5; mu2WTik = 10^2.5;
    case 2
        % muTV = 10^3.5; mu2TV = 10^2;
        % muWTV = 10^3; mu2WTV = 10^2.5;
        % muWTik = 10^4; mu2WTik = 10^2;
        muTV = 10^3.5; mu2TV = 10^2;
        muWTV = 10^2.5; mu2WTV = 10^2;
        muWTik = 10^4; mu2WTik = 10^1.5;
    case 3
        %muTV = 10^3; mu2TV = 10^0.5;
        muTV = 10^3; mu2TV = 10^0.5;
        %muWTV = 10^2.5; mu2WTV = 10^1;
        muWTV = 10^3; mu2WTV = 10^1.5;
        muWTik = 10^3.5; mu2WTik = 10^1;
end
muB = 10^3; muC = 10^0.5;

groundTruthTargets = [0.97,0.95,0.95,0.55];

%% Plotting B-mode
dynRange = [-50,0];
attRange = [0.4,1.1];
bsRange = [-2 2];


%% ROI selection
c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.2;
% roiL = 0.8; roiD = 0.6;

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);

rI = 0.6; rB = 1.2; % Both
% inc = (Xq - c1x).^2 + (Zq - c1z).^2 < rI^2;
% back = (Xq - c1x).^2 + (Zq - c1z).^2 > rB^2;

x0mask = c1x - roiL/2; 
z0mask = c1z - roiL/2;
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);
%figure, imagesc(x,z,inc|back)

%% RSLD-TV

b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muTV,mu2TV,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BR,Xq,Zq);
RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
rmseTV = sqrt((RmseInc + RmseBack)/2);

%% British Columbia Approach

% Weights
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
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );

% Regularization: Au = b
tol = 1e-3;
tic
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muWTV,mu2WTV,m,n,tol,mask(:),w);
%[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muTV,mu2TV,m,n,tol,mask(:),w);
toc
BRBC = (reshape(Bn*NptodB,m,n));
CRBC = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRBC,Xq,Zq);
RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
rmseSWTV = sqrt((RmseInc + RmseBack)/2);
%% NEW WEIGHTS
% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

ratioCutOff = 6;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% Weighting equation and regularizations
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

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muWTik,mu2WTik,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*NptodB,m,n));
CRWTik = (reshape(Cn,m,n));


AttInterp = interp2(X,Z,BRWTik,Xq,Zq);
RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
rmseWFR = sqrt((RmseInc + RmseBack)/2);
%%
AttInterp = interp2(X,Z,BR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsTV(iAcq) = r;


AttInterp = interp2(X,Z,BRBC,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsSWTV(iAcq) = r;


AttInterp = interp2(X,Z,BRWTik,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsWFR(iAcq) = r;
%%

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 4.2]);
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-TV')
subtitle(['RMSE:',num2str(rmseTV,2),', CNR:',num2str(MetricsTV(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

figure('Units','centimeters', 'Position',[5 5 5 4.2]);
imagesc(x_ACS,z_ACS,BRBC, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-SWTV')
subtitle(['RMSE:',num2str(rmseSWTV,2),', CNR:',num2str(MetricsSWTV(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off


figure('Units','centimeters', 'Position',[5 5 5 4.2]);
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-WFR')
subtitle(['RMSE:',num2str(rmseWFR,2),', CNR:',num2str(MetricsWFR(iAcq).cnr,2)])
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off
%%
save_all_figures_to_directory(figDir,['phantom',num2str(iAcq),'fig']);
close all
end

%%
for iAcq = 1:3
    fprintf("T%i\n",iAcq);
    r = MetricsTV(iAcq);
    fprintf("TV,   MPEinc: %.2f, MPEback: %.2f, CNR: %.2f, \n",...
        r.MPEInc, r.MPEBack, r.cnr)
    r = MetricsSWTV(iAcq);
    fprintf("SWTV, MPEinc: %.2f, MPEback: %.2f, CNR: %.2f, \n",...
        r.MPEInc, r.MPEBack, r.cnr)
    r = MetricsWFR(iAcq);
    fprintf("WFR,  MPEinc: %.2f, MPEback: %.2f, CNR: %.2f, \n\n",...
        r.MPEInc, r.MPEBack, r.cnr)
end

%%
%save_all_figures_to_directory(figDir,'phantom');
%close all

% ========================================================================
% ========================================================================
% ========================================================================
%% Clinical case
clear,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\CUELLO#3'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];
figDir = 'C:\Users\sebas\Pictures\ISBI2024\v5';

targetFiles = dir([targetDir,'\*.mat']);

blocksize = 15;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 1;
NptodB = log10(exp(1))*20;

% attRange = [0,1.7];
attRange = [0.3,1.8];
%% Loading case
iAcq = 8;
for iRoi = 1:2
load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = RF(:,:,1);
% dynRange = [-50,0];


BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));
% figure('Units','centimeters', 'Position',[5 5 15 15]),
% imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
% ylim([0.05 3])

if iRoi == 1
    rect = [1.03; 0.49; 1.6; 1.69]; % Previous rectangle
    % rect = [1.03; 0.4; 1.6; 1.8];
    muTV = 10^3; mu2TV = 10^0.5;
    muWTV = 10^3; mu2WTV = 10^0.5;
    muWTik = 10^3; mu2WTik = 10^0.5;
else
    rect = [2.63; 0.49; 1.6; 1.69]; % Previous rectangle
    % rect = [2.63; 0.4; 1.6; 1.8];
    muTV = 10^4; mu2TV = 10^2;
    muWTV = 10^4; mu2WTV = 10^2;
    muWTik = 10^4; mu2WTik = 10^2;
end
% hold on
% rectangle('Position',rect)
% hold off
%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);

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
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = attenuation_phantoms_Np(f, 4, []); % CAMBIAESTO
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
refFiles = dir([refDir,'\*.mat']);
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load([refDir,'\',refFiles(iRef).name]);
    samRef = out.RF;
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

%% RSLD-TV

b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muTV,mu2TV,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

%% British Columbia Approach

% Weights
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
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );

% Regularization: Au = b
tol = 1e-3;
tic
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muWTV,mu2WTV,m,n,tol,mask(:),w);
toc
BRBC = (reshape(Bn*NptodB,m,n));
CRBC = (reshape(Cn,m,n));

%% NEW WEIGHTS
tol = 1e-3;
mask = ones(m,n,p);
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muWTik,mu2WTik,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
% test = reshape(Bn,m,n)*NptodB;
% figure,
% imagesc(x_ACS,z_ACS,test);
% axis image
% colormap turbo
% colorbar

% w = 1./((bscMap/10).^2 + 1);

% Weight function
ratioCutOff = 6;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% figure, imagesc(w)

% Weighting equation and regularizations
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

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muWTik,mu2WTik,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*NptodB,m,n));
CRWTik = (reshape(Cn,m,n));

imageData.x = x_ACS;
imageData.z = z_ACS;
imageData.roi = roi;
imageData.TV = BR;
imageData.SWTV = BRBC;
imageData.WFR = BRWTik;
dataRoi{iRoi} = imageData;

end

%% PLOTTING IMAGES

alpha = 0.7;
dynRange = [-40 -10];

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(xFull,zFull,BmodeFull); axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('B-mode')
fontsize(gcf,8,'points')


figure('Units','centimeters', 'Position',[5 5 5 3.8]);
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.TV,dynRange,attRange,alpha,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.TV,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
iRoi = 2;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('RSLD-TV')
hColor.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')


figure('Units','centimeters', 'Position',[5 5 5 3.8]);
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.SWTV,dynRange,attRange,alpha,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.SWTV,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
iRoi = 2;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('RSLD-SWTV')
hColor.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
iRoi = 1;
[~,hB,hColor] = imOverlayInterp(BmodeFull,dataRoi{iRoi}.WFR,dynRange,attRange,alpha,...
    dataRoi{iRoi}.x,dataRoi{iRoi}.z,dataRoi{iRoi}.roi,xFull,zFull);
% Interpolation
iRoi = 2;
[X,Z] = meshgrid(dataRoi{iRoi}.x,dataRoi{iRoi}.z);
[Xq,Zq] = meshgrid(xFull,zFull);
imgInterp = interp2(X,Z,dataRoi{iRoi}.WFR,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion & dataRoi{iRoi}.roi;
% Overlap
hold on;
iRoi = 2;
hF = imagesc(dataRoi{iRoi}.x,dataRoi{iRoi}.z,imgInterp,attRange);
set(hF,'XData',get(hB,'XData'),'YData',get(hB,'YData'))
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
ylim([0.05 3])
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
title('RRSLD-WFR')
hColor.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')

%%
fprintf("Homogeneous results: \n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{2}.TV(:)),...
    std(dataRoi{2}.TV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{2}.SWTV(:)),...
    std(dataRoi{2}.SWTV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataRoi{2}.WFR(:)),...
    std(dataRoi{2}.WFR(:)))

[X,Z] = meshgrid(dataRoi{1}.x,dataRoi{1}.z);
maskThyroid = Z>1.3;
dataTV = dataRoi{1}.TV(maskThyroid);
dataSWTV = dataRoi{1}.SWTV(maskThyroid);
dataWFR = dataRoi{1}.WFR(maskThyroid);
fprintf("\nHeterogeneous results: \n BOTTOM\n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataTV(:)),std(dataTV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataSWTV(:)),std(dataSWTV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataWFR(:)),std(dataWFR(:)))

maskNodule = Z<1.1;
dataTV = dataRoi{1}.TV(maskNodule);
dataSWTV = dataRoi{1}.SWTV(maskNodule);
dataWFR = dataRoi{1}.WFR(maskNodule);
fprintf("\nHeterogeneous results: \n TOP\n")
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataTV(:)),std(dataTV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataSWTV(:)),std(dataSWTV(:)))
fprintf("Mean: %.2f, Std: %.2f\n",mean(dataWFR(:)),std(dataWFR(:)))

%%

save_all_figures_to_directory(figDir,'clinical');
close all