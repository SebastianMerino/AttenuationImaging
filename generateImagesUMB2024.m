% ====================================================================== %
% Script to generate images for the Journal Proposal of UMB. 
% Created on Dec, 2023
% ====================================================================== %

clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

figDir = 'C:\Users\sebas\Pictures\UMB2024\21-12';
if (~exist(figDir,"dir")), mkdir(figDir); end

% ========================================================================
% ========================================================================
% ========================================================================
%% SIMULATION ON LAYERED MEDIA

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered_14_11_23'];
croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

groundTruthTop = [0.6,0.6,0.6,1.2,1.2,1.2];
groundTruthBottom = [1.2,1.2,1.2,0.6,0.6,0.6];

%% Setting up
figure('Units','centimeters', 'Position',[5 5 25 8]);
tl = tiledlayout(2,5, "Padding","tight");

for iAcq = 1:2
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");

% Plotting
dynRange = [-40,0];
attRange = [0.4,1.4];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating reference
[~,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
attIdeal(Z<=2) = groundTruthTop(iAcq);
attIdeal(Z>2) = groundTruthBottom(iAcq);

%% TV
muB = 10^3; muC = 10^0.5;
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

axialTV{iAcq} = mean(BRTV,2);
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
top = Zq < 1.9;
bottom = Zq > 2.1;
AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;
%% SWTV
% Calculating SNR
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% Calculating weights
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% Method
muB = 10^2.5; muC = 10^0;
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB,muC,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

axialSWTV{iAcq} = mean(BRSWTV,2);
AttInterp = interp2(X,Z,BRSWTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% TVL1
muB = 10^3.5; muC = 10^1.25;
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

axialTVL1{iAcq} = mean(BRTVL1,2);
AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% WFR
% Computing weights
ratioCutOff = 6;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((CRTVL1/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% Setting up new system
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Method
muB = 10^4.5; muC = 10^2;
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB,muC,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);
CRWFR = reshape(Cn*NptodB,m,n);

axialWFR{iAcq} = mean(BRWFR,2);
AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;
%% Plotting

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1, 'westoutside')
title('Bmode')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('TV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t4,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

end

%% Axial profiles
figure('Units','centimeters', 'Position',[5 5 12 12])
tiledlayout(2,1)
nexttile,
plot(z_ACS, axialTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{1}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{1}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{1}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([z_ACS(1) z_ACS(end)])
title('Axial profiles')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

nexttile,
plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([z_ACS(1) z_ACS(end)])
title('Axial profiles')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 


%%
figure('Units','centimeters', 'Position',[5 5 15 5])
tiledlayout(1,2)

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';


%%
save_all_figures_to_directory(figDir,['sim',num2str(iAcq),'Figure']);
close all

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasTop)
disp(results2.biasTop)
disp(results3.biasTop)
disp(results4.biasTop)

disp('Bias Bottom')
disp(results1.biasBottom)
disp(results2.biasBottom)
disp(results3.biasBottom)
disp(results4.biasBottom)

disp('RMSE Top')
disp(results1.rmseTop)
disp(results2.rmseTop)
disp(results3.rmseTop)
disp(results4.rmseTop)

disp('RMSE Bottom')
disp(results1.rmseBottom)
disp(results2.rmseBottom)
disp(results3.rmseBottom)
disp(results4.rmseBottom)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)

%%
% ========================================================================
% ========================================================================
% ========================================================================
%% SIMULATION ON INCLUSIONS

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\Simulation_23_12_18'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

groundTruthTop = [0.6,0.6,0.6,1.2,1.2,1.2];
groundTruthBottom = [1.2,1.2,1.2,0.6,0.6,0.6];

%% Setting up
figure('Units','centimeters', 'Position',[5 5 25 12]);
tl = tiledlayout(3,5, "Padding","tight");

for iAcq = 1:3
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");

% Plotting
dynRange = [-40,0];
attRange = [0.4,1.4];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating reference
[X,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
rInc = 0.8;
inclusion = (X.^2 + (Z-2).^2)<= rInc^2;
attIdeal(~inclusion) = groundTruthTop(iAcq);
attIdeal(inclusion) = groundTruthBottom(iAcq); %incl = bottom

%% TV
switch iAcq
    case 1
        muB = 10^3.5; muC = 10^3;
    case 2
        muB = 10^3.5; muC = 10^2;
    case 3
        muB = 10^3.5; muC = 10^2;
end
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

axialTV{iAcq} = mean(BRTV(:,20:27),2);
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
inclusion = (Xq.^2 + (Zq-2).^2)<= rInc^2;
top = ~inclusion;
bottom = inclusion;
AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;
%% SWTV
% Calculating SNR
envelope = abs(hilbert(sam1));
SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

% Calculating weights
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% Method
switch iAcq
    case 1
        muB = 10^3; muC = 10^3;
    case 2
        muB = 10^3.5; muC = 10^2.5;
    case 3
        muB = 10^3.5; muC = 10^2.5;
end
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB,muC,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

axialSWTV{iAcq} = mean(BRSWTV(:,20:27),2);
AttInterp = interp2(X,Z,BRSWTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% TVL1
switch iAcq
    case 1
        muB = 10^3; muC = 10^1.5;
    case 2
        muB = 10^3; muC = 10^0.5;
    case 3
        muB = 10^3; muC = 10^0.5;
end
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

axialTVL1{iAcq} = mean(BRTVL1(:,20:27),2);
AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% WFR
% Computing weights
% ratioCutOff = 6;
% order = 5;
% reject = 0.1;
% extension = 3;
% w = (1-reject)*(1./((CRTVL1/ratioCutOff).^(2*order) + 1))+reject;
% w = movmin(w,extension);
w = db2mag(-sqrt(CRTVL1.^2 + 1));

% Setting up new system
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Method
switch iAcq
    case 1
        muB = 10^3.5; muC = 10^2;
    case 2
        muB = 10^3.5; muC = 10^1.5;
    case 3
        muB = 10^3.5; muC = 10^1.5;
end
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB,muC,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);
CRWFR = reshape(Cn*NptodB,m,n);

axialWFR{iAcq} = mean(BRWFR(:,20:27),2);
AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;
%% Plotting

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1, 'westoutside')
title('Bmode')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('TV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t4,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

end

%%
figure('Units','centimeters', 'Position',[5 5 15 5])
tiledlayout(1,2)
t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';

%% Axial profiles
figure('Units','centimeters', 'Position',[5 5 12 12])
tiledlayout(2,1)
nexttile,
plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal(:,20:27),2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([z_ACS(1) z_ACS(end)])
title('Axial profiles')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

nexttile,
plot(z_ACS, axialTV{3}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{3}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{3}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{3}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal(:,20:27),2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([z_ACS(1) z_ACS(end)])
title('Axial profiles')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

%%
save_all_figures_to_directory(figDir,['simInc',num2str(iAcq),'Figure']);
close all

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasTop)
disp(results2.biasTop)
disp(results3.biasTop)
disp(results4.biasTop)

disp('Bias Bottom')
disp(results1.biasBottom)
disp(results2.biasBottom)
disp(results3.biasBottom)
disp(results4.biasBottom)

disp('RMSE Top')
disp(results1.rmseTop)
disp(results2.rmseTop)
disp(results3.rmseTop)
disp(results4.rmseTop)

disp('RMSE Bottom')
disp(results1.rmseBottom)
disp(results2.rmseBottom)
disp(results3.rmseBottom)
disp(results4.rmseBottom)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)


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

%% For looping each phantom
figure('Units','centimeters', 'Position',[5 5 25 12]);
tl = tiledlayout(3,5, "Padding","tight");

for iAcq = 1:3
fprintf("Phantom no. %i, %s\n",iAcq,croppedFiles(iAcq+5).name);
load(fullfile(croppedDir,croppedFiles(iAcq+5).name));

switch iAcq
    case 1
        muTV = 10^3; mu2TV = 10^2.5;
        muWTV = 10^3; mu2WTV = 10^3;
        muTik = 10^3; mu2Tik = 10^1.5;
        muWTik = 10^3.5; mu2WTik = 10^2.5;
    case 2
        muTV = 10^3.5; mu2TV = 10^2;
        muWTV = 10^2.5; mu2WTV = 10^2;
        muTik = 10^3; mu2Tik = 10^0.5;
        muWTik = 10^4; mu2WTik = 10^1.5;
    case 3
        muTV = 10^3; mu2TV = 10^0.5;
        muWTV = 10^3; mu2WTV = 10^1.5;
        muWTik = 10^3.5; mu2WTik = 10^1;
        muTik = 10^3; mu2Tik = 10^0;
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
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
r.rmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsTV(iAcq) = r;

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
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
r.rmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsSWTV(iAcq) = r;

%% TVL1
tic
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muTik,mu2Tik,m,n,tol,mask(:));
toc
BRTik = (reshape(Bn*NptodB,m,n));
CRTik = (reshape(Cn,m,n));


AttInterp = interp2(X,Z,BRTik,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
r.rmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsTVL1(iAcq) = r;

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
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.rmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") ;
r.rmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsWFR(iAcq) = r;

%%
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';
fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('TV')
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

t3 = nexttile;
imagesc(x_ACS,z_ACS,BRBC, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title('SWTV')
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

t4 = nexttile;
imagesc(x_ACS,z_ACS,BRTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t4,turbo)
axis image
title('TVL1')
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('WFR')
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
end

%%
save_all_figures_to_directory(figDir,'phantom');
close all
%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Inc')
disp(results1.biasInc)
disp(results2.biasInc)
disp(results3.biasInc)
disp(results4.biasInc)

disp('Bias Back')
disp(results1.biasBack)
disp(results2.biasBack)
disp(results3.biasBack)
disp(results4.biasBack)

disp('RMSE Inc')
disp(results1.rmseInc)
disp(results2.rmseInc)
disp(results3.rmseInc)
disp(results4.rmseInc)

disp('RMSE Back')
disp(results1.rmseBack)
disp(results2.rmseBack)
disp(results3.rmseBack)
disp(results4.rmseBack)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)
%%
%save_all_figures_to_directory(figDir,'phantom');
%close all
