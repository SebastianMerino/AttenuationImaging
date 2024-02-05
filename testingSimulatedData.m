% ====================================================================== %
% Script for testing simulated data with an inclusion 
% Created on Dec, 2023
% ====================================================================== %

clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\simulations_processed\24_01_26'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\Simulation_23_12_18'];
croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

figDir = [baseDir,'\fig\24-01-31'];
if (~exist(figDir,"dir")), mkdir(figDir); end

% Thyroid
% groundTruthBack = [1.5,1.3,1.5,1.5,1.5];
% groundTruthInc = [0.8,0.8,1,0.8,0.8];

% normal inclusion
groundTruthBack = [0.6,0.6,0.6];
groundTruthInc = [1.2,1.2,1.2];

% Plotting
dynRange = [-40,0];
attRange = [0.4,1.4]; % 0.5-1.8
bsRange = [-15 15];
NptodB = log10(exp(1))*20;


rInc = 0.7; % 0.8 for other inclusions

%     cz = 20e-3; cx = 0;
%     r = 7e-3;
%% Loading data
for iAcq = 3:3 % length(croppedFiles)
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");


% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating reference
% [~,Z] = meshgrid(x_ACS,z_ACS);
% attIdeal = ones(size(Z));
% attIdeal(Z<=2) = groundTruthTop(iAcq);
% attIdeal(Z>2) = groundTruthBottom(iAcq);
[X,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
inc = (X.^2 + (Z-2).^2)<= rInc^2;
attIdeal(~inc) = groundTruthBack(iAcq);
attIdeal(inc) = groundTruthInc(iAcq); %incl = bottom

figure('Units','centimeters', 'Position',[5 5 20 6]);
tiledlayout(1,2)
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)])
ylim([z_ACS(1) z_ACS(end)])
colormap(t1,gray)
colorbar('westoutside')
title('Bmode')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
axis image
colormap(t2,turbo);
c = colorbar('westoutside');
c.Label.String = 'dB/cm/MHz';
title('Ideal ACS')

%%
c1x = 0; c1z = 1.95;
roiL = 0.7; roiD = 0.5;
roiLz = 1.2;
%roiL = 1.2; roiD = 0.6;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);
figure, 
imagesc(x,z,Bmode, [-50 0])
axis equal
xlim([x_ACS(1) x_ACS(end)])
ylim([z_ACS(1) z_ACS(end)])
colormap(gray)
colorbar('westoutside')
title('Bmode')
hold on
contour(x,z,back, 1, 'w--')
contour(x,z,inc, 1, 'w--')
hold off

%% Spectrum
% Heterogeneous
[X,Z] = meshgrid(x_ACS,z_ACS);
region1 = (X.^2 + (Z-2).^2) <= rInc^2;
region2 = (X.^2 + (Z-2).^2) >= rInc^2;
sld1 = squeeze(sum(sum(b.*region1,1),2))/sum(region1(:)) * NptodB /4/L;
acs1 = f\sld1;
fprintf('Attenuation is %.2f\n',acs1)
sld2 = squeeze(sum(sum(b.*region2,1),2))/sum(region2(:)) * NptodB /4/L;
acs2 = f\sld2;
fprintf('Attenuation is %.2f\n',acs2)
figure, plot(f,sld1)
hold on
plot(f,sld2)
plot(f,acs1*f, 'k--')
plot(f,acs2*f, 'k--')
hold off
grid on,
xlim([0,max(f)]), ylim([0 15]),
xlabel('Frequency [MHz]')
ylabel('Att. [dB/cm]')
title('Mean SLD')
legend('Inc','Back')


%% RSLD
muB = 10.^(3:0.5:3.5);
muC = 10.^(1:2);
% muB = 10^ 3.5;
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'Isotropic RSLD')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;

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
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
% aSNR = 1; bSNR = 0.1;
aSNR = 1; bSNR = 0.5;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

%%
muB = 10.^(2.5:0.5:3);
muC = 10.^(0:2);
% muB = 10.^2.5;
% muC = 10^0.5;

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
        m,n,tol,mask(:),w);
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);
        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

figure('Units','centimeters', 'Position',[5 5 22 6]);
tl = tiledlayout(1,3, "Padding","tight");
title(tl,'RSLD - SWTV by British Columbia')
t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Weights')
c = colorbar;

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% Minimizing BS log ratio
muB = 10.^(3:0.5:4);
muC = 10.^(0:1:2);
% muB = 10^3;
% muC = 10^0;
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
        AttInterp = interp2(X,Z,BR,Xq,Zq);

        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end
%%
figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'RSLD with TV(B)+||C||_1')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD-TVL1, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD-TVL1, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS
muB = 10^3; muC = 10^0;

[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

%%
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;
% 
muB = 10.^(2.5:0.5:4);
muC = 10.^(0:0.5:2);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        toc

        BR = reshape(Bn*NptodB,m,n);
        CR = (reshape(Cn*NptodB,m,n));
%        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
        AttInterp = interp2(X,Z,BR,Xq,Zq);

        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

figure('Units','centimeters', 'Position',[5 5 22 6]);
tl = tiledlayout(1,3, "Padding","tight");
title(tl,'Weighted Fidelity and Regularization')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Weights')
c = colorbar;
%c.Label.String = 'BS log ratio [dB]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD-WFR, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD-WFR, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");

r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;


%%
save_all_figures_to_directory(figDir,['sim',num2str(iAcq),'Figure']);
% close all


end


%% Attenuation 

results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasBack)
disp(results2.biasBack)
disp(results3.biasBack)
disp(results4.biasBack)

disp('Bias Bottom')
disp(results1.biasInc)
disp(results2.biasInc)
disp(results3.biasInc)
disp(results4.biasInc)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)