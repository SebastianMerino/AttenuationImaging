clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\Simulation_23_12_18'];
baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\process_simulation\23_12_31'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

figDir = [baseDir,'\fig\12-19'];
if (~exist(figDir,"dir")), mkdir(figDir); end

% groundTruthTop = [0.5,1,1,0.5,1,1];
% groundTruthBottom = [1,0.5,1,1,0.5,1];
groundTruthTop = [0.6,0.6,0.6,1.2,1.2,1.2];
groundTruthBottom = [1.2,1.2,1.2,0.6,0.6,0.6];
%% Loading data
for iAcq = 1:length(croppedFiles)
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
% [~,Z] = meshgrid(x_ACS,z_ACS);
% attIdeal = ones(size(Z));
% attIdeal(Z<=2) = groundTruthTop(iAcq);
% attIdeal(Z>2) = groundTruthBottom(iAcq);
[X,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
rInc = 0.8;
inclusion = (X.^2 + (Z-2).^2)<= rInc^2;
attIdeal(~inclusion) = groundTruthTop(iAcq);
attIdeal(inclusion) = groundTruthBottom(iAcq); %incl = bottom

figure('Units','centimeters', 'Position',[5 5 9 6]);
imagesc(x,z,Bmode,dynRange)
axis image
colormap(gray)
colorbar('westoutside')
title('Bmode')

figure('Units','centimeters', 'Position',[5 5 9 6]);
imagesc(x,z,attIdeal,attRange)
axis image
colormap(turbo)
c = colorbar('westoutside');
c.Label.String = 'dB/cm/MHz';
title('Ideal ACS')
%% Standard SLD
[u,~] = cgs(A'*A,A'*b(:),1e-6,20);
% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); CS = u(end/2+1:end);
BS = reshape(BS,m,n)*NptodB;
CS = reshape(CS,m,n)*NptodB;

figure('Units','centimeters', 'Position',[5 5 30 8]);
tl = tiledlayout(1,3);
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
c.Label.String = 'BS log ratio [dB]';

%% Spectrum
% [~,Z] = meshgrid(x_ACS,z_ACS);
% topLayer = Z < 1.9;
% bottomLayer = Z > 2.1;
% NpTodB = 20*log10(exp(1));
% 
% topSLD = squeeze(sum(sum(b.*topLayer,1),2))/sum(topLayer(:));
% slopeTop = f\topSLD;
% fprintf('Attenuation is %.2f\n',slopeTop*NpTodB)
% 
% bottomSLD = squeeze(sum(sum(b.*bottomLayer,1),2))/sum(bottomLayer(:));
% slopeBottom = f\bottomSLD;
% fprintf('Attenuation is %.2f\n',slopeBottom*NpTodB)
% plot(f,topSLD)
% hold on
% plot(f,bottomSLD)
% plot(f,slopeTop*f, 'k--')
% plot(f,slopeBottom*f, 'k--')
% hold off
% 
% legend('Top','Bottom')


%% RSLD
muB = 10.^(2.5:0.5:3.5);
muC = 10.^(0:1:3);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
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

axialTV = mean(BRopt,2);
%%
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
top = Zq < 1.9;
bottom = Zq > 2.1;

AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
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
muB = 10.^(2.5:0.5:3.5);
muC = 10.^(0.5:0.5:3);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
        m,n,tol,mask(:),w);
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
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

axialSWTV = mean(BRopt,2);
%%
AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% Minimizing BS log ratio
muB = 10.^(3:0.5:4);
muC = 10.^(0:0.5:2);
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);
        figure('Units','centimeters', 'Position',[5 5 15 6]);
        tl = tiledlayout(1,2, "Padding","tight");
        title(tl,'RSLD with TV(B)+||C||_1')
        
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR, attRange)
        colormap(t2,turbo)
        axis image
        title(['RSLD-TVL1, \mu=',num2str(muB(mmB),2)])
        c = colorbar;
        c.Label.String = 'Att. [db/cm/MHz]';
        
        t3 = nexttile; 
        imagesc(x_ACS,z_ACS,CR, bsRange)
        colormap(t3,parula)
        axis image
        title(['RSLD-TVL1, \mu=',num2str(muC(mmC),2)])
        c = colorbar;
        c.Label.String = 'BS log ratio [dB]';

        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
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

axialTVL1 = mean(BRopt,2);

%%
AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS
muB = 10^3; muC = 10^1;

[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 6;
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
muB = 10.^(3.5:0.5:4.5);
muC = 10.^(1:0.5:2);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        toc

        BR = reshape(Bn*NptodB,m,n);
        CR = (reshape(Cn*NptodB,m,n));
        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
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

axialWFR = mean(BRopt,2);

AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
% r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
%     groundTruthTop(iAcq),"omitnan") * 100;
% r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
%     groundTruthBottom(end),"omitnan") * 100;
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");

r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;

%%
figure('Units','centimeters', 'Position',[5 5 12 5])
plot(z_ACS, axialTV, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.4])
xlim([z_ACS(1) z_ACS(end)])
title('Axial profiles')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

%%
save_all_figures_to_directory(figDir,['sim',num2str(iAcq),'Figure']);
close all


end


%% Attenuation 

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

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)

% 
% figure('Units','centimeters', 'Position',[5 5 12 8])
% errorbar(results1.meanTop, results1.stdTop,'o', 'LineWidth',2)
% hold on,
% errorbar(results2.meanTop, results2.stdTop, 'x', 'LineWidth',2)
% errorbar(results3.meanTop, results3.stdTop, 's', 'LineWidth',2)
% plot(groundTruthTop, "_k"	, 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([0.3 1.6])
% grid on
% legend({'TV','TVL1','WFR'})
% title('Results top')
% xlabel('Target')
% ylabel('Attenuation [db/cm/MHz]')
% 
% figure('Units','centimeters', 'Position',[5 5 12 8])
% errorbar(results1.meanBottom, results1.stdBottom,'o', 'LineWidth',2)
% hold on,
% errorbar(results2.meanBottom, results2.stdBottom, 'x', 'LineWidth',2)
% errorbar(results3.meanBottom, results3.stdBottom, 's', 'LineWidth',2)
% plot(groundTruthBottom, '_k', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([0.3 1.6])
% grid on
% legend({'TV','TVL1','WFR'})
% title('Results bottom')
% xlabel('Target')
% ylabel('Attenuation [db/cm/MHz]')
% 
% %% CV
% figure('Units','centimeters', 'Position',[5 5 12 8])
% plot(results1.stdTop./results1.meanTop*100,'o', 'LineWidth',2)
% hold on,
% plot(results2.stdTop./results2.meanTop*100, 'x', 'LineWidth',2)
% plot(results3.stdTop./results3.meanTop*100, 's', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([0 15])
% grid on
% legend({'TV','TVL1','WFR'})
% title('Results top')
% xlabel('Target')
% ylabel('Coefficient of Variation [%]')
% 
% figure('Units','centimeters', 'Position',[5 5 12 8])
% plot(results1.stdBottom./results1.meanBottom*100,'o', 'LineWidth',2)
% hold on,
% plot(results2.stdBottom./results2.meanBottom*100, 'x', 'LineWidth',2)
% plot(results3.stdBottom./results3.meanBottom*100, 's', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([0 15])
% grid on
% legend({'TV','TVL1','WFR'})
% title('Results bottom')
% xlabel('Target')
% ylabel('Coefficient of Variation [%]')
% 
% %% Mean Percentage error
% figure('Units','centimeters', 'Position',[5 5 12 8])
% plot(results1.MPETop, 'o', 'LineWidth',2)
% hold on,
% plot(results2.MPETop, 'x', 'LineWidth',2)
% plot(results3.MPETop, 's', 'LineWidth',2)
% yline(0, 'k--', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([-10 50])
% grid on
% legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
% title('Results top')
% xlabel('Target')
% ylabel('MPE %')
% 
% figure('Units','centimeters', 'Position',[5 5 12 8])
% plot(results1.MPEBottom, 'o', 'LineWidth',2)
% hold on,
% plot(results2.MPEBottom, 'x', 'LineWidth',2)
% plot(results3.MPEBottom, 's', 'LineWidth',2)
% yline(0, 'k--', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([-10 50])
% grid on
% legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
% title('Results bottom')
% xlabel('Target')
% ylabel('MPE %')
% 
% %% cnr
% figure('Units','centimeters', 'Position',[5 5 12 8])
% plot(results1.cnr, 'o', 'LineWidth',2)
% hold on,
% plot(results2.cnr, 'x', 'LineWidth',2)
% plot(results3.cnr, 's', 'LineWidth',2)
% hold off
% xlim([0.5 6.5])
% ylim([-0.1 20])
% grid on
% legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
% title('Contrast-to-noise ratio')
% xlabel('Target')
% ylabel('CNR')
% 
% 
% %%
% save_all_figures_to_directory(figDir,'figure');

