clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\layeredNew'];
% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\layeredNew'];
% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\layered_14_11_23'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layered_14_11_23'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

figDir = [baseDir,'\fig\22-11'];
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
[~,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
attIdeal(Z<=2) = groundTruthTop(iAcq);
attIdeal(Z>2) = groundTruthBottom(iAcq);

% Creating meshgrid for interpolation and ROIs
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
top = Zq < 1.9;
bottom = Zq > 2.1;

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


%% Minimizing BS log ratio and WEIGHTS - VERSION 1 
% First estimation
muB = 10^4; muC = 10^1.5;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 6;
w = 1./((bscMap/ratioCutOff).^2 + 1);
%figure, imagesc(w)

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muB = 10.^(3.5:0.5:5);
muC = 10.^(1:0.5:3);
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
MetricsWFRv1(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS - VERSION 2
% First estimation
muB = 10^4; muC = 10^1.5;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 6;
w = 1./((bscMap/ratioCutOff).^4 + 1);
% figure, imagesc(w)

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muB = 10.^(3.5:0.5:5);
muC = 10.^(1:0.5:3);
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
MetricsWFRv2(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS - VERSION 3
% First estimation
muB = 10^4; muC = 10^1.5;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 6;
order = 5;
reject = 0.1;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
% figure, imagesc(w)

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muB = 10.^(3.5:0.5:5);
muC = 10.^(1:0.5:3);
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
MetricsWFRv3(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS - VERSION 4
% First estimation
muB = 10^4; muC = 10^1.5;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

ratioCutOff = 6;
order = 5;
reject = 0.1;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,3);
%figure, imagesc(w)

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muB = 10.^(3.5:0.5:5);
muC = 10.^(1:0.5:3);
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
MetricsWFRv4(iAcq) = r;

%%
save_all_figures_to_directory(figDir,['sim',num2str(iAcq),'Figure']);
close all


end


%% Attenuation 

results1 = struct2table(MetricsWFRv1);
results2 = struct2table(MetricsWFRv2);
results3 = struct2table(MetricsWFRv3);
results4 = struct2table(MetricsWFRv4);

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

