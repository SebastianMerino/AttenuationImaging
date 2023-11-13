clear,clc
close all
addpath('./functions_v7');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layeredNew'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

figDir = [baseDir,'\fig\10-11'];
if (~exist(figDir,"dir")), mkdir(figDir); end

groundTruthTop = [0.5,1,1,0.5,1,1];
groundTruthBottom = [1,0.5,1,1,0.5,1];

%% Loading data
for iAcq = 1:3
iAcq = 5;
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");

dynRange = [-50,0];
attRange = [0,1.5];
bsRange = [-2 2];


%% Standard SLD
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
[u,~] = cgs(A'*A,A'*b(:),1e-6,20);
% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);
CS = reshape(CS,m,n);

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
c.Label.String = 'BS log ratio (a.u.)';

%% 
[~,Z] = meshgrid(x_ACS,z_ACS);
topLayer = Z < 2.2;
bottomLayer = Z > 2.4;
NpTodB = 20*log10(exp(1));

topSLD = squeeze(sum(sum(b.*topLayer,1),2))/sum(topLayer(:));
slopeTop = f\topSLD;
fprintf('Attenuation is %.2f\n',slopeTop*NpTodB)

bottomSLD = squeeze(sum(sum(b.*bottomLayer,1),2))/sum(bottomLayer(:));
slopeBottom = f\bottomSLD;
fprintf('Attenuation is %.2f\n',slopeBottom*NpTodB)
plot(f,topSLD)
hold on
plot(f,bottomSLD)
plot(f,slopeTop*f, 'k--')
plot(f,slopeBottom*f, 'k--')
hold off

legend('Top','Bottom')

%% RSLD
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

%attRange = [0 1];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0,2,3);
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
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
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    nexttile;
    axis off
    
    for ii = 1:size(BR,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end
%%
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR(:,:,2),Xq,Zq);

top = Zq < 2.3;
bottom = Zq > 2.5;

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

% figure,
% imagesc(x_ACS,z_ACS,SNR)
% colorbar

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
a = 1; b = 0.1;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0,2,3);
BBC = zeros(m,n,length(mu));
CBC = zeros(m,n,length(mu));

for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),mu(mm),mu2(mm2),...
        m,n,tol,mask(:),w);
        toc
        BBC(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CBC(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BBC,3)+1);
    title(tl,'British Columbia Approach')
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BBC,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BBC(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    t3 = nexttile;
    imagesc(x_ACS,z_ACS,w,[0 1])
    axis image
    colormap(t3,parula)
    colorbar(t3,'westoutside')
    title('Weights')
    
    for ii = 1:size(BBC,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CBC(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end

%%
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BBC(:,:,2),Xq,Zq);

top = Zq < 2.3;
bottom = Zq > 2.5;

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

%% Minimizing BS log ratio
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(3,4,3);
mu2 = logspace(0,2,3);
BR2 = zeros(m,n,length(mu2));
CR2 = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR2(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR2(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR2,3)+1);
    title(tl,{'RSLD with isotropic TV and Tikhonov reg.',''})
    t1 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t1,gray)
    colorbar(t1,'westoutside')
    title('Bmode')
    
    for ii = 1:size(BR2,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,BR2(:,:,ii), attRange)
        colormap(t2,turbo)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    nexttile; axis off
    
    for ii = 1:size(BR2,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CR2(:,:,ii), bsRange)
        colormap(t2,parula)
        axis equal tight
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar;
    c.Label.String = 'BS log ratio (a.u.)';
end

%%
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR2(:,:,2),Xq,Zq);

top = Zq < 2.3;
bottom = Zq > 2.5;

r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
    groundTruthTop(iAcq),"omitnan") * 100;
r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
    groundTruthBottom(end),"omitnan") * 100;
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVTik(iAcq) = r;

%% NEW WEIGHTS
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );


% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
mu = 1e3;
mu2 = 1;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu,mu2,m,n,tol,mask(:));
bscMap = (reshape(Cn,m,n));

logBscRatio = bscMap*log10(exp(1))*20;
w = 1./((logBscRatio/6).^2 + 1);

% [~,Z] = meshgrid(x_ACS,z_ACS);
% w = double(Z < 2.3 | Z > 2.8)*0.9 + 0.1;
% w(Z>2.8) = 0.6;

%%

figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'New Weights',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,logBscRatio, [-20 20])
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('BS ratio')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar;
axis image
title('Weights')

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

mask = ones(m,n,p);
mu = logspace(3.5,4.5,3);
mu2 = logspace(0.5,2.5,3);
BRW = zeros(m,n,length(mu2));
CRW = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu(mm),mu2(mm2),m,n,tol,mask(:),w);
        toc
        BRW(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CRW(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BRW,3)+1);
    title(tl,{'TV, Tikhonov reg and weights',''})
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
        title(['RSLD, \mu=',num2str(mu(mm),2)])
    end
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    
    t3 = nexttile;
    imagesc(x_ACS,z_ACS,w,[0 1])
    axis image
    colormap(t3,parula)
    colorbar(t3,'westoutside')
    title('Weights')
    
    for ii = 1:size(BRW,3)
        t2 = nexttile; 
        imagesc(x_ACS,z_ACS,CRW(:,:,ii), bsRange)
        colormap(t2,parula)
        axis image
        title(['RSLD, \mu=',num2str(mu2(ii),2)])
    end
    c = colorbar(t2);
    c.Label.String = 'BS log ratio (a.u.)';
end

%%
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRW(:,:,2),Xq,Zq);

top = Zq < 2.3;
bottom = Zq > 2.5;

r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.MPETop = mean( (AttInterp(top) - groundTruthTop(iAcq)) /...
    groundTruthTop(iAcq),"omitnan") * 100;
r.MPEBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)) /...
    groundTruthBottom(end),"omitnan") * 100;
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTVTik(iAcq) = r;

%%
attTV = mean(BR(:,:,2),2);
attSWTV = mean(BBC(:,:,2),2);
attTVTik = mean(BR2(:,:,2),2);
attSWTVTik = mean(BRW(:,:,2),2);

bscTV = mean(CR(:,:,2),2);
bscSWTV = mean(CBC(:,:,2),2);
bscTVTik = mean(CR2(:,:,2),2);
bscSWTVTik = mean(CRW(:,:,2),2);

figure('Units','centimeters', 'Position',[5 5 20 15]);
tl = tiledlayout(2,1);
title(tl,'Mean vertical profile')

nexttile
plot(z_ACS,[attTV,attSWTV,attTVTik,attSWTVTik])
title('Attenuation')
grid on
xlim([z_ACS(1),z_ACS(end)])
xlabel('Depth [cm]'), ylabel('ACS [dB/cm/MHz]')
legend({'TV','SWTV','TV+Tik', 'SWTV+SWTik'}, 'Location','northeastoutside')

nexttile
plot(z_ACS,[bscTV,bscSWTV,bscTVTik,bscSWTVTik])
title('Backscatter')
grid on
xlim([z_ACS(1),z_ACS(end)])
xlabel('Depth [cm]'), ylabel('BCS log ratio [a.u.]')
legend({'TV','SWTV','TV+Tik', 'SWTV+SWTik'}, 'Location','northeastoutside')


% %%
newDir = fullfile(figDir,croppedFiles(iAcq).name(1:end-4));
if(~exist(newDir,"dir")), mkdir(newDir); end
save_all_figures_to_directory(newDir);
close all


end


%% Attenuation 

results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsSWTVTik);

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(results1.meanTop, results1.stdTop,'o', 'LineWidth',2)
hold on,
errorbar(results2.meanTop, results2.stdTop, 'x', 'LineWidth',2)
errorbar(results3.meanTop, results3.stdTop, 's', 'LineWidth',2)
plot(groundTruthTop, "_k"	, 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0.3 1.6])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'})
title('Results top')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(results1.meanBottom, results1.stdBottom,'o', 'LineWidth',2)
hold on,
errorbar(results2.meanBottom, results2.stdBottom, 'x', 'LineWidth',2)
errorbar(results3.meanBottom, results3.stdBottom, 's', 'LineWidth',2)
plot(groundTruthBottom, '_k', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0.3 1.6])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'})
title('Results bottom')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

%% CV
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(results1.stdTop./results1.meanTop*100,'o', 'LineWidth',2)
hold on,
plot(results2.stdTop./results2.meanTop*100, 'x', 'LineWidth',2)
plot(results3.stdTop./results3.meanTop*100, 's', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0 15])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'})
title('Results top')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(results1.stdBottom./results1.meanBottom*100,'o', 'LineWidth',2)
hold on,
plot(results2.stdBottom./results2.meanBottom*100, 'x', 'LineWidth',2)
plot(results3.stdBottom./results3.meanBottom*100, 's', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([0 15])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'})
title('Results bottom')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

%% Mean Percentage error
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(results1.MPETop, 'o', 'LineWidth',2)
hold on,
plot(results2.MPETop, 'x', 'LineWidth',2)
plot(results3.MPETop, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([-10 50])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
title('Results top')
xlabel('Target')
ylabel('MPE %')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(results1.MPEBottom, 'o', 'LineWidth',2)
hold on,
plot(results2.MPEBottom, 'x', 'LineWidth',2)
plot(results3.MPEBottom, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([-10 50])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
title('Results bottom')
xlabel('Target')
ylabel('MPE %')

%% cnr
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(results1.cnr, 'o', 'LineWidth',2)
hold on,
plot(results2.cnr, 'x', 'LineWidth',2)
plot(results3.cnr, 's', 'LineWidth',2)
hold off
xlim([0.5 3.5])
ylim([-10 20])
grid on
legend({'TV','TV+Tik','SWTV+SWTik'}, 'Location','northeast')
title('Contrast-to-noise ratio')
xlabel('Target')
ylabel('CNR')


%%
save_all_figures_to_directory(figDir);

