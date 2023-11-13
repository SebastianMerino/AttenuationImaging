clear,clc
close all
addpath('./functions_v7');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\layeredNew2'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

%%
iAcq = 4;
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
load(fullfile(baseDir,'raw',croppedFiles(iAcq).name),"medium");

b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );


dynRange = [-50,0];
attRange = [0,1.5];
bsRange = [-2 2];

[~,Z] = meshgrid(x_ACS,z_ACS);
if iAcq == 3
    topLayer = Z < 1.7;
    bottomLayer = Z > 2.3;
    alphaTop = mean(medium.alpha_coeff(1:500,:),"all")
    alphaBottom = mean(medium.alpha_coeff(600:end,:),"all")
    stdTop = std(medium.density(1:500,:),[],"all")/1000
    stdBottom = std(medium.density(600:end,:),[],"all")/1000
elseif iAcq == 4
    topLayer = Z < 2.2;
    bottomLayer = Z > 2.8;
    alphaTop = mean(medium.alpha_coeff(1:650,:),"all")
    alphaBottom = mean(medium.alpha_coeff(700:end,:),"all")
    stdTop = std(medium.density(1:650,:),[],"all")/1000
    stdBottom = std(medium.density(700:end,:),[],"all")/1000
end
NpTodB = 20*log10(exp(1));

topSLD = squeeze(sum(sum(b.*topLayer,1),2))/sum(topLayer(:));
slopeTop = f\topSLD;
fprintf('ACS Top: %.2f\n',slopeTop*NpTodB)

bottomSLD = squeeze(sum(sum(b.*bottomLayer,1),2))/sum(bottomLayer(:));
slopeBottom = f\bottomSLD;
fprintf('ACS Bottom: %.2f\n',slopeBottom*NpTodB)
figure,
plot(f,topSLD)
hold on
plot(f,bottomSLD)
plot(f,slopeTop*f, 'k--')
plot(f,slopeBottom*f, 'k--')
hold off

legend('Top','Bottom')
%%
axialSLD = squeeze(mean(b,2));
surf(f,z_ACS,movmean(axialSLD,10))
%plot(f,movmean(axialSLD,20)')
%% RSLD
% A = [A1 A2];

%attRange = [0 1];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(1,2,3)*10;
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
mu2 = logspace(1,2,3);
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

mask = ones(m,n,p);
mu = logspace(3,4,3);
mu2 = logspace(1,2,3);
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