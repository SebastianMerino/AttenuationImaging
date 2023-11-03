clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\CUELLO#3'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end

figDir = [baseDir,'\fig\01-11'];
if (~exist(figDir,"dir")), mkdir(figDir); end

% CASOS INTERESANTES: 2,4,6,8,9,13

%% Loading data

iAcq = 6;
fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
load(fullfile(croppedDir,croppedFiles(iAcq).name));
dynRange = [-40,-5];
attRange = [0.4,1.6];
%attRange = [0,1]; % Just for 13 acq
bsRange = [-2 2];

switch iAcq
    case 2
        muTV = 1e3; mu2TV = 10;
        muWTV = 3.2e2; mu2WTV = 10;
        muTik = 1.5E3; mu2Tik = 10;
        muWTik = 1.5E3; mu2WTik = 10;
    case 4
        muTV = 1e3; mu2TV = 10;
        muWTV = 3.2e2; mu2WTV = 10;
        muTik = 1.5E3; mu2Tik = 3.2;
        muWTik = 1.5E3; mu2WTik = 3.2;
    case 6
        muTV = 1e3; mu2TV = 10;
        muWTV = 3.2e2; mu2WTV = 10;
        muTik = 1.5E3; mu2Tik = 10;
        muWTik = 1.5E3; mu2WTik = 10;
    case 8
        muTV = 1e3; mu2TV = 10;
        muWTV = 3.2e2; mu2WTV = 10;
        muTik = 1.5E3; mu2Tik = 10;
        muWTik = 1.5E3; mu2WTik = 10;
    otherwise
        muTV = 1e3; mu2TV = 10;
        muWTV = 3.2e2; mu2WTV = 10;
        muTik = 1.5E3; mu2Tik = 10;
        muWTik = 1.5E3; mu2WTik = 10;
end
%% RSLD
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muTV,mu2TV,m,n,tol,mask(:));
toc
BR = (reshape(Bn*8.686,m,n));
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
a = 1; b = 0.1;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'Weights by SNR',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,db(SNR), [-1 7])
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')

t3 = nexttile;
imagesc(x_ACS,z_ACS,w,[0 1])
colormap(t3,parula)
colorbar;
axis image
title('Weights')

% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
tic
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muWTV,mu2WTV,m,n,tol,mask(:),w);
toc
BRBC = (reshape(Bn*8.686,m,n));
CRBC = (reshape(Cn,m,n));

%% Minimizing BS log ratio
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

tic
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muTik,mu2Tik,m,n,tol,mask(:));
toc
BRTik = (reshape(Bn*8.686,m,n));
CRTik = (reshape(Cn,m,n));

%% NEW WEIGHTS
b = (log(Sp) - log(Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );


% Regularization: Au = b
tol = 1e-3;
mask = ones(m,n,p);
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muWTik,mu2WTik,m,n,tol,mask(:));
bscMap = (reshape(Cn,m,n));

logBscRatio = bscMap*log10(exp(1))*20;
w = 1./((logBscRatio/10).^2 + 1);


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
imagesc(x_ACS,z_ACS,w, [0 1])
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

tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muWTik,mu2WTik,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*8.686,m,n));
CRWTik = (reshape(Cn,m,n));

%% Plotting three results
figure('Units','centimeters', 'Position',[5 5 30 10]);
tiledlayout(2,5);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BR, attRange)
colormap(t2,turbo)
axis equal tight
title(['TV, \mu=',num2str(muTV,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRBC, attRange)
colormap(t2,turbo)
axis equal tight
title(['SWTV, \mu=',num2str(muWTV,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTik, attRange)
colormap(t2,turbo)
axis equal tight
title(['TV and Tik, \mu=',num2str(muTik,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWTik, attRange)
colormap(t2,turbo)
axis equal tight
title(['SWTV and Tik, \mu=',num2str(muWTik,2)])

c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

nexttile;
axis off

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CR, bsRange)
colormap(t2,parula)
axis equal tight
title(['TV, \mu=',num2str(mu2TV,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRBC, bsRange)
colormap(t2,parula)
axis equal tight
title(['SWTV, \mu=',num2str(mu2WTV,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTik, bsRange)
colormap(t2,parula)
axis equal tight
title(['TV and Tik, \mu=',num2str(mu2Tik,2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRWTik, bsRange)
colormap(t2,parula)
axis equal tight
title(['SWTV and Tik, \mu=',num2str(mu2WTik,2)])

c = colorbar;
c.Label.String = 'Bsc. log ratio [a.u.]';

%%
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure,
[~,~,hColor] = imOverlayInterp(BmodeFull,BR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';

figure,
[~,~,hColor] = imOverlayInterp(BmodeFull,BRWTik,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';

%%
newDir = fullfile(figDir,croppedFiles(iAcq).name(1:end-4));
if(~exist(newDir,"dir")), mkdir(newDir); end
save_all_figures_to_directory(newDir);
% close all
