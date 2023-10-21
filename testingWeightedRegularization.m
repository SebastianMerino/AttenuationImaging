clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets' ...
%     '\Attenuation\Timana\data'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\Timana\ref'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig\20-10'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% Loading data
for iAcq = 1:8
iAcq = 7;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);
attRange = [0.4,1.1];
bsRange = [-2 2];
%% Spectrum
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

%% Standard SLD
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

[u,~] = cgs(A'*A,A'*b(:),1e-6,20);

% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); %CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);

figure('Units','centimeters', 'Position',[5 5 20 8]);
tl = tiledlayout(1,2);
title(tl,'Standard RSLD')
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BS, attRange)
colormap(t2,turbo)
axis equal tight
title('SLD')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';


%% RSLD
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0.5,1.5,3);
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
    
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t3,gray)
    colorbar(t3,'westoutside')
    title('Bmode')
    
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

%% British Columbia Approach
envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        
        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
a = 1; b = 0.1;
desvMin = 15;
w = a./(1 + exp(b.*(desvSNR - desvMin)));

% % Weights
% figure('Units','centimeters', 'Position',[5 5 30 8]),
% tl = tiledlayout(1,3);
% title(tl,{'Weights proposed by BC',''});
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% colormap(t1,gray)
% colorbar
% axis equal
% xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
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
% t3 = nexttile;
% imagesc(x_ACS,z_ACS,w)
% colormap(t3,parula)
% colorbar;
% axis image
% title('Weights')

% RSLD ANISOTROPIC AND BS WEIGHTED
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
mu2 = mu/100;
for mm = 1:length(mu)
    tic
    [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),mu(mm),mu2(mm),...
        m,n,tol,mask(:),w);
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 12]);
tl = tiledlayout(2,size(BR,3)+1);
title(tl,'British Columbia Approach')
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
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar(t3,'westoutside')
axis image
title('Weights')

for ii = 1:size(BR,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
    colormap(t2,parula)
    axis equal tight
    title(['RSLD, \mu=',num2str(mu2(ii),2)])
end
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';

%% Minimizing BS log ratio
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
mu2 = logspace(0.5,1.5,3);
BR = zeros(m,n,length(mu2));
CR = zeros(m,n,length(mu2));
for mm = 1:length(mu)
    for mm2 = 1:length(mu2)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),mu(mm),mu2(mm2),m,n,tol,mask(:));
        toc
        BR(:,:,mm2) = (reshape(Bn*8.686,m,n));
        CR(:,:,mm2) = (reshape(Cn,m,n));
    end
    
    % Plotting
    figure('Units','centimeters', 'Position',[5 5 30 12]);
    tl = tiledlayout(2,size(BR,3)+1);
    title(tl,'RSLD with isotropic TV and Tikhonov reg.')
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
    
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)]),
    ylim([z_ACS(1) z_ACS(end)]),
    colormap(t3,gray)
    colorbar(t3,'westoutside')
    title('Bmode')
    
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


%% PROPOSED WEIGHTS
envelope = abs(hilbert(sam1));

BSratio = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        BSratio(ii,jj) = mean(sub_block_p(:))/mean(sub_block_d(:));
    end
end
w = 1./((db(BSratio)/3).^2 + 1);

% Weights
figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,'Proposed Weights');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,db(BSratio))
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
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = logspace(2.5,3.5,3);
%mu2 = mu/100;
mu2 = [3.2 3.2 3.2];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    tic
    [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,mu(mm),mu2(mm),m,n,tol,mask(:),w);
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 12]);
tl = tiledlayout(2,size(BR,3)+1);
title(tl,'TV, Tikhonov reg and weights')
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
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile;
imagesc(x_ACS,z_ACS,w)
colormap(t3,parula)
colorbar(t3,'westoutside')
axis image
title('Weights')

for ii = 1:size(BR,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,CR(:,:,ii), bsRange)
    colormap(t2,parula)
    axis equal tight
    title(['RSLD, \mu=',num2str(mu2(ii),2)])
end
c = colorbar;
c.Label.String = 'BS log ratio (a.u.)';

%%
targetDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);
close all
end