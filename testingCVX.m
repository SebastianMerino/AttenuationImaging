clear,clc
close all
addpath('./functions_v7');

% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID316V2\06-08-2023-Generic'];

% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% Loading data
iAcq = 6;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);

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
tiledlayout(1,2);
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

%% Standard SLD with CVX
b = b(:);

cvx_begin
    variable u(2*m*n)
    minimize( norm( A * u - b, 2 ) )
cvx_end

BS = u(1:end/2); %CS = u(end/2+1:end);
BS = reshape(BS*8.686,m,n);    % [dB.cm^{-1}.MHz^{-1}]

figure('Units','centimeters', 'Position',[5 5 20 8]);
tiledlayout(1,2);
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
mu = 1e4*[0.4,1.2,3.6];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    tic
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
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
saveas(gcf,[figDir,'\equalizedT',num2str(iAcq),'.png'])

%% Regularized SLD with CVX -> ESTO NO SE PUEDE PIPIPIPI
% M = m; N = n;
% D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
% D(:,end) = [];
% D(M,M) = 0;
% Dx = kron(speye(N),D);
% 
% D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
% D(:,end) = [];
% D(N,N) = 0;
% Dy = kron(D,speye(M));
% 
% 
% b = (log(Sp) - log(Sd)) - (diffraction_compensation);
% 
% A1 = kron( 4*L*f , speye(m*n) );
% A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];
% 
% b = b(:);
% 
% cvx_begin
%     variable B(m*n)
%     variable C(m*n)
%     minimize( 0.5*norm( b - A1*B - A2*C, 2 ) + ...
%         norm( huber(square(Dx*B) + square(Dy*B) ) ,1) )
% cvx_end

%%

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];


% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
%mu = 1e4*[0.4,1.2,3.6];
mu = 4e3;
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    tic
    [Bn,Cn,Error] = AlterOpti_ADMM_v2(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));
    toc
    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,2);
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
saveas(gcf,[figDir,'\equalizedT',num2str(iAcq),'.png'])

