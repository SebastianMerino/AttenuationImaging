% Tests new optimization method and second derivatives

setup,

dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_04_inc';
refDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\24_04_25_ref';
resultsDir = 'C:\Users\sebas\Documents\Data\Attenuation\JournalResults\sim_inc';
% dataDir = 'P:\smerino\simulation_acs\rf_data\24_04_04_inc';
% refDir = 'P:\smerino\simulation_acs\rf_data\24_04_25_ref';
% resultsDir = 'P:\smerino\UFFC2024results\simulation';

[~,~] = mkdir(resultsDir);
targetFiles = dir([dataDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);
tableName = 'simuInc.xlsx';

%%
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6; % original 3.3-8.7s
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% New simu
referenceAtt    = 0.6;
groundTruthBack = [0.5,0.5,0.5];
groundTruthInc = [1,1,1];

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting
dynRange = [-40,0];
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

iAcq = 2;

%% Setting up

%for iAcq = 1:3
% Regularization
switch iAcq
    case 1
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBswift = 10^3.5; muCswift = 10^2;
    case 2
        muBtv = 10^2.5; muCtv = 10^2.5;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
    case 3
        muBtv = 10^3.5; muCtv = 10^3.5;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3.5; muCtvl1 = 10^1;
        muBswift = 10^3.5; muCswift = 10^1;
end

load(fullfile(dataDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
z = z-3.5*medium.sound_speed_ref/6.66e6*100*0.5;

sam1 = rf(:,:,1);
Bmode = db(hilbert(sam1));
dynRange = [-50,0];

%% Cropping and finding sample sizes
% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);
Bmode = Bmode(ind_z,ind_x);
Bmode = Bmode - max(Bmode(:));

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
nz = 2*round(blocksize*wl/dz /2 * ratio_zx); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

%% Spectrum
% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation
% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
windowing = windowing*ones(1,nx);

% For looping
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p);
Sd_ref = zeros(m,n,p);
compensation = zeros(m,n,p,Nref);

for iRef = 1:Nref %Nref
    load(fullfile(refDir,refFiles(iRef).name),"rf","medium");
    acs_mean = medium.alpha_coeff(1,1);
    att_ref = acs_mean*(f.^medium.alpha_power)/NptodB;
    att_ref_map = repmat(reshape(att_ref,[1 1 p]),m,n,1);

    samRef = rf;
    samRef = samRef(ind_z,ind_x); % Cropping
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:) = (tempSp(rang));
            Sd_ref(ii,jj,:) = (tempSd(rang));
        end
    end
    compensation(:,:,:,iRef) = log(Sp_ref) - log(Sd_ref) - 4*L*att_ref_map;
end

compensation = mean(compensation,4);
% compensation = repmat(mean(compensation,3),1,1,p);

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

%% Setting up
% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating masks and ideal map
rInc = 0.7; c1x = 0; c1z = 2;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
inclusion = (Xq.^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = (Xq.^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;
attIdeal = ones(size(Xq))*groundTruthBack(iAcq);
attIdeal((Xq.^2 + (Zq-c1z).^2)<= rInc^2) = groundTruthInc(iAcq);

inclusionACS = (X.^2 + (Z-c1z).^2)<= rInc^2;
attIdealACS = ones(size(X))*groundTruthBack(iAcq);
attIdealACS(inclusionACS) = groundTruthInc(iAcq); %incl = inclusion

%% TV
muBtv = 10^3;
muCtv = 10^3;
tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 22 6]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('RSLD')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRTV, bsRange)
colormap(t1,"parula")
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')


%% Second derivatives
w = ones(m,n); wExt = w;
Wdiag = spdiags(w(:),0,m*n,m*n);
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);


% muB = 10^1.5;
% lambdaB = 0;
% muC = 10^0.5;
% lambdaC = 0;

% muB = 0;
% lambdaB = 10^2;
% muC = 0;
% lambdaC = 10^1;

muB = 10^1.5;
lambdaB = 10^1.5;
muC = 10^0.5;
lambdaC = 10^0.5;

D = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1);
D(:,end) = [];
D(m,m) = 0;
Dy = kron(speye(n),D);

D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1);
D(:,end) = [];
D(n,n) = 0;
Dx = kron(D,speye(m));

Dyy = Dy*Dy;
Dxx = Dx*Dx;
Dxy = Dx*Dy;
DB = [muB*Wdiag*Dx;muB*Wdiag*Dy;...
    lambdaB*Wdiag*Dxx;lambdaB*Wdiag*Dyy;sqrt(2)*lambdaB*Wdiag*Dxy];
DC = [muC*speye(m*n);lambdaC*Dx;lambdaC*Dy];
Tik = blkdiag(DB,DC);


u = pcg(A'*A + Tik'*Tik,A'*b(:),tol,2000);
Bn = u(1:m*n);
Cn = u(1+m*n:end);

BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 22 6]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRTV, bsRange)
colormap(t1,"parula")
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

%% Second derivatives + norm 1
w = ones(m,n); wExt = w;
Wdiag = spdiags(w(:),0,m*n,m*n);
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);


% muB = 10^2.4;
% lambdaB = 0;
% muC = 10^1;
% lambdaC = 0;
% muB = 0;
% lambdaB = 10^2.5;
% muC = 0;
% lambdaC = 10^1.5;

muB = 10^1.5;
lambdaB = 10^2.5;
muC = 10^1;
lambdaC = 10^1.5;

D = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1);
D(:,end) = [];
D(m,m) = 0;
Dy = kron(speye(n),D);

D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1);
D(:,end) = [];
D(n,n) = 0;
Dx = kron(D,speye(m));

Dyy = Dy*Dy;
Dxx = Dx*Dx;
Dxy = Dx*Dy;
DB = [muB*Wdiag*Dx;muB*Wdiag*Dy;...
    lambdaB*Wdiag*Dxx;lambdaB*Wdiag*Dyy;sqrt(2)*lambdaB*Wdiag*Dxy];
DC = [muC*speye(m*n);lambdaC*Dx;lambdaC*Dy];
Tik = blkdiag(DB,DC);

% Optimization with Lavarello's implementation
% Params.alpha2 = 1; Params.tolerance = tol/1e5;
% Params.beta = 0.5; Params.k = 1;
% Params.operator = 'L'; Params.L = Tik;
% u = optimTikhonovReg(W*A,W*b(:),Params);

% My implementation
tol = 1e-5; rho = 0.01; maxIte = 500;
u = optimAdmmTikhonov(A, b(:), Tik, 0.1, rho, maxIte, tol);
Bn = u(1:m*n);
Cn = u(1+m*n:end);

BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 22 6]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRTV, bsRange)
colormap(t1,"parula")
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')
%% Weight map
muBnew = 10^1.5; 
muCnew = 10^0.2;
D = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1);
D(:,end) = [];
D(m,m) = 0;
Dy = kron(speye(n),D);
% Dy = Wdiag*Dy;

D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1);
D(:,end) = [];
D(n,n) = 0;
Dx = kron(D,speye(m));
% Dx = Wdiag*Dx;

D = [Dx;Dy];
D = blkdiag(muBnew*D,muCnew*speye(m*n));
% D = blkdiag(muBnew*D,muCnew*D);

% Params.alpha2 = 1; Params.tolerance = tol/1e5;
% Params.beta = 0; Params.k = 2;
% Params.operator = 'L'; Params.L = D;
% u = optimTikhonovReg(A,b(:),Params);
u = pcg(A'*A + D'*D,A'*b(:),tol,2000);
Cn = u(1+m*n:end);
CRTV = reshape(Cn*NptodB,m,n);

% Weight map
ratioCutOff = 7;
w = (1-reject)*(abs(CRTV)<ratioCutOff)+reject;
wExt = movmin(w,extension);

figure,
imagesc(x_ACS,z_ACS,w)
axis image

%% Second derivatives + weights
muB = 10^2;
lambdaB = 10^2;
muC = 10^0.5;
lambdaC = 10^0;

D = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1);
D(:,end) = [];
D(m,m) = 0;
Dy = kron(speye(n),D);

D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1);
D(:,end) = [];
D(n,n) = 0;
Dx = kron(D,speye(m));

Dyy = Dy*Dy;
Dxx = Dx*Dx;
Dxy = Dx*Dy;
DB = [muB*Wdiag*Dx;muB*Wdiag*Dy;...
    lambdaB*Wdiag*Dxx;lambdaB*Wdiag*Dyy;lambdaB*Wdiag*Dxy];
DC = [muC*speye(m*n);lambdaC*Dx;lambdaC*Dy];
Tik = blkdiag(DB,DC);

W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
% Params.alpha2 = 1; Params.tolerance = tol/1e5;
% Params.beta = 0; Params.k = 2;
% Params.operator = 'L'; Params.L = Tik;
% u = optimTikhonovReg(W*A,W*b(:),Params);
WA = W*A; Wb = W*b(:);
u = pcg(WA'*WA + Tik'*Tik,WA'*Wb(:),tol,2000);
Bn = u(1:m*n);
Cn = u(1+m*n:end);

BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 22 6]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('Weighted')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRTV, bsRange)
colormap(t1,"parula")
axis image
title('Weighted')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

%% Second derivatives + weights + norm 1

% muB = 10^2.4;
% lambdaB = 0;
% muC = 10^1;
% lambdaC = 0;
% muB = 0;
% lambdaB = 10^2.5;
% muC = 0;
% lambdaC = 10^1.5;

muB = 10^2;
lambdaB = 10^2;
muC = 10^1;
lambdaC = 10^1.5;

D = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1);
D(:,end) = [];
D(m,m) = 0;
Dy = kron(speye(n),D);

D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1);
D(:,end) = [];
D(n,n) = 0;
Dx = kron(D,speye(m));

Dyy = Dy*Dy;
Dxx = Dx*Dx;
Dxy = Dx*Dy;
DB = [muB*Wdiag*Dx;muB*Wdiag*Dy;...
    lambdaB*Wdiag*Dxx;lambdaB*Wdiag*Dyy;sqrt(2)*lambdaB*Wdiag*Dxy];
DC = [muC*speye(m*n);lambdaC*Dx;lambdaC*Dy];
Tik = blkdiag(DB,DC);

WA = W*A; Wb = W*b(:);

% My implementation
tol = 1e-5; rho = 0.1; maxIte = 500;
u = optimAdmmTikhonov(WA, Wb, Tik, 0.1, rho, maxIte, tol);
Bn = u(1:m*n);
Cn = u(1+m*n:end);

BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 22 6]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
title('Ideal')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRTV, bsRange)
colormap(t1,"parula")
axis image
title('Test')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')