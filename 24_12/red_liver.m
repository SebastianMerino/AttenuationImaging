setup,

dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\liver_RED';
refDir = fullfile(dataDir,'ref');

refFiles = dir(fullfile(refDir,'*.mat'));

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
dynRange = [-50,0];
attRange = [0.3,1.5];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 1.7; z_sup = 4.5;


%% Setting up
load(fullfile(dataDir,'rf_qus_livernew202410_AC_test'));

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
xBm = x;
zBm = z;

sam1 = rf(:,:,1);
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

%% Cropping and finding sample sizes
% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);
% Bmode = Bmode(ind_z,ind_x);

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
NFFT = 2^(nextpow2(nz/2));
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
    load(fullfile(refDir,refFiles(iRef).name),"rf");
    acs_mean = 0.4;
    att_ref = acs_mean*(f)/NptodB;
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
[Xq,Zq] = meshgrid(xBm,zBm);
attIdeal = interp2(kgrid.y*100, kgrid.x*100 - kgrid.x(1)*100, attenuation_map, Xq, Zq);

%% TV
muBtv = 10^3;
muCtv = 10^3;
tic
[Bn,Cn,ite] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
exTime = toc;
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);


roi = ind_x.*ind_z';
figure('Units','centimeters', 'Position',[5 5 24 8]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(xBm,zBm,Bmode,dynRange)
axis image
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(xBm,zBm,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
axis image
title('Ideal')

t3 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,BRTV,dynRange,attRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

t4 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,CRTV,dynRange,bsRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

colormap(t1,gray)
colormap(t3,turbo)
colormap(t4,parula)

%% SWIFT
muBswift = 10^3; muCswift = 10^2;

% First iteration
tic
[Bn,Cn,ite] = optimSwift(A1,A2,b(:),muBswift,muCswift,...
    m,n,tol,ones(m,n));
exTime = toc;
fprintf('\nExecution time: %.4f\n',exTime)
fprintf('Number of iterations: %d\n',ite)

BSWIFT = reshape(Bn*NptodB,m,n);
CSWIFT = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 8]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(xBm,zBm,Bmode,dynRange)
axis image
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(xBm,zBm,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
axis image
title('Ideal')

t3 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,BSWIFT,dynRange,attRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

t4 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,CSWIFT,dynRange,bsRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

colormap(t1,gray)
colormap(t3,turbo)
colormap(t4,parula)


% Weight map
w = (1-reject)*(abs(CSWIFT)<ratioCutOff)+reject;
wExt = movmin(w,extension);

%% SWIFT: second iteration
% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
fprintf("Inner terations:\n")
tic
[Bn,Cn,ite] = optimSwift(A1w,A2w,bw,muBswift,muCswift,m,n,tol,w);
exTime = toc;
% fprintf('\nADMM Execution time: %.4f\n',exTime)
fprintf('\nADMM Number of iterations: %d\n\n',ite)
BSWIFT = reshape(Bn*NptodB,m,n);
CSWIFT = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 8]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(xBm,zBm,Bmode,dynRange)
axis image
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(xBm,zBm,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
colormap(t2,turbo)
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
axis image
title('Ideal')

t3 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,BSWIFT,dynRange,attRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

t4 = nexttile;
[~,hB,hColor] = imOverlayInterp(Bmode,CSWIFT,dynRange,bsRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

colormap(t1,gray)
colormap(t3,turbo)
colormap(t4,parula)


%% RED no weigths
mu_aux = 10^8;
tic
[err_fp2,nfp2 ,u2]  =  admm_red_median(A'*A,A'*b(:),mu_aux,0.001,2*m*n,1500,4,1,5,m,n,mu_aux/1);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

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
imagesc(x_ACS,z_ACS,BRED, attRange)
colormap(t1,turbo)
axis image
title('RED')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRED, bsRange)
colormap(t1,"parula")
axis image
title('RED')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')


%% RED weigthed
Aw = [A1w,A2w];

mu_aux = 10^7;
tic
[err_fp2,nfp2 ,u2]  =  admm_red_median(Aw'*Aw,Aw'*bw(:),mu_aux,0.001,2*m*n,1500,4,1,5,m,n,mu_aux/1);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

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
imagesc(x_ACS,z_ACS,BRED, attRange)
colormap(t1,turbo)
axis image
title('RED')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CRED, bsRange)
colormap(t1,"parula")
axis image
title('RED')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')


