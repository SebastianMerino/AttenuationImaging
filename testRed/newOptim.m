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

%% SWIFT
muBswift = 10^3.5; muCswift = 10^0.5;

% First iteration
tic
[Bn,Cn,ite] = optimSwift(A1,A2,b(:),muBswift,muCswift,...
    m,n,tol,ones(m,n));
exTime = toc;
fprintf('\nExecution time: %.4f\n',exTime)
fprintf('Number of iterations: %d\n',ite)

BSWIFT = reshape(Bn*NptodB,m,n);
CSWIFT = reshape(Cn*NptodB,m,n);

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
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(t1,turbo)
axis image
title('SWIFT ACS')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CSWIFT, bsRange)
colormap(t1,"parula")
axis image
title('SWIFT BSC')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

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
imagesc(x_ACS,z_ACS,BSWIFT, attRange)
colormap(t1,turbo)
axis image
title('SWIFT ACS')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,CSWIFT, bsRange)
colormap(t1,"parula")
axis image
title('SWIFT BSC')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

%% Iterative solvers v2
cgsIter = [70,50,25,25,3,1,62,22,15,0,65,11,0,60,5,71,1,56,1,53,0,52,0,52,0,48,0,46,0,42,0,34,0,31,0,39,0,35,0,37,0,34,0];
pcgIter = [137,95,54,15,4,1,124,40,2,112,12,0,105,4,98,1,93,0,89,0,82,0,78,0,76,0,72,0,68,0,67,0,61,0,57,0,55,0,52,0,50,0];
minresIter = [125,83,39,15,5,1,109,28,2,100,11,0,93,4,86,1,80,1,76,0,71,0,65,0,61,0,58,0,55,0,55,0,50,0,45,0,42,0,40,0,37,0];
bicgIter = [137,95,54,15,4,1,124,40,2,112,12,0,105,4,98,1,93,0,89,0,82,0,78,0,76,0,72,0,68,0,67,0,61,0,57,0,55,0,52,0,50,0];
symmlq = [136,94,53,14,3,0,123,39,1,111,11,0,104,3,97,0,92,0,88,0,81,0,77,0,75,0,71,0,67,0,66,0,60,0,56,0,54,0,51,0,49,0];
figure, histogram(cgsIter,0:10:200)
title("CGS")
xlabel("Iterations")
ylabel("Count")
xlim([0 150])

figure, histogram(pcgIter,0:10:200)
title("PCG")
xlabel("Iterations")
ylabel("Count")
xlim([0 150])

figure, histogram(minresIter,0:10:200)
title("Minres")
xlabel("Iterations")
ylabel("Count")
xlim([0 150])

figure, histogram(bicgIter,0:10:200)
title("BICG")
xlabel("Iterations")
ylabel("Count")
xlim([0 150])

figure, histogram(symmlq,0:10:200)
title("SYMMLQ")
xlabel("Iterations")
ylabel("Count")
xlim([0 150])

%% Iterative solvers v1
cgsIter = [80,81,80,80,80,80,93,97,97,101,101,101,116,116,112,112,104,104,102,102,100,100,100,100,99,99,99,99,99,99,100,100,98,98,98,98,98,98,98,98,97,97];
pcgIter = [126,130,130,130,130,130,141,141,141,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142];
minresIter = [115,119,119,119,119,119,128,128,128,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,128,128,128,128,128,128,128,128,128,128,128,128,128,128];
bicgIter = [126,130,130,130,130,130,141,141,141,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142,142];
symmlq = [125,129,129,129,129,129,140,140,140,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141,141];
figure, histogram(cgsIter)
title("CGS")
xlabel("Iterations")
ylabel("Count")
xlim([50 200])

figure, histogram(pcgIter)
title("PCG")
xlabel("Iterations")
ylabel("Count")
xlim([50 200])

figure, histogram(minresIter)
title("Minres")
xlabel("Iterations")
ylabel("Count")
xlim([50 200])

figure, histogram(bicgIter)
title("BICG")
xlabel("Iterations")
ylabel("Count")
xlim([50 200])

figure, histogram(symmlq)
title("SYMMLQ")
xlabel("Iterations")
ylabel("Count")
xlim([50 200])


%% CGS solver
iterTol1e3 = [62,72,24,7,1,44,42,1,65,1,56,1,48,13,1,49,4,3,1,49,3,3,1,51,3,1,48,3,1,42,3,1,39,3,2,1,37,3,1,39,3,1,22,5,1,18,4,1,16,4,1,16,4,1,14,4,1];
iterTol1e6 = [110,76,57,16,5,1,67,41,1,69,26,3,55,35,3,57,16,4,1,67,10,2,52,26,4,57,11,1,50,10,4,51,10,3,52,4,2,59,1,49,4,1,44,4,3,53,3,52,3,53,1,50,3];
iterTik1 = [15,10,6,4,5,4,8,7,3,7,5,7,5,7,3,7,3,7,3,7,3,7,3,6,2,5,2,5,2,5,2,5,2,5,2,5,2,5,2,5,1,5,5,5,5,5,5];
iterTik2 = [17,10,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2,16,9,7,2];

figure, histogram(iterTik2,0:2:50)
title("Tolerance 1e-3")
xlabel("CGS Iterations")
ylabel("Count")
xlim([0 50])
% figure, histogram(iterTol1e3,0:10:200)
% title("Tolerance 1e-3")
% xlabel("CGS Iterations")
% ylabel("Count")
% xlim([0 150])

% figure, histogram(iterTol1e6,0:10:200)
% title("Tolerance 1e-6")
% xlabel("CGS Iterations")
% ylabel("Count")
% xlim([0 150])

