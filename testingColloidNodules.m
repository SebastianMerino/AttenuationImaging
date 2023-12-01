% Script to explore colloid nodules
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');

%% Clinical case
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\ThyroidSelected\CUELLO#3'];
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\ThyroidSelected\CUELLO#2'];
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\ThyroidSelected\coloide'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\ThyroidSelected\carcinoma'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];
figDir = [baseDir,'\fig\01-12'];

targetFiles = dir([targetDir,'\*.mat']);

blocksize = 10;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8e6;
overlap_pc      = 0.8;
ratio_zx        = 1;

%% Loading case FULL VERSION
iAcq = 5;
load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

sam1 = RF(:,:,1);

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));

%% Cropping and finding sample sizes

% Manual cropping
dynRange = [-50,-10];
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
confirmation = '';
while ~strcmp(confirmation,'Yes')
    rect = getrect;
    confirmation = questdlg('Sure?');
end
close,
x_inf = rect(1); x_sup = rect(1)+rect(3);
z_inf = rect(2); z_sup = rect(2)+rect(4);

% Actual crop
% switch iAcq
%     case 1
%         limitsX = [0.31    2.16];
%         limitsZ = [0.38    2.26];
%         % rect = [1.3611    0.7675    1.8377    0.8917]
%     case 3
%         limitsX = [0.64    2.59];
%         limitsZ = [0.87    2.92];
%         % rect = [2.1598    1.3335    0.7676    1.1631]
%     case 9
%         limitsX = [0.98    3.62];
%         limitsZ = [1.42    3.01];
% end
% x_inf = limitsX(1); x_sup = limitsX(2);
% z_inf = limitsZ(1); z_sup = limitsZ(2);

% Limits for ACS estimation
ind_x = x_inf <= xFull & xFull <= x_sup;
ind_z = z_inf <= zFull & zFull <= z_sup;
roi = ind_x.*ind_z';
x = xFull(ind_x);
z = zFull(ind_z);
sam1 = sam1(ind_z,ind_x);

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
nz = 2*round(blocksize*wl/dz /2); % Window size
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);

% Frequency samples
NFFT = 2^(nextpow2(nz/2)+2);
band = (0:NFFT-1)'/NFFT * fs;   % [Hz] Band of frequencies
rang = band > freq_L & band < freq_H ;   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(f);

% Plot region of interest B-mode image
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation

% Generating references
att_ref = attenuation_phantoms_Np(f, 3, []);
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
windowing = tukeywin(nz/2,0.25);
windowing = windowing*ones(1,nx);

% For looping
refFiles = dir([refDir,'\*.mat']);
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load([refDir,'\',refFiles(iRef).name]);
    samRef = out.RF;
    samRef = samRef(ind_z,ind_x); % Cropping
    % figure,imagesc(db(hilbert(samRef)))
    for jj=1:n
        for ii=1:m
            xw = x0(jj) ;   % x window
            zp = z0p(ii);
            zd = z0d(ii);

            sub_block_p = samRef(zp:zp+nz/2-1,xw:xw+nx-1);
            sub_block_d = samRef(zd:zd+nz/2-1,xw:xw+nx-1);
            [tempSp,~] = spectra(sub_block_p,windowing,0,nz/2,NFFT);
            [tempSd,~] = spectra(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
            Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
        end
    end
end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref

%% Setting up
% Spectrum
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

% System of eq
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
b = (log(Sp) - log(Sd)) - (compensation);

% Optimization constants
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Plotting constants
dynRange = [-40,-5];
attRange = [0.3,1.7];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% RSLD

% muB = 10.^(2.5:0.5:3.5);
% muC = 10.^(0:1:2);
% 
% for mmB = 1:length(muB)
%     for mmC = 1:length(muC)
%         tic
%         [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
%         toc
%         BR = reshape(Bn*NptodB,m,n);
%         CR = reshape(Cn*NptodB,m,n);
%         figure('Units','centimeters', 'Position',[5 5 30 8]);
%         tl = tiledlayout(1,3);
%         title(tl,'RSLD-TV')
%         subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
%         t1 = nexttile;
%         imagesc(x,z,Bmode,dynRange)
%         axis equal
%         xlim([x_ACS(1) x_ACS(end)]),
%         ylim([z_ACS(1) z_ACS(end)]),
%         colormap(t1,gray)
%         colorbar(t1,'westoutside')
%         title('Bmode')
% 
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,BR, attRange)
%         colormap(t2,turbo)
%         axis image
%         title(['RSLD, \mu=',num2str(muB(mmB),2)])
%         c = colorbar;
%         c.Label.String = 'Att. [db/cm/MHz]';
% 
%         t3 = nexttile; 
%         imagesc(x_ACS,z_ACS,CR, bsRange)
%         colormap(t3,parula)
%         axis image
%         title(['RSLD, \mu=',num2str(muC(mmC),2)])
%         c = colorbar;
%         c.Label.String = 'BS log ratio (a.u.)';
%     end
% end

%% British Columbia WEIGHTS
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
aSNR = 1; bSNR = 0.1;
desvMin = 15;
wBC = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

figure('Units','centimeters', 'Position',[5 5 30 8]),
tl = tiledlayout(1,3);
title(tl,{'New Weights',''});
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,SNR)
colormap(t2,parula)
c = colorbar;
% ylabel(c,'dB')
axis image
title('SNR')

t3 = nexttile;
imagesc(x_ACS,z_ACS,wBC)
colormap(t3,parula)
colorbar;
axis image
title('Weights')
%% British columbia approach
% muB = 10.^(2.5:0.5:3.5);
% muC = 10.^(-0.5:0.5:0.5);
% 
% for mmB = 1:length(muB)
%     for mmC = 1:length(muC)
%         tic
%         [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
%         m,n,tol,mask(:),wBC);
%         toc
%         BR = reshape(Bn*NptodB,m,n);
%         CR = reshape(Cn*NptodB,m,n);
%         figure('Units','centimeters', 'Position',[5 5 30 8]);
%         tl = tiledlayout(1,3);
%         title(tl,'British Columbia Approach')
%         subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
%         t1 = nexttile;
%         imagesc(x,z,Bmode,dynRange)
%         axis equal
%         xlim([x_ACS(1) x_ACS(end)]),
%         ylim([z_ACS(1) z_ACS(end)]),
%         colormap(t1,gray)
%         colorbar(t1,'westoutside')
%         title('Bmode')
% 
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,BR, attRange)
%         colormap(t2,turbo)
%         axis image
%         title(['RSLD, \mu=',num2str(muB(mmB),2)])
%         c = colorbar;
%         c.Label.String = 'Att. [db/cm/MHz]';
% 
%         t3 = nexttile; 
%         imagesc(x_ACS,z_ACS,CR, bsRange)
%         colormap(t3,parula)
%         axis image
%         title(['RSLD, \mu=',num2str(muC(mmC),2)])
%         c = colorbar;
%         c.Label.String = 'BS log ratio (a.u.)';
%     end
% end


%% Minimizing BS log ratio

% % muB = 10.^(2:0.5:3);
% muC = 10.^(0.5:0.5:1.5);
% muB = 10.^(2:0.5:3);
% % muC = 10.^(1:0.5:2);
% for mmB = 1:length(muB)
%     for mmC = 1:length(muC)
%         tic
%         [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
%         toc
%         BR = reshape(Bn*NptodB,m,n);
%         CR = reshape(Cn*NptodB,m,n);
%         figure('Units','centimeters', 'Position',[5 5 30 8]);
%         tl = tiledlayout(1,3);
%         title(tl,'TV-L1')
%         subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
%         t1 = nexttile;
%         imagesc(x,z,Bmode,dynRange)
%         axis equal
%         xlim([x_ACS(1) x_ACS(end)]),
%         ylim([z_ACS(1) z_ACS(end)]),
%         colormap(t1,gray)
%         colorbar(t1,'westoutside')
%         title('Bmode')
% 
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,BR, attRange)
%         colormap(t2,turbo)
%         axis image
%         title(['RSLD, \mu=',num2str(muB(mmB),2)])
%         c = colorbar;
%         c.Label.String = 'Att. [db/cm/MHz]';
% 
%         t3 = nexttile; 
%         imagesc(x_ACS,z_ACS,CR, bsRange)
%         colormap(t3,parula)
%         axis image
%         title(['RSLD, \mu=',num2str(muC(mmC),2)])
%         c = colorbar;
%         c.Label.String = 'BS log ratio (a.u.)';
%     end
% end


%% NEW WEIGHTS
% First estimation
muB0 = 1e3; muC0 = 10^0;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB0,muC0,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

% Weight function
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% New equations
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

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
imagesc(x_ACS,z_ACS,bscMap, [-20 20])
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
% muB = 10.^(2:0.5:3)*10;
% muC = 10.^(-0.5:0.5:0.5)*10;
% 
% for mmB = 1:length(muB)
%     for mmC = 1:length(muC)
%         tic
%         [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
%         toc
%         BR = reshape(Bn*NptodB,m,n);
%         CR = reshape(Cn*NptodB,m,n);
%         figure('Units','centimeters', 'Position',[5 5 30 8]);
%         tl = tiledlayout(1,3);
%         title(tl,'TV-L1')
%         subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})
%         t1 = nexttile;
%         imagesc(x,z,Bmode,dynRange)
%         axis equal
%         xlim([x_ACS(1) x_ACS(end)]),
%         ylim([z_ACS(1) z_ACS(end)]),
%         colormap(t1,gray)
%         colorbar(t1,'westoutside')
%         title('Bmode')
% 
%         t2 = nexttile; 
%         imagesc(x_ACS,z_ACS,BR, attRange)
%         colormap(t2,turbo)
%         axis image
%         title(['RSLD, \mu=',num2str(muB(mmB),2)])
%         c = colorbar;
%         c.Label.String = 'Att. [db/cm/MHz]';
% 
%         t3 = nexttile; 
%         imagesc(x_ACS,z_ACS,CR, bsRange)
%         colormap(t3,parula)
%         axis image
%         title(['RSLD, \mu=',num2str(muC(mmC),2)])
%         c = colorbar;
%         c.Label.String = 'BS log ratio (a.u.)';
%     end
% end

% ======================================================================
% ======================================================================
% ======================================================================
%% ACUMULADO
% Zona heterogenea
% switch iAcq
%     case 1
%         muBtv = 10^3; muCtv = 10^1;
%         muBswtv = 10^2.5; muCswtv = 10^1;
%         muBtvl1 = 10^3; muCtvl1 = 10^0.5;
%         muBwfr = 10^3.5; muCwfr = 10^0.5;
%     case 3
%         muBtv = 10^3; muCtv = 10^1;
%         muBswtv = 10^2.5; muCswtv = 10^1;
%         muBtvl1 = 10^3; muCtvl1 = 10^-0.5;
%         muBwfr = 10^3.3; muCwfr = 10^-0.5;
%     case 9
%         muBtv = 10^3; muCtv = 10^1;
%         muBswtv = 10^2.5; muCswtv = 10^1;
%         muBtvl1 = 10^2.5; muCtvl1 = 10^-0.5;
%         muBwfr = 10^3; muCwfr = 10^-0.5;
% end

% Inclusion
% switch iAcq
%     case 1
%         muBtv = 10^4; muCtv = 10^2;
%         muBswtv = 10^3.5; muCswtv = 10^2;
%         muBtvl1 = 10^4; muCtvl1 = 10^1.5;
%         muBwfr = 10^4.5; muCwfr = 10^1.5;
%     case 3
        % muBtv = 10^4.5; muCtv = 10^2;
        % muBswtv = 10^4; muCswtv = 10^2;
        % muBtvl1 = 10^4.5; muCtvl1 = 10^0.5;
        % muBwfr = 10^4; muCwfr = 10^1;
%     case 9
%         muBtv = 10^4; muCtv = 10^1;
%         muBswtv = 10^3.5; muCswtv = 10^1;
%         muBtvl1 = 10^3.5; muCtvl1 = 10^-0.5;
%         muBwfr = 10^4; muCwfr = 10^-0.5;
% end

% % Casos coloides cuello 3
% switch iAcq
%     case 3
%         muBtv = 10^2.5; muCtv = 10^0;
%         muBswtv = 10^2; muCswtv = 10^0;
%         muBtvl1 = 10^2.5; muCtvl1 = 10^0.5;
%         muBwfr = 10^2.5; muCwfr = 10^0.5;
%     case 4
%         muBtv = 10^3.5; muCtv = 10^0;
%         muBswtv = 10^3; muCswtv = 10^0;
%         muBtvl1 = 10^3.5; muCtvl1 = 10^0.5;
%         muBwfr = 10^4; muCwfr = 10^0.5;
%         % muBtv = 10^2.5; muCtv = 10^1;
%         % muBswtv = 10^2; muCswtv = 10^1;
%         % muBtvl1 = 10^2.5; muCtvl1 = 10^0.5;
%         % muBwfr = 10^3; muCwfr = 10^1;
% end

muBtv = 10^3.5; muCtv = 10^0;
muBswtv = 10^3; muCswtv = 10^0;
muBtvl1 = 10^3.5; muCtvl1 = 10^0.5;
muBwfr = 10^4; muCwfr = 10^0.5;

% RSLD-TV
[Bn,~] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);

% British Columbia
[Bn,~] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
m,n,tol,mask(:),wBC);
BRSWTV = reshape(Bn*NptodB,m,n);

% TV + L1 (no weights)
[Bn,~] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);

% WFR   
[Bn,~] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);

% Plots

figure('Units','centimeters', 'Position',[5 5 24 12]);
tl = tiledlayout(2,2);
title(tl,{'Comparison'})
subtitle(tl,{['Patient ',targetFiles(iAcq).name(1:end-4)],''})

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t2,turbo)
axis equal tight
title('TV')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t2,turbo)
axis equal tight
title('SWTV')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t2,turbo)
axis equal tight
title('TV-L1')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t2,turbo)
axis equal tight
title('WFR')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%% mean
fileID = fopen('results.txt','a');
fprintf(fileID,'Patient %s\n',targetFiles(iAcq).name(1:end-4));
fprintf(fileID,'TV   ACS = %.2f\n',mean(BRTV,'all'));
fprintf(fileID,'SWTV ACS = %.2f\n',mean(BRSWTV,'all'));
fprintf(fileID,'TVL1 ACS = %.2f\n',mean(BRTVL1,'all'));
fprintf(fileID,'WFR  ACS = %.2f\n',mean(BRWFR,'all'));
fclose(fileID);

%%
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure('Units','centimeters', 'Position',[5 5 10 15])
tiledlayout(2,1)
t2 = nexttile;
imagesc(xFull,zFull,BmodeFull,[-50 0])
axis image
title('B-mode')
ylim([0.1, 3.5])
c = colorbar;
c.Label.String = 'dB';
hold on
contour(xFull,zFull,roi,1,'w--')
hold off

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BRWFR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.1, 3.5])
hold on
contour(xFull,zFull,roi,1,'w--')
hold off

colormap(t2,gray)


%%
save_all_figures_to_directory(figDir,[targetFiles(iAcq).name(1:end-4),'fig']);
close all