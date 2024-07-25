% ======================================================================
% ======================================================================
%% aaaaaa
setup,
baseDir = 'D:\CarotidPlaque';
resultsDir = fullfile(baseDir,'results');
refDir = baseDir;
[~,~,~] = mkdir(resultsDir);

%% Constants
blocksize = 8;     % Block size in wavelengths
freq_L = 7e6; freq_H = 21e6;
freq_C = (freq_L + freq_H)/2;

overlap_pc      = 0.8;
ratio_zx        = 12/8;
% ratio_zx        = 1;
x_inf = -0.5; x_sup = 1.1;
z_inf = 1.5; z_sup = 2.6;
NptodB = log10(exp(1))*20;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;


% Plotting constants
dynRange = [-50,0];
attRange = [0,1.6];

tol = 1e-3;

% c1x = 1.95; c1z = 1.93;
% roiL = 1; roiD = 0.6;
% roiLz = 1.5;

% For looping each phantom
muBtv = 10^3.5; muCtv = 10^3.5;
% muBswtv = 10^3; muCswtv = 10^0;
muBswtv = 10^2.5; muCswtv = 10^-0.5;
muBtvl1 = 10^3.5; muCtvl1 = 10^1;
% muBwfr = 10^3.5; muCwfr = 10^1;
muBwfr = 10^3; muCwfr = 10^0.5;


%%

load(fullfile(baseDir,'DFG_270623','L228_Trans_6mm_bfsumed.mat'))

sam1 = bf(1:2031,:);
Bmode = db(hilbert(bf));
Bmode = Bmode - max(Bmode(:));
x = grid.x;
z = grid.z;
fs = param.fs;

dynRange = [-70,0];
figure, imagesc(x*100,z*100,Bmode, dynRange)
axis image
colormap gray

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

%% Cropping and finding sample sizes

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
roi = ind_x.*ind_z';
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Wavelength size
c0 = 1540;
wl = c0/mean(freq_C);   % Wavelength (m)

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
%%
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% meanSpectrum = db(meanSpectrum);
% meanSpectrum = meanSpectrum - max(meanSpectrum(:));
% figure,plot(fpxx/1e6,meanSpectrum)
% xline([freq_L,freq_H]/1e6)
% xlim([0 30])
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% grid on

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
alpha0_ref = 0.0014*10*NptodB;
att_ref = alpha0_ref*f/NptodB; % From phantom especifications
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Windows for spectrum
% windowing = tukeywin(nz/2,0.25);
windowing = hamming(nz/2);
windowing = windowing*ones(1,nx);

% For looping
refFiles = dir(fullfile(refDir,'L228_ref_BA_DFG_MP_*_bfsumed.mat'));
Nref = length(refFiles);
%%
% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    out = load(fullfile(refDir,refFiles(iRef).name));
    samRef = out.bf;
    samRef = samRef(ind_z,ind_x); % Cropping
    % figure,imagesc(db(hilbert(samRef)))

    %%
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

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
% rInc = 0.95;
% inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
% back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;

% x0mask = c1x - roiL/2;
% z0mask = c1z - roiLz/2;
% [back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

%% RSLD-TV

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

% AttInterp = interp2(X,Z,BR,Xq,Zq);
% r.meanInc = mean(AttInterp(inc),"omitnan");
% r.stdInc = std(AttInterp(inc),"omitnan");
% r.meanBack = mean(AttInterp(back),"omitnan");
% r.stdBack = std(AttInterp(back),"omitnan");
% r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
% r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
% r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
%     "omitnan") );
% r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
%     "omitnan") );
% r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
% r.method = 'TV';
% MetricsTV(iAcq) = r;

%% British Columbia Approach

% Weights
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
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% computation
tic
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,m,n,tol,mask(:),w);
toc
BRBC = (reshape(Bn*NptodB,m,n));
CRBC = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRBC,Xq,Zq);
% r.meanInc = mean(AttInterp(inc),"omitnan");
% r.stdInc = std(AttInterp(inc),"omitnan");
% r.meanBack = mean(AttInterp(back),"omitnan");
% r.stdBack = std(AttInterp(back),"omitnan");
% r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
% r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
% r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
%     "omitnan") );
% r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
%     "omitnan") );
% r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
% r.method = 'SWTV';
% MetricsSWTV(iAcq) = r;

%% TVL1
tic
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
toc
BRTik = (reshape(Bn*NptodB,m,n));
CRTik = (reshape(Cn,m,n));


AttInterp = interp2(X,Z,BRTik,Xq,Zq);
% r.meanInc = mean(AttInterp(inc),"omitnan");
% r.stdInc = std(AttInterp(inc),"omitnan");
% r.meanBack = mean(AttInterp(back),"omitnan");
% r.stdBack = std(AttInterp(back),"omitnan");
% r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
% r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
% r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
%     "omitnan") );
% r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
%     "omitnan") );
% r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
% r.method = 'TVL1';
% MetricsTVL1(iAcq) = r;

%% NEW WEIGHTS
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBwfr,muCwfr,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;

w = (1-reject)* (1./((bscMap/ratioCutOff).^(2*order) + 1)) + reject;
% %%
% [X,Z] = meshgrid(x_ACS, z_ACS);
% inside = X>0.15 & X<0.44 & Z>2.156 & Z<2.3445;
% figure, imagesc(x_ACS,z_ACS,inside)
% 
% figure, imagesc(x,z,Bmode),
% colorbar
% 
% w(inside) = 1;
% figure, imagesc(x_ACS,z_ACS,w),
% colorbar

wExt = movmin(w,extension);


%%
% New equations
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BRWFR = (reshape(Bn*NptodB,m,n));

AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
% r.meanInc = mean(AttInterp(inc),"omitnan");
% r.stdInc = std(AttInterp(inc),"omitnan");
% r.meanBack = mean(AttInterp(back),"omitnan");
% r.stdBack = std(AttInterp(back),"omitnan");
% r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
% r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
% r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
%     "omitnan") );
% r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
%     "omitnan") );
% r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
% r.method = 'SWIFT';
% MetricsWFR(iAcq) = r;

%%

figure('Units','centimeters', 'Position',[5 2 25 8]);
tl = tiledlayout(1,4, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar;
c.Label.String = 'dB';
% hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('RSLD')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
% hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t3 = nexttile;
imagesc(x_ACS,z_ACS,BRBC, attRange)
% xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title('SWTV-ACE')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
% hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% t4 = nexttile;
% imagesc(x_ACS,z_ACS,BRTik, attRange)
% xlabel('Lateral [cm]'), ylabel('Axial [cm]')
% colormap(t4,turbo)
% axis image
% title('TVL1')
% % fontsize(gcf,8,'points')
% hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWFR, attRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('SWIFT')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';

% fontsize(gcf,8,'points')
% hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off


%%
save_all_figures_to_directory(resultsDir,'fig');
close all
% % %%
% results1 = struct2table(MetricsTV);
% results2 = struct2table(MetricsSWTV);
% results3 = struct2table(MetricsTVL1);
% results4 = struct2table(MetricsWFR);
%
% disp('Bias Inc')
% disp(results1.biasInc)
% disp(results2.biasInc)
% disp(results3.biasInc)
% disp(results4.biasInc)
%
% disp('Bias Back')
% disp(results1.biasBack)
% disp(results2.biasBack)
% disp(results3.biasBack)
% disp(results4.biasBack)
%
% disp('RMSE Inc')
% disp(results1.rmseInc)
% disp(results2.rmseInc)
% disp(results3.rmseInc)
% disp(results4.rmseInc)
%
% disp('RMSE Back')
% disp(results1.rmseBack)
% disp(results2.rmseBack)
% disp(results3.rmseBack)
% disp(results4.rmseBack)
%
% disp('CNR')
% disp(results1.cnr)
% disp(results2.cnr)
% disp(results3.cnr)
% disp(results4.cnr)
%
% T = [results1;results2;results3;results4];
% writetable(T,fullfile(resultsDir,tableName),...
%      'WriteRowNames',true);
%%
