% ====================================================================== %
% Script for finding optimal reg for data with an inclusion 
% Created on Jan 31, 2024
% ====================================================================== %
clc, clear,

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\24_04_04_layered'];
refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\24_04_25_ref'];

resultsDir = [targetDir,'\results\opt-reg'];
mkdir(resultsDir);

targetFiles = dir([targetDir,'\rf*.mat']);
% targetFiles = targetFiles(2:3);
refFiles = dir([refDir,'\rf*.mat']);

%%
% SETTING PARAMETERS
blocksize = 8;     % Block size in wavelengths
freq_L = 3.5e6; freq_H = 8.5e6; % freq_L = 3.3e6; freq_H = 8.7e6;

overlap_pc      = 0.8;
ratio_zx        = 12/8;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.2;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% Plotting
dynRange = [-40,0];
% attRange = [0.4,1.4];
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% GT
groundTruthTop = [0.5,0.5,0.5];
groundTruthBottom = [1,1,1];

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;

%% Setting up

for iAcq = 1:length(targetFiles)
load(fullfile(targetDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);
dynRange = [-50,0];
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));

%% Cropping and finding sample sizes

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);
Bmode = Bmode(ind_z,ind_x);

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
% nz = 2*round(blocksize*wl/dz /2);
L = (nz/2)*dz*100;   % (cm)
z0p = 1:wz:length(z)-nz;
z0d = z0p + nz/2;
z_ACS = z(z0p+ nz/2);
m  = length(z0p);


%% Spectrum
% BW from spectrogram
[pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
meanSpectrum = mean(pxx,2);
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
figure,plot(fpxx/1e6,db(meanSpectrum./max(meanSpectrum)))
xline([freq_L,freq_H]/1e6)
xlabel('Frequency [MHz]')
ylabel('Magnitude')
xlim([0 15])
grid on

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
fprintf('Diff x: %.2f mm, z: %.2f mm\n',wx*dx*1e3,wz*dz*1e3)
fprintf('Region of interest rows: %i, col: %i\n\n',m,n);

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

for iRef = 1:Nref
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

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
% [top,bottom] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

attIdeal = ones(size(Zq))*groundTruthTop(iAcq);
attIdeal(Zq>2) = groundTruthBottom(iAcq);
top = Zq < 1.9; % 0.1 cm interface
bottom = Zq > 2.1;

figure('Units','centimeters', 'Position',[5 5 20 6]);
tiledlayout(1,2)
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)])
ylim([z_ACS(1) z_ACS(end)])
colormap(t1,gray)
colorbar('westoutside')
title('Bmode')

t2 = nexttile;
imagesc(x,z,attIdeal,attRange)
axis image
colormap(t2,turbo);
c = colorbar('westoutside');
c.Label.String = 'dB/cm/MHz';
title('Ideal ACS')

% %% Spectrum
% % Heterogeneous
% region1 = (X.^2 + (Z-2).^2) <= rBottom^2;
% region2 = (X.^2 + (Z-2).^2) >= rBottom^2;
% sld1 = squeeze(sum(sum(b.*region1,1),2))/sum(region1(:)) * NptodB /4/L;
% acs1 = f\sld1;
% fprintf('Attenuation is %.2f\n',acs1)
% sld2 = squeeze(sum(sum(b.*region2,1),2))/sum(region2(:)) * NptodB /4/L;
% acs2 = f\sld2;
% fprintf('Attenuation is %.2f\n',acs2)
% figure, plot(f,sld1)
% hold on
% plot(f,sld2)
% plot(f,acs1*f, 'k--')
% plot(f,acs2*f, 'k--')
% hold off
% grid on,
% xlim([0,max(f)]), ylim([0 15]),
% xlabel('Frequency [MHz]')
% ylabel('Att. [dB/cm]')
% title('Mean SLD')
% legend('bottom','top')


%% RSLD
muB = 10.^(1:0.5:3.5);
% muC = 10.^(1:0.5:3);
minRMSE = 100;
for mmB = 1:length(muB)
    % for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muB(mmB),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,...
            "omitnan") ;
        RmseTop = mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseBottom + RmseTop)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muB(mmB);
            BRopt = BR;
            CRopt = CR;
        end
    % end
end

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'Isotropic RSLD')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu_b=10^{',num2str(log10(muBopt),2),'}'])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu_c=10^{',num2str(log10(muCopt),2),'}'])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
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
        
        temp = [sub_block_p(:);sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


muB = 10.^(2:0.5:3);
muC = 10.^(0:0.5:3);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
        m,n,tol,mask(:),w);
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);
        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,...
            "omitnan") ;
        RmseTop = mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseBottom + RmseTop)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

figure('Units','centimeters', 'Position',[5 5 22 6]);
tl = tiledlayout(1,3, "Padding","tight");
title(tl,'RSLD - SWTV by British Columbia')
t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Weights')
c = colorbar;

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu_b=10^{',num2str(log10(muBopt),2),'}'])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu_c=10^{',num2str(log10(muCopt),2),'}'])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%%
% %% Minimizing BS log ratio
% muB = 10.^(2.5:0.5:4);
% muC = 10.^(0:0.5:3);
% minRMSE = 100;
% for mmB = 1:length(muB)
%     for mmC = 1:length(muC)
%         tic
%         [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
%         toc
%         BR = reshape(Bn*NptodB,m,n);
%         CR = reshape(Cn*NptodB,m,n);
% 
%         % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
%         AttInterp = interp2(X,Z,BR,Xq,Zq);
% 
%         RmseBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,...
%             "omitnan") ;
%         RmseTop = mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,...
%             "omitnan");
%         RMSE = sqrt((RmseBottom + RmseTop)/2);
% 
%         if RMSE<minRMSE
%             minRMSE = RMSE;
%             muBopt = muB(mmB);
%             muCopt = muC(mmC);
%             BRopt = BR;
%             CRopt = CR;
%         end
%     end
% end

% figure('Units','centimeters', 'Position',[5 5 15 6]);
% tl = tiledlayout(1,2, "Padding","tight");
% title(tl,'RSLD with TV(B)+||C||_1')
% 
% t2 = nexttile; 
% imagesc(x_ACS,z_ACS,BRopt, attRange)
% colormap(t2,turbo)
% axis image
% title(['RSLD-TVL1, \mu_b=10^{',num2str(log10(muBopt),2),'}'])
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';
% 
% t3 = nexttile; 
% imagesc(x_ACS,z_ACS,CRopt, bsRange)
% colormap(t3,parula)
% axis image
% title(['RSLD-TVL1, \mu_c=10^{',num2str(log10(muCopt),2),'}'])
% c = colorbar;
% c.Label.String = 'BS log ratio [dB]';
% 
% AttInterp = interp2(X,Z,BRopt,Xq,Zq);
% r.meanTop = mean(AttInterp(top),"omitnan");
% r.stdTop = std(AttInterp(top),"omitnan");
% r.meanBottom = mean(AttInterp(bottom),"omitnan");
% r.stdBottom = std(AttInterp(bottom),"omitnan");
% r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
% r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
% r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
% MetricsTVL1(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS
muB = 10.^(2.5:0.5:4);
muC = 10.^(-0.5:0.5:2);
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        % First iteration
        [~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        bscMap = reshape(Cn*NptodB,m,n);
        
        % Weight map
        w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
        w = movmin(w,extension);

        % Weight matrices and new system
        W = repmat(w,[1 1 p]);
        W = spdiags(W(:),0,m*n*p,m*n*p);
        bw = W*b(:);        
        A1w = W*A1;
        A2w = W*A2;

        % Second iteration
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        BR = reshape(Bn*NptodB,m,n);
        CR = (reshape(Cn*NptodB,m,n));

        % Interp and RMSE
        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseBottom = mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,...
            "omitnan") ;
        RmseTop = mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseBottom + RmseTop)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end
%%
figure('Units','centimeters', 'Position',[5 5 22 6]);
tl = tiledlayout(1,3, "Padding","tight");
title(tl,'Weighted Fidelity and Regularization') 

t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Weights')
c = colorbar;
%c.Label.String = 'BS log ratio [dB]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD-SWIFT, \mu_b=10^{',num2str(log10(muBopt),2),'}'])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD-SWIFT, \mu_c=10^{',num2str(log10(muCopt),2),'}'])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");

r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWIFT(iAcq) = r;

%%
save_all_figures_to_directory(resultsDir,['sim',num2str(iAcq),'Figure']);
close all

end

% results1 = struct2table(MetricsTV);
% results2 = struct2table(MetricsSWTV);
% results3 = struct2table(MetricsTVL1);
% results4 = struct2table(MetricsSWIFT);
% 
% disp('Bias Top')
% disp(results1.biasTop)
% disp(results2.biasTop)
% disp(results3.biasTop)
% disp(results4.biasTop)
% 
% disp('Bias Bottom')
% disp(results1.biasBottom)
% disp(results2.biasBottom)
% disp(results3.biasBottom)
% disp(results4.biasBottom)
% 
% disp('CNR')
% disp(results1.cnr)
% disp(results2.cnr)
% disp(results3.cnr)
% disp(results4.cnr)