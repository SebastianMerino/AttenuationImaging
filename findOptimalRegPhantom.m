% ====================================================================== %
% Script for finding optimal reg for phantom data 
% Created on Feb 7, 2024
% ====================================================================== %
clc, clear,

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
rawFiles = dir([targetDir,'\*.rf']);
refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];

resultsDir = fullfile(targetDir,'results','24-02-13','BS_8_10-BW_cut15');
if ~exist(resultsDir,"dir"); mkdir(resultsDir); end

targetFiles = dir([targetDir,'\*.mat']);
targetFiles = targetFiles(1:end);

blocksize = 8;     % Block size in wavelengths
freq_L = 3e6; freq_H = 6e6; 
% freq_C = mean([freq_L freq_H]);
overlap_pc      = 0.8;
ratio_zx        = 12/8;

% G.T.
% groundTruthInc = [0.97,0.95,0.95];
% groundTruthBack = [0.55,0.55,0.55];
groundTruthInc = [0.52,0.55,0.74,0.61,0.75,0.97,0.95,0.95];
groundTruthBack = [0.55,0.55,0.55,0.55,0.55,0.55,0.55,0.55];

% Weights SWTV
aSNR = 1; bSNR = 0.1;
desvMin = 15;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% Plotting
dynRange = [-40,0];
bsRange = [-15 15];
attRange = [0.4,1.1];
NptodB = log10(exp(1))*20;

x_inf = 0.1; x_sup = 3.8;
z_inf = 0.2; z_sup = 3.5;


%% Setting up

for iAcq = 1:length(targetFiles)
fprintf("Phantom no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(targetDir,targetFiles(iAcq).name));

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = RF(:,:,1);

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
wl = c0/freq_C;   % Wavelength (m)

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
if iAcq == 1
    % Generating references
    att_ref = 0.53*f/8.686; % From phantom especifications
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
end

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

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

% figure,
% imagesc(x,z,Bmode,dynRange)
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
% xlabel('Lateral [cm]'), ylabel('Axial [cm]')
% axis image
% colormap(gray)
% title('B-mode')
% %subtitle(' ')
% c = colorbar('Location', 'westoutside');
% c.Label.String = 'dB';
% fontsize(gcf,8,'points')
% hold on 
% rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
%     'LineStyle','--', 'LineWidth',1)
% hold off


%% RSLD
muB = 10.^(3:0.5:3.5);
muC = 10.^(1:3);
% muB = 10^ 3.5;
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'Isotropic RSLD')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
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
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


muB = 10.^(2.5:0.5:3);
muC = 10.^(0:2);

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
        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);

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
title(['RSLD, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% Minimizing BS log ratio
muB = 10.^(3:0.5:4);
muC = 10.^(0:1:2);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        % RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
        AttInterp = interp2(X,Z,BR,Xq,Zq);

        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);

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
figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'RSLD with TV(B)+||C||_1')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRopt, attRange)
colormap(t2,turbo)
axis image
title(['RSLD-TVL1, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD-TVL1, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% Minimizing BS log ratio and WEIGHTS
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),10^3,1,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);


W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;


muB = 10.^(2.5:0.5:4);
muC = 10.^(-0.5:0.5:2);
% muB = 10.^(4);
% muC = 10.^(1.5);

minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        %tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        %toc

        BR = reshape(Bn*NptodB,m,n);
        CR = (reshape(Cn*NptodB,m,n));
%        RMSE = sqrt(mean((BR-attIdeal).^2,'all'));
        AttInterp = interp2(X,Z,BR,Xq,Zq);

        RmseInc = mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,...
            "omitnan");
        % disp('RMSE INC BACK')
        % disp(sqrt(RmseInc))
        % disp(sqrt(RmseBack))
        RMSE = sqrt((RmseInc + RmseBack)/2);

        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end

        % figure('Units','centimeters', 'Position',[5 5 22 6]);
        % tl = tiledlayout(1,3, "Padding","tight");
        % title(tl,'Weighted Fidelity and Regularization')
        % 
        % t1 = nexttile; 
        % imagesc(x_ACS,z_ACS,w, [0 1])
        % colormap(t1,parula)
        % axis image
        % title('Weights')
        % c = colorbar;
        % c.Label.String = 'BS log ratio [dB]';
        % 
        % t2 = nexttile; 
        % imagesc(x_ACS,z_ACS,BRopt, attRange)
        % colormap(t2,turbo)
        % axis image
        % title(['RSLD-WFR, \mu=',num2str(muBopt,2)])
        % c = colorbar;
        % c.Label.String = 'Att. [db/cm/MHz]';
        % 
        % t3 = nexttile; 
        % imagesc(x_ACS,z_ACS,CRopt, bsRange)
        % colormap(t3,parula)
        % axis image
        % title(['RSLD-WFR, \mu=',num2str(muCopt,2)])
        % c = colorbar;
        % c.Label.String = 'BS log ratio [dB]';
        % pause(0.1)
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
title(['RSLD-WFR, \mu=',num2str(muBopt,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRopt, bsRange)
colormap(t3,parula)
axis image
title(['RSLD-WFR, \mu=',num2str(muCopt,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


AttInterp = interp2(X,Z,BRopt,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdBottom = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");

r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;

%%
save_all_figures_to_directory(resultsDir,['sim',num2str(iAcq),'Figure']);
close all

end

results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasBack)
disp(results2.biasBack)
disp(results3.biasBack)
disp(results4.biasBack)

disp('Bias Bottom')
disp(results1.biasInc)
disp(results2.biasInc)
disp(results3.biasInc)
disp(results4.biasInc)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)