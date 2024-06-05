% ======================================================================
% ======================================================================
%% PHANTOMSSS
clear, clc
targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% resultsDir = 'C:\Users\sebas\Pictures\Journal2024\24-05-17';

% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\phantoms\ID316V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\phantoms\ID544V2\06-08-2023-Generic'];
% resultsDir = 'C:\Users\smerino.C084288\Pictures\JOURNAL\24-02-20\BS_8_12';

rawFiles = dir([targetDir,'\*.rf']);
targetFiles = dir([targetDir,'\*.mat']);
targetFiles = targetFiles(end-2:end);
resultsDir = fullfile(targetDir,'results','24-05-20');
if ~exist("resultsDir","dir"); mkdir(resultsDir); end
tableName = 'phantoms.xlsx';

%% Constants
blocksize = 8;     % Block size in wavelengths
freq_L = 2.5e6; freq_H = 7.5e6;
freq_C = 5e6;

overlap_pc      = 0.8;
ratio_zx        = 12/8;
x_inf = 0.1; x_sup = 3.8;
z_inf = 0.2; z_sup = 3.5;
NptodB = log10(exp(1))*20;

% Weights SWTV
aSNR = 1; bSNR = 0.1;
desvMin = 15;

% Weight map
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

groundTruthTargets = [0.97,0.95,0.95,0.55];

% Plotting constants
dynRange = [-50,0];
attRange = [0.4,1.1];

tol = 1e-3;

c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%% For looping each phantom

iAcq = 2;

switch iAcq
    % Optimal reg for BS 8x12, rectangular ROIs
    case 1
        muBtv = 10^3.5; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^2;
        muBwfr = 10^3.5; muCwfr = 10^2;
    case 2
        muBtv = 10^3.5; muCtv = 10^2;
        muBswtv = 10^3; muCswtv = 10^0;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^4; muCwfr = 10^1;
    case 3
        muBtv = 10^3.5; muCtv = 10^1.5;
        muBswtv = 10^2.5; muCswtv = 10^0.5;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^3.5; muCwfr = 10^1;

end

switch iAcq
    case 1
        c1x = 1.8; c1z = 1.9;
    case 2
        c1x = 1.95; c1z = 1.95;
    case 3
        c1x = 1.85; c1z = 1.9;

end

%%
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

% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% figure,plot(fpxx/1e6,meanSpectrum)
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, 0.1);
% xline([freq_L,freq_H]/1e6)
% xlim([0 15])
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
if true %iAcq == 1
    % Generating references
    att_ref = 0.53*f/8.686; % From phantom especifications
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

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
rInc = 0.95;
inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;

x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
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

AttInterp = interp2(X,Z,BR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'TV';
MetricsTV(iAcq) = r;


%% NEW WEIGHTS
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn,m,n)*NptodB;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Regularization: Au = b
tic
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
toc
BRWTik = (reshape(Bn*NptodB,m,n));
CRWTik = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BRWTik,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'WFR';
MetricsWFR(iAcq) = r;

%% Iterative WFR
bscRange = [-20,20]/2;
nIter = 5;

% Initialization
muB = 10^3.5;
muC = 10^1;
bk = b; 
Bant = zeros(m*n,1);
w = ones(m,n);

costC = zeros(1,nIter);
diffB = zeros(1,nIter);

figure('Units','centimeters', 'Position',[5 5 25 8]);
tl = tiledlayout(2,nIter+1, "Padding","tight");
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';

nexttile(2+nIter), axis off,

for ii=1:nIter
% Weighted system of eq
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);
A1w = W*A1;
A2w = W*A2;

% Solving
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB,muC,m,n,tol,mask(:),w);

% Computing weight map
bscMap = reshape(Cn,m,n)*NptodB;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

costC(ii) = sum(abs(Cn));
diffB(ii) = norm(Bn - Bant);
Bant = Bn;

Cn = reshape(Cn,m,n);
Biter = reshape(Bn,m,n)*NptodB;
t2 = nexttile(ii+1);
imagesc(x_ACS,z_ACS,Biter, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title(['ACS Iter',num2str(ii)])

t2 = nexttile(ii+2+nIter);
imagesc(x_ACS,z_ACS,Cn*NptodB, bscRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,parula)
axis image
title(['BSC Iter',num2str(ii)])

end

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(costC)
xlabel('Iteration')
ylabel('||C||_1')
grid on
xlim([1 nIter])

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(diffB)
xlabel('Iteration')
ylabel('||B_{k} - B_{k-1}||_2')
grid on
xlim([1 nIter])

save_all_figures_to_directory(resultsDir,'iterWFRfig');
close all
%% Iterative Tikhonov
bscRange = [-20,20]/2;
% Initialization
muB = 1e3;
muC = 1e1;
bk = b; 
Bant = zeros(m*n,1);

nIter = 8;
costC = zeros(1,nIter);
diffB = zeros(1,nIter);

figure('Units','centimeters', 'Position',[5 5 25 8]);
tl = tiledlayout(2,nIter+1, "Padding","tight");
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';

nexttile(2+nIter), axis off,

for ii=1:nIter
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,bk(:),muB,muC,m,n,tol,mask(:));
costC(ii) = sum(abs(Cn));
diffB(ii) = norm(Bn - Bant);
Bant = Bn;

Cn = reshape(Cn,m,n);
Biter = reshape(Bn,m,n)*NptodB;
t2 = nexttile(ii+1);
imagesc(x_ACS,z_ACS,Biter, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title(['ACS Iter',num2str(ii)])

t2 = nexttile(ii+2+nIter);
imagesc(x_ACS,z_ACS,Cn*NptodB, bscRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,parula)
axis image
title(['BSC Iter',num2str(ii)])

bk = bk - Cn;
end

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(costC)
xlabel('Iteration')
ylabel('||C||_1')
grid on
xlim([1 nIter])

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(diffB)
xlabel('Iteration')
ylabel('||B_{k} - B_{k-1}||_2')
grid on
xlim([1 nIter])

save_all_figures_to_directory(resultsDir,'iterTVL1fig');
close all
%% Iterative TV
bscRange = [-20,20]/2;
% Initialization
muB = 1e3;
muC = 1e1;
bk = b; 
Bant = zeros(m*n,1);

nIter = 5;
costC = zeros(1,nIter);
diffB = zeros(1,nIter);

figure('Units','centimeters', 'Position',[5 5 25 8]);
tl = tiledlayout(2,nIter+1, "Padding","tight");
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';

nexttile(2+nIter), axis off,

for ii=1:nIter
[Bn,Cn] = AlterOpti_ADMM(A1,A2,bk(:),muB,muC,m,n,tol,mask(:));
costC(ii) = sum(abs(Cn));
diffB(ii) = norm(Bn - Bant);
Bant = Bn;

Cn = reshape(Cn,m,n);
Biter = reshape(Bn,m,n)*NptodB;

t2 = nexttile(ii+1);
imagesc(x_ACS,z_ACS,Biter, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title(['ACS Iter',num2str(ii)])

t2 = nexttile(ii+2+nIter);
imagesc(x_ACS,z_ACS,Cn*NptodB, bscRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,parula)
axis image
title(['BSC Iter',num2str(ii)])

bk = bk - Cn;
end

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(costC)
xlabel('Iteration')
ylabel('||C||_1')
grid on

figure('Units','centimeters', 'Position',[5 5 8 8]),
plot(diffB)
xlabel('Iteration')
ylabel('||B_{k} - B_{k-1}||_2')
grid on

save_all_figures_to_directory(resultsDir,'iterTV');
close all
%%
figure('Units','centimeters', 'Position',[5 5 25 4]);
tl = tiledlayout(1,4, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(t1,gray)
title('B-mode')
%subtitle(' ')
c = colorbar('Location', 'westoutside');
c.Label.String = 'dB';
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

% fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('TV')
% fontsize(gcf,8,'points')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
% fontsize(gcf,8,'points')
% hold on 
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
% hold off

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Inc')
disp(results1.biasInc)
disp(results2.biasInc)
disp(results3.biasInc)
disp(results4.biasInc)

disp('Bias Back')
disp(results1.biasBack)
disp(results2.biasBack)
disp(results3.biasBack)
disp(results4.biasBack)

disp('RMSE Inc')
disp(results1.rmseInc)
disp(results2.rmseInc)
disp(results3.rmseInc)
disp(results4.rmseInc)

disp('RMSE Back')
disp(results1.rmseBack)
disp(results2.rmseBack)
disp(results3.rmseBack)
disp(results4.rmseBack)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)

T = [results1;results2;results3;results4];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);
%%
