% ========================================================================
% ========================================================================
% ========================================================================
%% SIMULATION ON INCLUSIONS

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\Simulation_23_12_18'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
croppedDir = [baseDir,'\cropped'];

if (~exist(croppedDir,"dir")), mkdir(croppedDir); end
targetFiles = dir([targetDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);

for iAcq = 1:length(croppedFiles)
    fprintf("Acquisition no. %i, patient %s\n",iAcq,croppedFiles(iAcq).name);
end 

blocksize = 10;     % Block size in wavelengths
freq_L = 3.3e6; freq_H = 8.7e6; % original 3.3-8.7s
overlap_pc      = 0.8;
ratio_zx        = 1;
referenceAtt    = 0.7;

figDir = 'C:\Users\sebas\Pictures\Journal2024';

groundTruthTop = [0.6,0.6,0.6];
groundTruthBottom = [1.2,1.2,1.2];

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
attRange = [0.4,1.4];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%% Setting up

for iAcq = 1:3
% Regularization
switch iAcq
    case 1
        muBtv = 10^3.5; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3; muCtvl1 = 10^1.5;
        muBwfr = 10^3.5; muCwfr = 10^2;
    case 2
        muBtv = 10^3.5; muCtv = 10^2;
        muBswtv = 10^3.5; muCswtv = 10^2.5;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^3.5; muCwfr = 10^1.5;
    case 3
        muBtv = 10^3.5; muCtv = 10^2;
        muBswtv = 10^3.5; muCswtv = 10^2.5;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^3.5; muCwfr = 10^1.5;
end

load(fullfile(targetDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = rf(:,:,1);
dynRange = [-50,0];

%% Cropping and finding sample sizes
% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.5; z_sup = 3.5;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
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

%% Spectrum
% BW from spectrogram
[pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
meanSpectrum = mean(pxx,2);
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
figure,plot(fpxx/1e6,meanSpectrum)
xline([freq_L,freq_H]/1e6)
xlabel('Frequency [MHz]')
ylabel('Magnitude')
xlim([0 15])

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
% Generating references
att_ref = referenceAtt*(f.^1.05)/(20*log10(exp(1)));
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
Nref = length(refFiles);

% Memory allocation
Sp_ref = zeros(m,n,p,Nref);
Sd_ref = zeros(m,n,p,Nref);
for iRef = 1:Nref
    load(fullfile(refDir,refFiles(iRef).name),"rf","medium");
    % disp(mean(medium.alpha_coeff(:)))
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

%% Setting up
% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% Creating reference
[X,Z] = meshgrid(x_ACS,z_ACS);
attIdeal = ones(size(Z));
rInc = 0.8;
inclusion = (X.^2 + (Z-2).^2)<= rInc^2;
attIdeal(~inclusion) = groundTruthTop(iAcq);
attIdeal(inclusion) = groundTruthBottom(iAcq); %incl = bottom

%% TV
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

axialTV{iAcq} = mean(BRTV(:,20:27),2);
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
inclusion = (Xq.^2 + (Zq-2).^2)<= rInc^2;
top = ~inclusion;
bottom = inclusion;
AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTV(iAcq) = r;
%% SWTV
% Calculating SNR
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

% Calculating weights
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
wSNR = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));

% Method
[Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muBswtv,muCswtv,...
m,n,tol,mask(:),wSNR);
BRSWTV = reshape(Bn*NptodB,m,n);
CRSWTV = reshape(Cn*NptodB,m,n);

axialSWTV{iAcq} = mean(BRSWTV(:,20:27),2);
AttInterp = interp2(X,Z,BRSWTV,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsSWTV(iAcq) = r;

%% TVL1
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

axialTVL1{iAcq} = mean(BRTVL1(:,20:27),2);
AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsTVL1(iAcq) = r;

%% WFR
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB,muC,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Computing weights
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% Setting up new system
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Method
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);
CRWFR = reshape(Cn*NptodB,m,n);

axialWFR{iAcq} = mean(BRWFR(:,20:27),2);
AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanTop = mean(AttInterp(top),"omitnan");
r.stdTop = std(AttInterp(top),"omitnan");
r.meanBottom = mean(AttInterp(bottom),"omitnan");
r.stdBottom = std(AttInterp(bottom),"omitnan");
r.biasTop = mean( AttInterp(top) - groundTruthTop(iAcq),"omitnan");
r.biasBottom = mean( AttInterp(bottom) - groundTruthBottom(iAcq),"omitnan");
r.rmseTop = sqrt(mean( (AttInterp(top) - groundTruthTop(iAcq)).^2,"omitnan"));
r.rmseBottom = sqrt(mean( (AttInterp(bottom) - groundTruthBottom(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanBottom - r.meanTop)/sqrt(r.stdTop^2 + r.stdBottom^2);
MetricsWFR(iAcq) = r;
%% Plotting
figure('Units','centimeters', 'Position',[5 5 25 4]);
tl = tiledlayout(1,5, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1, 'westoutside')
title('Bmode')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTV, attRange)
colormap(t1,turbo)
axis image
title('TV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
%c = colorbar;
%c.Label.String = 'Att. [db/cm/MHz]';

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t4,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

end

%%
figure('Units','centimeters', 'Position',[5 5 15 5])
tiledlayout(1,2)
t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = 'BS log ratio [dB]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';

% %% Axial profiles
% figure('Units','centimeters', 'Position',[5 5 12 12])
% tiledlayout(2,1)
% nexttile,
% plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
% hold on
% plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
% plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
% plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
% plot(z_ACS,mean(attIdeal(:,20:27),2), 'k--')
% hold off
% grid on
% ylim([0.4 1.4])
% xlim([z_ACS(1) z_ACS(end)])
% title('Axial profiles')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 
% 
% nexttile,
% plot(z_ACS, axialTV{3}, 'r:', 'LineWidth',1.5),
% hold on
% plot(z_ACS, axialSWTV{3}, 'r', 'LineWidth',1),
% plot(z_ACS, axialTVL1{3}, 'b:', 'LineWidth',1.5),
% plot(z_ACS, axialWFR{3}, 'b', 'LineWidth',1),
% plot(z_ACS,mean(attIdeal(:,20:27),2), 'k--')
% hold off
% grid on
% ylim([0.4 1.4])
% xlim([z_ACS(1) z_ACS(end)])
% title('Axial profiles')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

%%
save_all_figures_to_directory(figDir,['simInc',num2str(iAcq),'Figure']);
close all

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Top')
disp(results1.biasTop)
disp(results2.biasTop)
disp(results3.biasTop)
disp(results4.biasTop)

disp('Bias Bottom')
disp(results1.biasBottom)
disp(results2.biasBottom)
disp(results3.biasBottom)
disp(results4.biasBottom)

disp('RMSE Top')
disp(results1.rmseTop)
disp(results2.rmseTop)
disp(results3.rmseTop)
disp(results4.rmseTop)

disp('RMSE Bottom')
disp(results1.rmseBottom)
disp(results2.rmseBottom)
disp(results3.rmseBottom)
disp(results4.rmseBottom)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)

