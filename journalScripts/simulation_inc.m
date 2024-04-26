% ========================================================================
% ========================================================================
% ========================================================================
%% SIMULATION ON INCLUSIONS
clear,clc
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\Simulation_23_12_18'];
baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\simulations_processed\24_04_04_inc'];

targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref2'];
resultsDir = 'C:\Users\smerino.C084288\Pictures\JOURNAL\24-04-26';
[~,~,~] = mkdir(resultsDir);

targetFiles = dir([targetDir,'\rf*.mat']);
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
attRange = [0.4,1.1];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

% Region for attenuation imaging
x_inf = -1.5; x_sup = 1.5;
z_inf = 0.4; z_sup = 3.7;
%% Setting up

for iAcq = 1:3
% Regularization
switch iAcq
    case 1
        muBtv = 10^3.5; muCtv = 10^3;
        muBswtv = 10^3; muCswtv = 10^3;
        muBtvl1 = 10^3.5; muCtvl1 = 10^3;
        muBwfr = 10^3.5; muCwfr = 10^3;
    case 2
        muBtv = 10^3; muCtv = 10^1;
        muBswtv = 10^2.5; muCswtv = 10^0;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^3.5; muCwfr = 10^1;
    case 3
        muBtv = 10^3; muCtv = 10^1;
        muBswtv = 10^2.5; muCswtv = 10^0.5;
        muBtvl1 = 10^3; muCtvl1 = 10^0.5;
        muBwfr = 10^4; muCwfr = 10^2;
end

load(fullfile(targetDir,targetFiles(iAcq).name));

fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

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
% BW from spectrogram
[pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
meanSpectrum = mean(pxx,2);
% [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
% figure,plot(fpxx/1e6,meanSpectrum)
% xline([freq_L,freq_H]/1e6)
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% xlim([0 15])

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
rInc = 0.7;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
inclusion = (Xq.^2 + (Zq-2).^2)<= (rInc-0.1)^2;
back = (Xq.^2 + (Zq-2).^2) >= (rInc+0.1)^2;
attIdeal = ones(size(Xq))*groundTruthBack(iAcq);
attIdeal((Xq.^2 + (Zq-2).^2)<= rInc^2) = groundTruthInc(iAcq);

inclusionACS = (X.^2 + (Z-2).^2)<= rInc^2;
attIdealACS{iAcq} = ones(size(X))*groundTruthBack(iAcq);
attIdealACS{iAcq}(inclusionACS) = groundTruthInc(iAcq); %incl = inclusion

%% TV
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
BRTV = reshape(Bn*NptodB,m,n);
CRTV = reshape(Cn*NptodB,m,n);

axialTV{iAcq} = mean(BRTV(:,20:27),2);

AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
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
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
MetricsSWTV(iAcq) = r;

%% TVL1
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

axialTVL1{iAcq} = mean(BRTVL1(:,20:27),2);
AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
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
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inclusion),"omitnan");
r.stdInc = std(AttInterp(inclusion),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inclusion) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inclusion) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
MetricsWFR(iAcq) = r;
%% Plotting
figure('Units','centimeters', 'Position',[5 5 22 4]);
tiledlayout(1,6, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'westoutside');
c.Label.String = '[db]';
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
title('TV')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRSWTV, attRange)
colormap(t1,turbo)
axis image
title('SWTV')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t4 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t4,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [db/cm/MHz]';
% ylabel('Axial [cm]')
xlabel('Lateral [cm]')

fontsize(gcf,8,'points')

%%
figure('Units','centimeters', 'Position',[5 5 10 8])
tl = tiledlayout(2,2, "Padding","tight", 'TileSpacing','tight');
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
c = colorbar(t1, 'eastoutside');
c.Label.String = '[dB]';
title('Bmode')
xlabel('Lateral [cm]')
ylabel('Axial [cm]')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('\DeltaBSC')
c = colorbar(t2, 'eastoutside');
c.Label.String = '[dB/cm]';
xlabel('Lateral [cm]')
%ylabel('Axial [cm]')

t3 = nexttile([1,2]); 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar(t3, 'eastoutside');
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
c.Label.String = '[a.u.]';

fontsize(gcf,8,'points')

%%
axialTV{iAcq} = mean(BRTV(:,41:49),2);
axialSWTV{iAcq} = mean(BRSWTV(:,41:49),2);
axialTVL1{iAcq} = mean(BRTVL1(:,41:49),2);
axialWFR{iAcq} = mean(BRWFR(:,41:49),2);

lateralTV{iAcq} = mean(BRTV(24:26,:),1);
lateralSWTV{iAcq} = mean(BRSWTV(24:26,:),1);
lateralTVL1{iAcq} = mean(BRTVL1(24:26,:),1);
lateralWFR{iAcq} = mean(BRWFR(24:26,:),1);
%%
save_all_figures_to_directory(resultsDir,['simInc',num2str(iAcq),'Figure']);
close all

end


%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Back')
disp(results1.biasBack)
disp(results2.biasBack)
disp(results3.biasBack)
disp(results4.biasBack)

disp('Bias Inc')
disp(results1.biasInc)
disp(results2.biasInc)
disp(results3.biasInc)
disp(results4.biasInc)

disp('RMSE Back')
disp(results1.rmseBack)
disp(results2.rmseBack)
disp(results3.rmseBack)
disp(results4.rmseBack)

disp('RMSE Inc')
disp(results1.rmseInc)
disp(results2.rmseInc)
disp(results3.rmseInc)
disp(results4.rmseInc)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)



T = [results1;results2;results3;results4];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);

%% Lateral and axial profiles
layered = load('simuLayered.mat');
figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile([1,2]),
plot(layered.z_ACS, layered.axialTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(layered.z_ACS, layered.axialSWTV{1}, 'r', 'LineWidth',1),
plot(layered.z_ACS, layered.axialTVL1{1}, 'b:', 'LineWidth',1.5),
plot(layered.z_ACS, layered.axialWFR{1}, 'b', 'LineWidth',1),
plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([layered.z_ACS(1) layered.z_ACS(end)])
%title('Axial profile')
xlabel('Axial [cm]')
ylabel('ACS [dB/cm/MHz]')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northwest') 

figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile([1,2]),
plot(layered.z_ACS, layered.axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(layered.z_ACS, layered.axialSWTV{2}, 'r', 'LineWidth',1),
plot(layered.z_ACS, layered.axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(layered.z_ACS, layered.axialWFR{2}, 'b', 'LineWidth',1),
plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([layered.z_ACS(1) layered.z_ACS(end)])
%title('Axial profile')
xlabel('Axial [cm]')
ylabel('ACS [dB/cm/MHz]')
legend({'TV','SWTV','TVL1','WFR'}, 'Location','northwest') 

figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile([1,2]),
plot(layered.z_ACS, layered.axialTV{3}, 'r:', 'LineWidth',1.5),
hold on
plot(layered.z_ACS, layered.axialSWTV{3}, 'r', 'LineWidth',1),
plot(layered.z_ACS, layered.axialTVL1{3}, 'b:', 'LineWidth',1.5),
plot(layered.z_ACS, layered.axialWFR{3}, 'b', 'LineWidth',1),
plot(layered.z_ACS,mean(layered.attIdeal,2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([layered.z_ACS(1) layered.z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')
%legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

%%
figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
% figure('Units','centimeters', 'Position',[5 5 12 12])
nexttile,
plot(z_ACS, axialTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{1}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{1}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{1}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdealACS{1}(:,41:49),2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')

nexttile,
plot(x_ACS, lateralTV{1}, 'r:', 'LineWidth',1.5),
hold on
plot(x_ACS, lateralSWTV{1}, 'r', 'LineWidth',1),
plot(x_ACS, lateralTVL1{1}, 'b:', 'LineWidth',1.5),
plot(x_ACS, lateralWFR{1}, 'b', 'LineWidth',1),
plot(x_ACS,mean(attIdealACS{1}(24:26,:),1), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([x_ACS(1) x_ACS(end)])
%title('Lateral profile')
xlabel('Lateral [cm]')


figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile,
plot(z_ACS, axialTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{2}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{2}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{2}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdealACS{2}(:,41:49),2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')

nexttile,
plot(x_ACS, lateralTV{2}, 'r:', 'LineWidth',1.5),
hold on
plot(x_ACS, lateralSWTV{2}, 'r', 'LineWidth',1),
plot(x_ACS, lateralTVL1{2}, 'b:', 'LineWidth',1.5),
plot(x_ACS, lateralWFR{2}, 'b', 'LineWidth',1),
plot(x_ACS,mean(attIdealACS{2}(24:26,:),1), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([x_ACS(1) x_ACS(end)])
% title('Lateral profile')
xlabel('Lateral [cm]')


figure('Units','centimeters', 'Position',[5 5 10 4])
tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
nexttile,
plot(z_ACS, axialTV{3}, 'r:', 'LineWidth',1.5),
hold on
plot(z_ACS, axialSWTV{3}, 'r', 'LineWidth',1),
plot(z_ACS, axialTVL1{3}, 'b:', 'LineWidth',1.5),
plot(z_ACS, axialWFR{3}, 'b', 'LineWidth',1),
plot(z_ACS,mean(attIdealACS{3}(:,41:49),2), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([z_ACS(1) z_ACS(end)])
%title('Axial profile')
ylabel('ACS [dB/cm/MHz]')
xlabel('Axial [cm]')

nexttile,
plot(x_ACS, lateralTV{3}, 'r:', 'LineWidth',1.5),
hold on
plot(x_ACS, lateralSWTV{3}, 'r', 'LineWidth',1),
plot(x_ACS, lateralTVL1{3}, 'b:', 'LineWidth',1.5),
plot(x_ACS, lateralWFR{3}, 'b', 'LineWidth',1),
plot(x_ACS,mean(attIdealACS{3}(24:26,:),1), 'k--')
hold off
grid on
ylim([0.4 1.1])
xlim([x_ACS(1) x_ACS(end)])
% title('Lateral profile')
xlabel('Lateral [cm]')
% legend({'TV','SWTV','TVL1','WFR'}, 'Location','northeastoutside') 

fontsize(gcf,8,'points')

save_all_figures_to_directory(resultsDir,'SimAxialLat');

close all