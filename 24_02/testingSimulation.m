clear,clc
addpath('./functions_v7');
addpath('./AttUtils');

% baseDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\simulations_processed\24_01_26'];
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\24_02_05'];
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\24_01_30'];
targetDir = [baseDir,'\raw'];
refDir = [baseDir,'\ref'];
resultsDir = fullfile(baseDir,'results','24-02-07','BS_12_12');
if ~exist(resultsDir,'dir'); mkdir(resultsDir); end


targetFiles = dir([targetDir,'\rf*.mat']);
refFiles = dir([refDir,'\rf*.mat']);

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 12;     % Block size in wavelengths
freq_L = 3e6; freq_H = 8.5e6; % GOOD

overlap_pc      = 0.8;
ratio_zx        = 1;
referenceAtt    = 0.6;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% groundTruthBack = [0.6,1.5];
% groundTruthInc = [1.2,0.8];
% groundTruthBack = [1.5,1.5,0.8,0.8];
% groundTruthInc = [0.8,0.8,1.5,1.5];
groundTruthBack = [1.2,1.2,0.7,0.7];
groundTruthInc = [0.7,0.7,1.2,1.2];

% Plotting
dynRange = [-40,0];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;
attRange = [0.5 1.4];
% attRange = [0.6 1.7];

%% For looping
for iAcq = 5:8
iAcq = 5;

load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("Acquisition no. %i, patient %s\n",iAcq,targetFiles(iAcq).name);
iAcq = iAcq-4;

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
sam1 = rf(:,:,1);

% Region for attenuation imaging
BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));
xFull = x; zFull = z;
dynRange = [-50,-10];
figure('Units','centimeters', 'Position',[5 5 15 15]),
imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
ylim([0.5 3.5])
confirmation = '';
while ~strcmp(confirmation,'Yes')
    rect = getrect;
    confirmation = questdlg('Sure?');
    if strcmp(confirmation,'Cancel')
        disp(rect)
        break
    end
end
close,
x_inf = rect(1); x_sup = rect(1)+rect(3); 
z_inf = rect(2); z_sup = rect(2)+rect(4);

x_inf = -0.5; x_sup = 0.5;
z_inf = 1.5; z_sup = 2.5;

% x_inf = -1.5; x_sup = 1.5;
% z_inf = 0.5; z_sup = 3.5;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
ind_z = z_inf <= z & z <= z_sup;
x = x(ind_x);
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);
%%
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


% BW from spectrogram
% ratio = db2mag(-30);
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs);
% meanSpectrum = mean(pxx,2);
% % [freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, ratio);
% figure,plot(fpxx/1e6,meanSpectrum)
% xline([freq_L,freq_H]/1e6)
% xlabel('Frequency [MHz]')
% ylabel('Magnitude')
% xlim([0 15])
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
% figure, imagesc(x,z,Bmode);
% axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


fprintf('\nFrequency range: %.2f - %.2f MHz\n',freq_L*1e-6,freq_H*1e-6)
fprintf('Blocksize in wavelengths: %i\n',blocksize)
fprintf('Blocksize x: %.2f mm, z: %.2f mm\n',nx*dx*1e3,nz*dz*1e3)
fprintf('Blocksize in pixels nx: %i, nz: %i\n',nx,nz);
fprintf('Region of interest columns: %i, rows: %i\n\n',m,n);

%% Generating Diffraction compensation
if iAcq ==1
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
    tic
    load(fullfile(refDir,refFiles(iRef).name),"rf");
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
            % [~,tempSp] = spectra2(sub_block_p,windowing,0,nz/2,NFFT);
            % [~,tempSd] = spectra2(sub_block_d,windowing,0,nz/2,NFFT);

            Sp_ref(ii,jj,:,iRef) = (tempSp(rang));
            Sd_ref(ii,jj,:,iRef) = (tempSd(rang));
        end
    end
    toc
%     figure,imagesc(Sd_ref(:,:,40,iRef))

end

Sp = mean(Sp_ref,4); Sd = mean(Sd_ref,4);
compensation = ( log(Sp) - log(Sd) ) - 4*L*att_ref_map;
% compensation = ( Sp - Sd ) - 4*L*att_ref_map;

% Liberating memory to avoid killing my RAM
clear Sp_ref Sd_ref
end
% %%
% diffraction_xz = mean(compensation,3);
% diffraction_zf = squeeze(mean(compensation,2));
% figure, tiledlayout(1,2)
% nexttile,
% imagesc(x_ACS,z_ACS,diffraction_xz, [-1 1]);
% title('Diffraction compensation'),
% xlabel('x [cm]'), ylabel('z [cm]'),
% colorbar
% nexttile,
% imagesc(f,z_ACS,diffraction_zf, [-1 1]);
% title('Diffraction compensation'),
% xlabel('f [MHz]'), ylabel('z [cm]'),
% colorbar


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
        % [~,tempSp] = spectra2(sub_block_p,windowing,0,nz/2,NFFT);
        % [~,tempSd] = spectra2(sub_block_d,windowing,0,nz/2,NFFT);
        
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

%% Setting Up

% System of equations
b = (log(Sp) - log(Sd)) - (compensation);
% b = ((Sp) - (Sd)) - (compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];
tol = 1e-3;
clear mask
mask = ones(m,n,p);

% 

%% Ideal
rInc = 0.7;
cz = 2; cx = 0;
% 
% [X,Z] = meshgrid(x_ACS,z_ACS);
% attIdeal = ones(size(Z));
% inc = ((X-cx).^2 + (Z-cz).^2)<= rInc^2;
% attIdeal(~inc) = groundTruthBack(iAcq);
% attIdeal(inc) = groundTruthInc(iAcq); %incl = bottom
% 
% figure('Units','centimeters', 'Position',[5 5 20 6]);
% tiledlayout(1,2)
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis equal
% xlim([x_ACS(1) x_ACS(end)])
% ylim([z_ACS(1) z_ACS(end)])
% colormap(t1,gray)
% colorbar('westoutside')
% title('Bmode')
% 
% t2 = nexttile;
% imagesc(x,z,attIdeal,attRange)
% axis image
% colormap(t2,turbo);
% c = colorbar('westoutside');
% c.Label.String = 'dB/cm/MHz';
% title('Ideal ACS')
% 
% 
[X,Z] = meshgrid(x,z);
inc = ((X-cx).^2 + (Z-cz).^2)<= (rInc-0.05)^2;
back = ((X-cx).^2 + (Z-cz).^2) >= (rInc+0.05)^2;
% 
% %% Spectrum
% % % Heterogeneous
% [X,Z] = meshgrid(x_ACS,z_ACS);
% region1 = (X.^2 + (Z-2).^2) <= (rInc-0.1)^2;
% region2 = (X.^2 + (Z-2).^2) >= (rInc+0.1)^2;
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
% legend('Inc','Back')

%% RSLD-TV
muBtv = 10^3; muCtv = 10^1;
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));


% Overlay
[X,Z] = meshgrid(xFull,zFull);
roi = X >= x_ACS(1) & X <= x_ACS(end) & Z >= z_ACS(1) & Z <= z_ACS(end);
%figure, imagesc(roi);

figure('Units','centimeters', 'Position',[5 5 10 15])
tiledlayout(2,1)
t2 = nexttile;
imagesc(xFull,zFull,BmodeFull,[-50 0])
axis image
title('B-mode')
ylim([0.5, 3.5])
c = colorbar;
c.Label.String = 'dB';
hold on
% contour(xFull,zFull,roi,1,'w--')
% contour(x,z,maskThyroid,1,'w--')
hold off

nexttile,
[~,hB,hColor] = imOverlayInterp(BmodeFull,BR,[-50 0],attRange,0.5,...
    x_ACS,z_ACS,roi,xFull,zFull);
title('B-mode and attenuation map')
hColor.Label.String = 'dB/cm/MHz';
ylim([0.5, 3.5])
hold on
% contour(xFull,zFull,roi,1,'w--')
% contour(x,z,maskThyroid,1,'w--')
hold off

colormap(t2,gray)


%% TVL1
muBtvl1 = 1e3; muCtvl1 = 1;
% muBtvl1 = 3e3; muCtvl1 = 3;
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
BRTVL1 = reshape(Bn*NptodB,m,n);
CRTVL1 = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 8]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1)
title('Bmode')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
colormap(t1,turbo)
axis image
title('TVL1')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';

t4 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t4,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = '\DeltaBSC [dB/cm]';

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
MetricsTVL1(iAcq) = r;

%% WFR
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Computing weights
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

% w = ones(size(CRTVL1));
% borderMask = (Z-cz).^2 + (X-cx).^2 <= (rInc-0.1)^2 | ...
%     (Z-cz).^2 + (X-cx).^2 >= (rInc+0.1)^2 ;
% w(borderMask) = 1;
% w(~borderMask) = 0.01;

figure('Units','centimeters', 'Position',[5 5 24 8])
tl = tiledlayout(1,3, "Padding","tight");
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1)
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,CRTVL1, bsRange)
colormap(t2,parula)
axis image
title('TVL1')
c = colorbar;
c.Label.String = 'BS log ratio [dB/cm]';

t3 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
axis image
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';
%% -------
% Setting up new system
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muBwfr = 3e3; muCwfr = 3;
% muBwfr = 1e4; muCwfr = 10;
% muBwfr = 1e3; muCwfr = 1;
% Method
[Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBwfr,muCwfr,m,n,tol,mask(:),w);
BRWFR = reshape(Bn*NptodB,m,n);
CRWFR = reshape(Cn*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 8]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis equal
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
colormap(t1,gray)
colorbar(t1)
title('Bmode')

t1 = nexttile; 
imagesc(x_ACS,z_ACS,BRWFR, attRange)
colormap(t1,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';


t4 = nexttile; 
imagesc(x_ACS,z_ACS,CRWFR, bsRange)
colormap(t4,parula)
axis image
title('WFR')
c = colorbar;
c.Label.String = '\DeltaBSC [dB/cm]';

AttInterp = interp2(X,Z,BRWFR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
MetricsWFR(iAcq) = r;




end
%%
tableName = 'results.xlsx';
T = [struct2table(MetricsTVL1);struct2table(MetricsWFR)];
writetable(T,fullfile(resultsDir,tableName),...
     'WriteRowNames',true);

save_all_figures_to_directory(resultsDir,'fig')
close all,
% 

