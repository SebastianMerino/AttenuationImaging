clc, clear,

baseDir = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\DeepLearning';

targetDir = fullfile(baseDir,'raw');
targetFiles = dir(fullfile(targetDir,'*.mat'));
for iAcq = 1:length(targetFiles)
    fprintf("Acquisition no. %i, %s\n",iAcq,targetFiles(iAcq).name);
end 

refDir = fullfile(baseDir,'ref');
refFiles = dir(fullfile(refDir,'*.mat'));

%% Generating cropped data
% SETTING PARAMETERS
blocksize = 8;     % Block size in wavelengths
freq_L = 3e6; freq_H = 8.5e6; % GOOD
freqCutOff = db2mag(-30);

overlap_pc      = 0.8;
ratio_zx        = 12/8;
referenceAtt    = 0.5;

x_inf = -2; x_sup = 2;
z_inf = 0.2; z_sup = 3.5;

% Weight parameters
muB = 10^3; muC = 10^0;
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;


% Plotting
dynRange = [-40,0];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;
attRange = [0.5 1.4];
% attRange = [0.6 1.7];

%% For looping
iAcq = 6;

load(fullfile(targetDir,targetFiles(iAcq).name));
fprintf("Acquisition no. %i, %s\n",iAcq,targetFiles(iAcq).name);

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
sam1 = rf(:,:,1);
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

% selecting BW 
windowing = tukeywin(nz/2,0.25);
NFFT = 2^(nextpow2(nz/2)+2);
% [pxx,fpxx] = pwelch(sam1-mean(sam1),500,400,500,fs); % round(nz*overlap_pc)
[pxx,fpxx] = pwelch(sam1-mean(sam1),windowing,round(nz/4),NFFT,fs);
meanSpectrum = mean(pxx,2);
meanSpectrum = meanSpectrum./max(meanSpectrum);
figure,plot(fpxx/1e6,db(meanSpectrum))
[freq_L,freq_H] = findFreqBand(fpxx, meanSpectrum, freqCutOff);
xline([freq_L,freq_H]/1e6)
xlim([0 10])
ylim([-40 0])
xlabel('Frequency [MHz]')
ylabel('Magnitude')
grid on

% Frequency samples
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
if true
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

%% ROI selection
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);

% Setting up
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);


% %% SLD
% 
% sld1 = squeeze(sum(sum(b.*incAcs,1),2))/sum(incAcs(:)) * NptodB /4/L;
% acs1 = [f,ones(size(f))]\sld1;
% fprintf('Attenuation is %.2f\n',acs1(1))
% sld2 = squeeze(sum(sum(b.*backAcs,1),2))/sum(backAcs(:)) * NptodB /4/L;
% acs2 = [f,ones(size(f))]\sld2;
% fprintf('Attenuation is %.2f\n',acs2(1))
% figure('Units','centimeters', 'Position',[5 5 15 15])
% plot(f,sld1)
% hold on
% plot(f,sld2)
% plot(f,acs1(1)*f + acs1(2), 'k--')
% plot(f,acs2(1)*f + acs2(2), 'k--')
% hold off
% grid on,
% xlim([0,max(f)]), ylim([-3 10]),
% xlabel('Frequency [MHz]')
% ylabel('Att. [dB/cm]')
% title('Mean SLD')
% legend('Inc','Back')



%% RSLD-TV
muBtv = 10^3.5; muCtv = 10^3;

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn*NptodB,m,n));

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'Isotropic RSLD')
t2 = nexttile; 
imagesc(x_ACS,z_ACS,BR, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBtv,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
t3 = nexttile; 
imagesc(x_ACS,z_ACS,CR, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCtv,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';



%% TVL1
muBtvl1 = 10^3.5; muCtvl1 = 10^2;

tic
[Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBtvl1,muCtvl1,m,n,tol,mask(:));
toc
BRTik = (reshape(Bn*NptodB,m,n));
CRTik = (reshape(Cn*NptodB,m,n));

figure('Units','centimeters', 'Position',[5 5 15 6]);
tl = tiledlayout(1,2, "Padding","tight");
title(tl,'Isotropic RSLD')
t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRTik, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBtvl1,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRTik, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCtvl1,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';
 
%% NEW WEIGHTS
muBwfr = 10^3.5; muCwfr = 10^2;

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
CRWTik = (reshape(Cn*NptodB,m,n));

figure('Units','centimeters', 'Position',[5 5 20 6]);
tl = tiledlayout(1,3, "Padding","tight");
title(tl,'WFR')
t1 = nexttile; 
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t1,parula)
axis image
title('Weights')
c = colorbar;
c.Label.String = '[a.u.]';

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRWTik, attRange)
colormap(t2,turbo)
axis image
title(['RSLD, \mu=',num2str(muBwfr,2)])
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';
t3 = nexttile; 
imagesc(x_ACS,z_ACS,CRWTik, bsRange)
colormap(t3,parula)
axis image
title(['RSLD, \mu=',num2str(muCwfr,2)])
c = colorbar;
c.Label.String = 'BS log ratio [dB]';


%%
figure('Units','centimeters', 'Position',[5 5 20 4]);
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
fontsize(gcf,8,'points')

t2 = nexttile;
imagesc(x_ACS,z_ACS,BR, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title('TV')
fontsize(gcf,8,'points')

t4 = nexttile;
imagesc(x_ACS,z_ACS,BRTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t4,turbo)
axis image
title('TVL1')
fontsize(gcf,8,'points')

t5 = nexttile;
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(t5,turbo)
axis image
title('WFR')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
