clear,clc
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\', ...
    'Attenuation\DataQUS_4_Merino'];rfDir = [baseDir,'\Hashimoto'];
rfFiles = dir([rfDir,'\*.mat']);
%%
for iFile = 1:length(rfFiles)
    %%
    iFile = 8;
    load([rfDir,'\',rfFiles(iFile).name]);
    
    % Frequency response rf - 3 to 9 MHz
    %freqResponse = median(abs(fft(RF)),2);
    %f = (1:size(RF,1)) /size(RF,1) * fs;
    %plot(f/1e6,freqResponse)
    
    % Plotting Bmode
    dx = x(2)-x(1);
    dz = z(2)-z(1);
    x = x*1e2; % [cm]
    z = z*1e2; % [cm]
    dynRange = [-60,-10];
    
    figure, imagesc(x,z,Bmode,dynRange)
    colormap gray
    colorbar
    axis equal tight
end
%% LOADING PARAMETERS
blocksize = 20;     % Block size in wavelengths
c0 = 1540;
freq_L = 3; freq_H = 9;
overlap_pc      = 0.8;
winsize         = 0.5;

% Region for attenuation imaging
x_inf = 0.5; x_sup = 3;
z_inf = 0.5; z_sup = 3;

sam1 = RF;

% Limits for ACS estimation
ind_x = x_inf <= x & x <= x_sup;
x = x(ind_x);
ind_z = z_inf <= z & z <= z_sup;
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Freq limits in Hz
freq_L = freq_L*1e6;   % (Hz)
freq_H = freq_H*1e6;   % (Hz)

wl = c0/mean([freq_L freq_H]);   % Wavelength (m)

% Number of repetitions of datablocks within a single datablock
rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
rpt = round(rpt);
overlap = 1 - (1/rpt);


% Lateral samples between windows
wx  = round((blocksize*wl)/(dx*rpt));

% Number of lines laterally = repetitions * block without repetition
nx  = rpt*wx;
new_blocksize = round(nx*dx/(wl));

L2   = size(sam1,2);                % RF data columns
ncol = floor((L2-(rpt-1)*wx)/wx);   % Blocksize colums
sam1 = sam1(:,1:wx*(ncol+rpt-1));
L2 = size(sam1,2);                  % Actual RF data columns
x  = x(1:L2);
xi = 1; xf = L2;

x0 = (xi:wx:xf+1-nx);
x_ACS = x(1,x0+round(nx/2));
n  = length(x0);


% Axial samples between windows
wz = floor(nx*dx/(dz*rpt));
nz = rpt*wz;

% winsize: Percentage of window (max 0.5)
% nw: Samples of each window axially
nw = 2*floor(winsize*nx*dx/(2*dz)) - 1 ;
L = (nz - nw)*dz*100;   % (cm)

NFFT = 2^(nextpow2(nw)+2);
band = fs*linspace(0,1,NFFT)';   % [Hz] Band of frequencies
rang = (floor(freq_L/fs*NFFT)+1:round(freq_H/fs*NFFT));   % useful frequency range
f  = band(rang)*1e-6; % [MHz]
p = length(rang);

L1   = size(sam1,1);                % RF data: rows
nrow = floor((L1-(rpt-1)*wz)/wz);   % Blocksize rows
sam1 = sam1(1:wz*(nrow+rpt-1),:);
L1   = size(sam1,1);                % RF data: rows
z    = z(1:L1);
zi = 1; zf = L1;

z0 = (zi:wz:zf+1-nz);
z_ACS = z(z0+round(nz/2));
m  = length(z0);

z0p = z0 + (nw-1)/2;
z0d = z0 + (nz-1) - (nw-1)/2;


% Plot region of interest B-mode image
Im=abs(hilbert(sam1));   % envelope calculation
Bmode=20*log10(Im/max(Im(:)));   % log scale
figure, imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
hb2=colorbar; ylabel(hb2,'dB')
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');

disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);


%% CLUSTERING
figure('Units','centimeters', 'Position',[5 5 30 8]), 
tiledlayout(1,3)
nexttile,
imagesc(x,z,Bmode, dynRange)
colormap gray
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Bmode')

h = fspecial("average",[50 5]);
blurred = imfilter(Bmode,h,"symmetric"); 
nexttile, imagesc(x,z,blurred,dynRange);
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Blurred image')

[X,Z] = meshgrid(x,z);
Data = normalize([X(:) Z(:) blurred(:)]);
Data = Data.*[1 1 1];
[mBm,nBm] = size(Bmode);
nClusters = 7; % number of clusters
idx = kmeans(Data,nClusters);
ID = reshape(idx,[mBm,nBm]);

nexttile;
segmented = labeloverlay(reshape(normalize(Bmode(:),'range'),[mBm,nBm]),ID);
imagesc(x,z,segmented)
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Segmented image')

%% Equalizing
mask = ID==1;
energy1 = std(sam1(mask));
factor = ones(size(sam1));

for iCluster = 1:nClusters
    mask = ID==iCluster;
    factor(mask) = energy1/std(sam1(mask));
end
h = fspecial("average",[50 5]);
factor = imfilter(factor,h,"symmetric");

%
figure('Units','centimeters', 'Position',[5 5 30 8]), 
tiledlayout(1,2)
t2 = nexttile;
imagesc(x,z,factor)
colormap(t2,parula)
colorbar
axis image
title('Factor')

samEnhanced = sam1.*factor;
%samEnhanced = sam1;
Bmode2 = db(hilbert(samEnhanced));
Bmode2 = Bmode2 - max(Bmode2(:));
t3 = nexttile;
imagesc(x,z,Bmode2,dynRange)
colormap(t3,gray)
colorbar
axis image
title('Equalized B-mode')

%% Calculating spectra
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);
Nw = 400;
Noverlap = 200;
Sp = zeros(m,n,length(f));
Sd = zeros(m,n,length(f));
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = samEnhanced(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = samEnhanced(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

jj = floor(3*n/6);
ii = ceil(m/6);
figure(3)
plot(f,10*log10(squeeze(Sp(ii,jj,:)/max(Sd(ii,jj,:)))),'k');
hold on
plot(f,10*log10(squeeze(Sd(ii,jj,:)/max(Sd(ii,jj,:)))),'r');
hold off
title('Spectrum'); xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
axis([f(1) f(end) -30 10]);

%% Au = b
load([baseDir,'\References\P4_CUELLO_3.mat']);
b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

[u,~] = cgs2(A'*A,A'*b(:),1e-6,20);

% Standard SLD
% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); %CS = u(end/2+1:end);
BS = 8.686*BS;   % [dB.cm^{-1}.MHz^{-1}]
BS = reshape(BS,m,n);

imagesc(x_ACS,z_ACS,BS, [0 1]);
xlabel('x[mm]'), ylabel('z[mm]'),
axis equal tight
colormap turbo
c = colorbar;
c.Label.String = 'Attenuation [db/cm/MHz]';
title('SLD')

%% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = 1e4*[0.4,1.2,3.6];

for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


% Plotting
for ii = 1:size(BR,3)
    figure('Position',[100 100 900 300]), 
    tiledlayout(1,3)
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis equal
    xlim([x_ACS(1) x_ACS(end)])
    ylim([z_ACS(1) z_ACS(end)])
    colormap(t3,gray)
    colorbar

    t1 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), [0 1.5])
    colormap(t1,turbo)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Attenuation';
    title(['RSLD, \mu=',num2str(mu(ii),2)])


    t2 = nexttile;
    imagesc(x_ACS,z_ACS,CR(:,:,ii))
    colormap(t2,parula)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Backscatter';
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end