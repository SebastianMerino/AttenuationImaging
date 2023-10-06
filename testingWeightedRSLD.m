%%
% Background 0.55
% T1: 0.52
% T4: 0.81
% T8: 0.95
clear,clc
addpath('./functions_att');

rfDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% rfDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];
rawFiles = dir([rfDir,'\*.rf']);

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];
refFiles = dir([refDir,'\*.rf']);

disp('Background 0.55, T1: 0.52, T4: 0.81, T8: 0.95')
% %% Generating .mat data
% for iAcq = 1:length(rawFiles)
%     [out]=lectura_OK([rfDir,'\',rawFiles(iAcq).name]);
%     save([rfDir,'\ref',num2str(iAcq),'.mat'],'-struct','out');
% end

%% Selecting ROI limits and freq limits
iAcq = 4;
disp(['Target: ', rawFiles(iAcq).name]);
load([rfDir,'\ref',num2str(iAcq),'.mat'])
attRange = [0.3,1.2];

% Inclusion
% cx = 1.92; cz = 1.98; r = 0.98; % T1
cx = 1.90; cz = 1.87; r = 0.96; % T4
% cx = 1.85; cz = 1.93; r = 0.93; % T8

dx = x(2)-x(1);
dz = z(2)-z(1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]

sam1 = RF(:,:,1);
dynRange = [-50,0];

% Plotting
% [pSpectrum,f] = pwelch(RF,hamming(500),300,1024,fs);
% meanSpectrum = mean(pSpectrum,2);
% figure, plot(f,meanSpectrum)
% xlim([0 fs/2])
% title('Mean power spectrum')
% 
% Bmode = db(hilbert(sam1));
% Bmode = Bmode - max(Bmode(:));
% figure, imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');


%% SETTING PARAMETERS
blocksize = 20;     % Block size in wavelengths
c0 = 1540;
freq_L = 2; freq_H = 9;
overlap_pc      = 0.8;
winsize         = 0.5;
% Region for attenuation imaging
x_inf = 0.2; x_sup = 3.8;
z_inf = 0.2; z_sup = 3.5;

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
Bmode = db(hilbert(sam1));
Bmode = Bmode - max(Bmode(:));
% figure, imagesc(x,z,Bmode); axis image; colormap gray; clim(dynRange);
% hb2=colorbar; ylabel(hb2,'dB')
% xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');

disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);

%% Generating references
x1 = x; z1 = z;
for iRef = 1:length(refFiles)
    load([refDir,'\ref',num2str(iRef),'.mat'],'RF','x','z');
    x = x*1e2; z = z*1e2; % [cm]
    ind_x = x_inf <= x & x <= x_sup;
    x = x(ind_x);
    ind_z = z_inf <= z & z <= z_sup;
    z = z(ind_z);
    x = x(1:L2);
    z = z(1:L1);
    fprintf("Checking difference, x: %f, z: %f\n",norm(x1-x), norm(z1-z))
    samRef = RF;
    samRef = samRef(ind_z,ind_x);
    samRef = samRef(:,1:wx*(ncol+rpt-1));
    samRef = samRef(1:wz*(nrow+rpt-1),:);
    ref{iRef} = samRef;
end

%% Compensation
disp('Calculating diffraction compensation from references...')

windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
windowing = windowing*ones(1,nx);

%att_ref = attenuation_phantoms_Np(f, 3, []);
att_ref = 0.53*f/8.686; % From phantom especifications

att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

% Four 'samples' of the reference phantom
Sp_ref = zeros(m,n,p);
Sd_ref = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);
        % Reference 1
        sub_block_p = ref{1}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{1}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp1,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd1,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 2
        sub_block_p = ref{2}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{2}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp2,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd2,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 3
        sub_block_p = ref{3}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{3}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp3,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd3,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        % Reference 4
        sub_block_p = ref{4}(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = ref{4}(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        [tempSp4,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd4,~] = spectra(sub_block_d,windowing,0,nw,NFFT);

        tempSp = 1/4*(tempSp1 + tempSp2 + tempSp3 + tempSp4);
        tempSd = 1/4*(tempSd1 + tempSd2 + tempSd3 + tempSd4);
        Sp_ref(ii,jj,:) = (tempSp(rang));
        Sd_ref(ii,jj,:) = (tempSd(rang));
    end
end

diffraction_compensation = ( log(Sp_ref) - log(Sd_ref) ) - 4*L*att_ref_map;
disp('Done')
%% Calculating spectra
windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25

% Windowing neccesary before Fourier transform
windowing = windowing*ones(1,nx);
Sp = zeros(m,n,length(f));
Sd = zeros(m,n,length(f));
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,0,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,0,nw,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
    end
end

jj = round(n/2);
ii = round(m/4);
figure('Units','centimeters', 'Position',[5 5 10 10])
plot(f,10*log10(squeeze(Sp(ii,jj,:)/max(Sd(ii,jj,:)))),'k');
hold on
plot(f,10*log10(squeeze(Sd(ii,jj,:)/max(Sd(ii,jj,:)))),'r');
hold off
title('SLD at interface');
xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
axis([f(1) f(end) -30 10]);

%% Au = b
load('T4.mat')

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];


% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = 1e4*[0.4,1.2,3.6];
%mu = 1e4*[0.4,2,10];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


%% Plotting
for ii = 1:size(BR,3)
    figure('Units','centimeters', 'Position',[5 5 30 8]), 
    tiledlayout(1,3)
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis image
    colormap(t3,gray)
    colorbar
    title('Bmode')

    t1 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
    colormap(t1,turbo)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    title(['RSLD, \mu=',num2str(mu(ii),2)])


    t2 = nexttile;
    imagesc(x_ACS,z_ACS,CR(:,:,ii))
    colormap(t2,parula)
    axis equal tight
    c = colorbar;
    c.Label.String = 'BS constant';
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end

%% ----------------- WEIGHTS FROM LOCAL SNR -----------------
% FUNCTION
SNR = 0:0.01:6;
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = (SNR/SNRopt);
legends = {};
figure('Units','centimeters', 'Position',[5 5 20 10])
for gamma = 1:10
    w = desvSNR.^gamma.*exp(1-desvSNR.^gamma);
    plot(SNR,w)
    hold on
    legends{gamma} = ['\gamma = ',num2str(gamma)];
end
xline(1.91, 'k--')
xlabel('SNR (\mu/\sigma)')
ylabel('Weights')
legend([legends,'SNR_{opt}'])
hold off

%% Getting local weights
envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        
        %temp = [sub_block_p(:) sub_block_d(:)];
        %SNR(ii,jj) = mean(temp)/std(temp);
        SNR(ii,jj) = min([mean(sub_block_p)/std(sub_block_p),...
            mean(sub_block_d)/std(sub_block_d)]);
    end
end

% Weights
figure('Units','centimeters', 'Position',[5 5 30 8]),
tiledlayout(1,3)
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
colormap(t1,gray)
colorbar
axis equal
xlim([x_ACS(1) x_ACS(end)]), ylim([z_ACS(1) z_ACS(end)]);
% axis image
title('B-mode')

t2 = nexttile;
imagesc(x_ACS,z_ACS,db(SNR))
colormap(t2,parula)
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')


SNRopt = sqrt(1/(4/pi - 1));
desvSNR = (SNRopt./SNR);
gamma = 10;
w = desvSNR.^gamma.*exp(1-desvSNR.^gamma);

t3 = nexttile;
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
c = colorbar;
axis image
title(['Weights, order=',num2str(gamma)])



%% Au = b
W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);

A1w = W*A1;
A2w = W*A2;
Aw = [A1w A2w];

% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
mu = 1e4*[0.4,1.2,3.6];
BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1w,A2w,bw,mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


%% Plotting
for ii = 1:size(BR,3)
    figure('Units','centimeters', 'Position',[5 5 30 8]), 
    tiledlayout(1,3)
    t3 = nexttile;
    imagesc(x,z,Bmode,dynRange)
    axis image
    colormap(t3,gray)
    colorbar
    title('B-mode')

    t1 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
    colormap(t1,turbo)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Att. [db/cm/MHz]';
    title(['RSLD, \mu=',num2str(mu(ii),2)])


    t2 = nexttile;
    imagesc(x_ACS,z_ACS,CR(:,:,ii))
    colormap(t2,parula)
    axis equal tight
    c = colorbar;
    c.Label.String = 'BS constant';
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end



%% ------------------------ UTILITY FUNCTIONS ------------------------
% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_ANIS_TV(b,A,mu,M,N,tol,mask,minimask)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol && ite_irls < 20

    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];

    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);

    P = Dh.^2 + Dv.^2;
    eps = 0.1;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    %W = kron(speye(2),omega);

    Px = abs(Dh + eps);
    Px = 1./Px;
    Px = Px(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Px,0,omega);
    Wx = kron(speye(1),omega);

    Py = abs(Dv + eps);
    Py = 1./Py;
    Py = Py(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(Py,0,omega);
    Wy = kron(speye(1),omega);

    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u] = cgs2( AtA + mu*Dx'*Wx*Dx + mu*Dy'*Wy*Dy , A'*b, 1e-6 , 20, u );

    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc2(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));

end
end


function [B,C] = AlterOpti_ADMM(A1,A2,b,mu1,mu2,m,n,tol,mask)

p = length(mask)/(m*n);
minimask = reshape(mask,[m n p]);
minimask = minimask(:,:,1);
minimask = minimask(:);
b = mask.*b;
b(isnan(b)) = 0;
A = [A1 A2];
[u,~] = cgs2(A'*A,A'*b,1e-6,20);
B = reshape(u(1:end/2),m,n);
C = reshape(u(end/2+1:end),m,n);
%figure(109); imagesc(8.686*reshape(B,m,n)); colormap pink; caxis([0 1.2])
B = B(:);
C = C(:);
D = 0;
v = 0;

F(1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

ite  = 0;
error = 1;

while abs(error) > tol && ite < 20
    ite = ite + 1;

    rho = 1;
    % First part of ADMM algorithm: B
    % B = IRLS_ANIS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
    B = IRLS_TV(b-A2*C-D-v,A1,mu1/rho,m,n,tol,mask,minimask);
    %F(2*ite,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite,1) - F(2*ite-1,1);

    % Second part of ADMM algorithm: C
    % C = IRLS_ANIS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    C = IRLS_TV(b-A1*B-D-v,A2,mu2/rho,m,n,tol,mask,minimask);
    %F(2*ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);
    %error = F(2*ite+1,1) - F(2*ite,1);

    % Third part of ADMM algorithm: D
    % Least squares: 1/2*||D||_2^2 + rho/2*||D-w||_2^2
    w = b - A1*B - A2*C - v;
    D = (rho/(rho+1))*w;

    % Fourth part of ADMM algorithm: v
    v = v + A1*B + A2*C + D - b;
    F(ite+1,1) = 1/2*(norm( b - A1*B - A2*C ))^2 + mu1*TVcalc2(B,m,n,minimask) + mu2*TVcalc2(C,m,n,minimask);

end

end


% TV Andres Leonel Coila
function [TV] = TVcalc(B,M,N,mask)

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
P = sqrt(P);
TV = norm(P(:).*mask,1);

end



% TV Andres Leonel Coila - ANISOTROPIC
function [TV] = TVcalc2(B,M,N,mask)

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = abs(Dh) + abs(Dv);
%P = sqrt(P);
TV = norm(P(:).*mask,1);

end



% Total Variation: 0.5*||A*u(:)-b||_2^2 + lambda*TV(u)
function u = IRLS_TV(b,A,mu,M,N,tol,mask,minimask)

[u,~] = cgs2(A'*A,A'*b,1e-6,20);
%figure(109); imagesc(8.686*reshape(u,M,N)); colormap pink; caxis([0 1.2])

G(1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dx = kron(speye(N),D);

D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dy = kron(D,speye(M));

D = [Dx' Dy']';

ite_irls = 0;
error = 1;

while error > tol

    X = reshape(u,M,N);
    ite_irls = ite_irls + 1;
    Dh = diff(X,[],1);
    Dh = [Dh;zeros(1,N)];
    Dv = diff(X,[],2);
    Dv = [Dv zeros(M,1)];

    %Dx*X(:) - Dh(:);
    %Dy*X(:) - Dv(:);

    P = Dh.^2 + Dv.^2;
    eps = 0.01;
    P = 2*sqrt(P.^2 + eps^2);
    P = P.^(-0.5);
    P = P(:).*minimask;
    omega = speye(M*N);
    omega = spdiags(P,0,omega);
    W = kron(speye(2),omega);

    AtA = A'*A;
    %mu=5000;
    %[u] = cgs(AtA + mu*D'*W*D, A'*b,1e-6,200);
    [u] = cgs2( AtA + mu*D'*W*D, A'*b, 1e-6 , 20, u );

    G(ite_irls+1,1) = 1/2*(norm( (b - A*u) ))^2 + mu*TVcalc(u,M,N,minimask);
    error = abs(G(ite_irls+1) - G(ite_irls));

end

end


% A*x = b
function [x,ite_cgs] = cgs2(A,b,tol,maxit,varargin)

if length(varargin) == 1
    x = varargin{1};
else
    x= zeros(size(A,2),1);
end
r = A*x - b;
p = -r;
ite_cgs = 0;

while norm(r,2) > tol && ite_cgs < maxit

    alpha = (r'*r)/(p'*(A*p));
    x = x + alpha*p;
    rn = r + alpha*(A*p);
    beta = (rn'*rn)/(r'*r);
    p = -rn + beta*p;
    r = rn;
    ite_cgs = ite_cgs + 1;
end

end