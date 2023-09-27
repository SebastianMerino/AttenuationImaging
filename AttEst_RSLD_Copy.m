%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pontifica Universidad Católica del Perú
% Attenuation Estimator RSLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code starts here
function [] = AttEst_RSLD_Copy()
%% File management
clear; close all; clc;

addpath('./functions_att');

set(0,'DefaultTextFontName','Arial');
set(0,'DefaultTextFontSize',12);
font = 15;

blocksize = 20; % Block size in wavelengths

%% Loading data
load('parameters.mat')
dyn_range = dynamic_range_bmode;

load sam1.mat rf fs c0 x z
sam1 = rf;

% Image limits
ind_x =  overall_limits(1) <= x & x <= overall_limits(2); % 8th
x = x(ind_x);
ind_z =  overall_limits(3) <= z & z <= overall_limits(4); % 8th
z = z(ind_z);
sam1 = sam1(ind_z,ind_x);

% Plotting Bmode
dx = (x(end)-x(1))/(length(x)-1);
dz = (z(end)-z(1))/(length(z)-1);
x = x*1e2; % [cm]
z = z*1e2; % [cm]
Im = abs(hilbert(sam1));
Im_db = 20*log10(Im/max(Im(:)));

figure(1);
set(gca,'FontSize',font);
imagesc(x,z,Im_db); axis image; colormap gray; clim([0-dyn_range 0]);
hb1 = colorbar;
ylabel(hb1,'dB','FontSize', font)
%title('B-mode image');
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',font);

%% Loading parameters
winsize         = 0.5;
window_type     = 5;
saran_layer     = 0;

% parameters
mu = mu_values;

% Region for attenuation imaging
% x_inf = lateral_limits(1);
% x_sup = lateral_limits(2);
% z_inf = axial_limits(1);
% z_sup = axial_limits(2);
x_inf = 2.3; x_sup = 6.2;
z_inf = 4.7; z_sup = 10;

% Homogeneous regions for stats
c1x = center_inclusion(1);
c1y = center_inclusion(2);
% c2x = center_background(1);
% c2y = center_background(2);
c2x = c1x; c2y = 5.5;

% Other parameters
freq_L = frequency_limits(1);
freq_H = frequency_limits(2);
overlap_pc = overlap_percentage;

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
disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);

% Number of repetitions of datablocks within a single datablock
rpt = 1/(1-overlap_pc);   % r = rep = 1/(1-overlap)
rpt = round(rpt);
overlap = 1 - (1/rpt);

% Lateral samples between windows
wx  = round((blocksize*wl)/(dx*rpt));

% Number of lines laterally = repetitions * block without repetition
nx  = rpt*wx;
new_blocksize = round(nx*dx/(wl));

% RF data columns
L2   = size(sam1,2);

% Blocksize colums
ncol = floor((L2-(rpt-1)*wx)/wx);
sam1 = sam1(:,1:wx*(ncol+rpt-1));

% Actual rf data columns
L2 = size(sam1,2);
x  = x(1:L2);

xi = 1;
xf = L2;
x0 = (xi:wx:xf+1-nx);
x_ACS = x(1,x0+round(nx/2));

n  = length(x0);

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
L3 = length(rang);
p = L3;

L1   = size(sam1,1);   % RF data: rows
nrow = floor((L1-(rpt-1)*wz)/wz);        % Blocksize rows
sam1 = sam1(1:wz*(nrow+rpt-1),:);
L1   = size(sam1,1);   % RF data: rows
z    = z(1:L1);

zi = 1;
zf = L1;
z0 = (zi:wz:zf+1-nz);
m  = length(z0);

z_ACS = z(z0+round(nz/2));

z0p = z0 + (nw-1)/2;
z0d = z0 + (nz-1) - (nw-1)/2;

% Plot region of interest B-mode image
Im=abs(hilbert(sam1));   % envelope calculation
Im_db=20*log10(Im/max(Im(:)));   % log scale
figure(2); 
set(gca,'FontSize',font);
imagesc(x,z,Im_db); axis image; colormap gray; clim([0-dyn_range 0]);
hb2=colorbar; ylabel(hb2,'dB','FontSize', font)
xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
set(gca,'FontSize',font);

%% Calculating spectra
disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);


switch window_type
    case 5
        windowing = tukeywin(nw,0.25);   % Tukey Window. Parameter 0.25
    case 6
        windowing = hamming(nw);   % Hamming
    case 7
        windowing = rectwin(nw);   % Boxcar

end

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

        sub_block_p = sam1(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = sam1(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);

        [tempSp,~] = spectra(sub_block_p,windowing,saran_layer,nw,NFFT);
        [tempSd,~] = spectra(sub_block_d,windowing,saran_layer,nw,NFFT);
        Sp(ii,jj,:) = (tempSp(rang));
        Sd(ii,jj,:) = (tempSd(rang));
        
        % Sp(ii,jj,:) = spectra2(sub_block_p,Nw,Noverlap,f*1e6,fs);
        % Sd(ii,jj,:) = spectra2(sub_block_d,Nw,Noverlap,f*1e6,fs);
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
axis([0 10 -40 0]);

%% Diffraction compensation
load sam1.mat transducer
for ii=1:m
    zp = z0p(ii);
    zd = z0d(ii);
    zap = z(zp)*1e-2;
    zad = z(zd)*1e-2;
    D_p(ii,1,:) = diffraction_JR(transducer.fn,transducer.fl,band(rang),c0,zap);
    D_d(ii,1,:) = diffraction_JR(transducer.fn,transducer.fl,band(rang),c0,zad);
end

for ii=1:n
    D_p(:,ii,:) = D_p(:,1,:);
    D_d(:,ii,:) = D_d(:,1,:);
end

diffraction_compensation = log(D_p) - log(D_d);

%% Au = b
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


ROI = square_center(x_ACS,z_ACS,c1x,c1y,square_L) | ...
    square_center(x_ACS,z_ACS,c2x,c2y,square_L);
ROI = ones(size(BS));
[~,~,hAtt] = imoverlay2(Im_db,BS,[-50 -10],[0 1],1,x_ACS,z_ACS,ROI,x,z);
ylabel(hAtt,'Attenuation')

disp(['Frequency range: ',num2str(freq_L*1e-6,'%3.1f'),' - ',num2str(freq_H*1e-6,'%3.1f'),' MHz. c: ',...
    num2str(c0,'%4.1f'),' m/s. Wavelength: ',num2str(wl*1e3,'%2.2f'),' mm.']);
disp(['Blocksize. x: ',num2str(nx*dx*1e3,'%4.2f'),'mm, z: ',num2str(nz*dz*1e3,'%4.2f'),'mm, overlap: ',num2str(overlap*1e2,'%4.0f'),'%']);
disp(['Blocksize in wavelengths: ',num2str(new_blocksize,'%3.1f')]);
disp(['Blocksize in pixels. nf: ',num2str(p,'%i'),' nx: ',num2str(nx,'%i'),', nz: ',num2str(nz,'%i'),', nw: ',num2str(nw,'%i')]);
disp(['Region of interest. columns: ',num2str(ncol,'%i'),', rows: ',num2str(nrow,'%i')]);

%% Multi-frequency. Julien denoising
% bj = b. For Julien
jr.tol = 1e-3;
mask = ones(m,n,p);
jr.minimask = ones(m,n);

jr.factor = [0.1 0.3 1];
for qq = 1:length(jr.factor)

    for kk = 1:p
        jr.bNoise = b(:,:,kk);   % a single frequency
        jr.mean(kk,1) = abs(nanmean(jr.bNoise(:)));
        jr.sd(kk,1) = nanstd(jr.bNoise(:));
        jr.mu(kk,qq) = jr.factor(qq)*jr.sd(kk,1)/jr.mean(kk,1);
        jr.bDenoise = IRLS_ANIS_TV(jr.bNoise(:),speye(length(jr.bNoise(:))),jr.mu(kk,qq),m,n,jr.tol,mask,jr.minimask(:));
        jr.b(:,:,kk) = reshape(jr.bDenoise,m,n);

    end

    [jr.u{qq},~] = cgs2(A'*A,A'*jr.b(:),1e-6,20);

    % Standard SLD
    % BS: Beta. Attenuation coefficient slopes of blocks.
    % CS: Constants of blocks.
    BJtemp = jr.u{qq}(1:end/2); 
    CStemp = jr.u{qq}(end/2+1:end);
    BJtemp = 8.686*BJtemp;   % [dB.cm^{-1}.MHz^{-1}]
    BJ(:,:,qq) = reshape(BJtemp,m,n);
    CJ(:,:,qq) = reshape(CStemp,m,n);

    figure('Position',[100 100 600 300]), 
    tiledlayout(1,2)
    t1 = nexttile; 
    imagesc(x_ACS,z_ACS,BJ(:,:,qq), [0 1])
    colormap(t1,turbo)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Attenuation';
    title(['Julien denoising, Factor=',num2str(jr.factor(qq),2)])


    t2 = nexttile;
    imagesc(x_ACS,z_ACS,CJ(:,:,qq))
    colormap(t2,parula)
    axis equal tight
    c = colorbar;
    c.Label.String = 'Backscatter';
    title(['Julien denoising, Factor=',num2str(jr.factor(qq),2)])
end


%% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
%mu = 1e3*[0.5623,1.7783];
mu = 1e3*[0.4,1.2,3.6];
% mu = 3.3e1*[0.5623,1.7783];
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = reshape(Cn*8.686,m,n);
end


%% Plotting
for ii = 1:size(BR,3)
    figure('Position',[100 100 600 300]), 
    tiledlayout(1,2)
    t1 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), [0 1])
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

end


%% ------------------------ UTILITY FUNCTIONS ------------------------
function mask = square_center(x,y,c1x,c1y,L)
    [X,Y] = meshgrid(x,y);
    mask = (abs(X - c1x)<(L/2)) & (abs(Y - c1y)<(L/2));
end

function D = diffraction_JR(fn,fl,band,c0,za)

% Calcula el factor de compensacion por difraccion en la profundidad
% de analisis za a la frecuencia fa (X.Chen)
%
% Entradas:     fn: #focal del transductor
%               fl: Longitud focal del tx (m)
%               band: Frecuencia de analisis (Hz)
%               c: Velocidad del sonido en el medio (m/s)
%               za: Profundidad de analisis (m)
%
% Salidas:      D: Factor de compensacion por difraccion

a  = (fl/(fn*2));                        % Radio del transductor.
Gp = ( (pi*band/c0)*(a^2) ) / fl;           % Preasure gain factor.

for i = 1 : length(Gp)

    if 1/(1+pi/Gp(i))<=za/fl && za/fl<= 1/(1-pi/Gp(i))
        D(i) = ((pi*a^2)/(za^2))* 0.46*exp(-(0.46/pi)*Gp(i).^2*((fl/za)-1).^2);
    else
        %D(i) = ((pi*a^2)/(za.^2))* 1.07*(Gp(i)*((fl/za)-1)).^-2;
        D(i) = ((pi*a^2)/(za.^2))* 1.00*(Gp(i)*((fl/za)-1)).^-2;
    end
    % To PLOT
    %D(i)= D(i)/((pi*a^2)/(za^2));

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
