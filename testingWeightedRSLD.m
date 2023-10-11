clear,clc
close all
addpath('./functions');
%%
% targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID316V2\06-08-2023-Generic'];
targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID316V2\06-08-2023-Generic'];

% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%     '\ID544V2\06-08-2023-Generic'];
refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% FOR LOOPING
% for iAcq = 1:8
iAcq = 4;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);

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

% jj = round(n/2);
% ii = round(m/4);
% figure('Units','centimeters', 'Position',[5 5 10 10])
% plot(f,10*log10(squeeze(Sp(ii,jj,:)/max(Sd(ii,jj,:)))),'k');
% hold on
% plot(f,10*log10(squeeze(Sd(ii,jj,:)/max(Sd(ii,jj,:)))),'r');
% hold off
% title('SLD at interface');
% xlabel('\bfFrequency (MHz)'); ylabel('\bfIntensity Norm. (dB)');
% axis([f(1) f(end) -30 10]);

%% Au = b

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

%% Standard SLD
[u,~] = cgs(A'*A,A'*b(:));

% BS: Beta. Attenuation coefficient slopes of blocks.
% CS: Constants of blocks.
BS = u(1:end/2); %CS = u(end/2+1:end);
BS = reshape(BS*8.686,m,n);    % [dB.cm^{-1}.MHz^{-1}]

% figure('Units','centimeters', 'Position',[5 5 20 8]);
% tiledlayout(1,2);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% t2 = nexttile; 
% imagesc(x_ACS,z_ACS,BS, attRange)
% colormap(t2,turbo)
% axis equal tight
% title(['RSLD, \mu=',num2str(mu(ii),2)])
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';


%% SOLVING WITH CVX
% b = (log(Sp) - log(Sd)) - (diffraction_compensation);
% 
% A1 = kron( 4*L*f , speye(m*n) );
% A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];
% 
% b = b(:);
% 
% cvx_begin
%     variable u(2*m*n)
%     minimize( norm( A * u - b, 2 ) )
% cvx_end
% 
% BS = u(1:end/2); %CS = u(end/2+1:end);
% BS = reshape(BS*8.686,m,n);    % [dB.cm^{-1}.MHz^{-1}]
% 
% figure('Units','centimeters', 'Position',[5 5 20 8]);
% tiledlayout(1,2);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% t2 = nexttile; 
% imagesc(x_ACS,z_ACS,BS, attRange)
% colormap(t2,turbo)
% axis equal tight
% title(['RSLD, \mu=',num2str(mu(ii),2)])
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';


%% Regularized SLD
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

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

for ii = 1:size(BR,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
    colormap(t2,turbo)
    axis equal tight
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

%saveas(gcf,[figDir,'\T',num2str(iAcq),'.png'])

%% SOLVING WITH CVX
% b = (log(Sp) - log(Sd)) - (diffraction_compensation);
% 
% A1 = kron( 4*L*f , speye(m*n) );
% A2 = kron( ones(size(f)) , speye(m*n) );
% A = [A1 A2];
% 
% b = b(:);
% 
% cvx_begin
%     variable u(2*m*n)
%     minimize( 0.5*norm( A * u - b, 2 )^2 )
% cvx_end


%% ----------------- WEIGHTS FROM LOCAL SNR -----------------
% FUNCTION
% SNR = 0:0.01:6;
% SNRopt = sqrt(1/(4/pi - 1));
% desvSNR = (SNR/SNRopt);
% legends = {};
% figure('Units','centimeters', 'Position',[5 5 20 10])
% for gamma = 1:10
%     w = desvSNR.^gamma.*exp(1-desvSNR.^gamma);
%     plot(SNR,w)
%     hold on
%     legends{gamma} = ['\gamma = ',num2str(gamma)];
% end
% xline(1.91, 'k--')
% xlabel('SNR (\mu/\sigma)')
% ylabel('Weights')
% legend([legends,'SNR_{opt}'])
% hold off

%% Getting local weights
%----------
gamma = 6;
%-------------

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
w = desvSNR.^gamma.*exp(1-desvSNR.^gamma);

t3 = nexttile;
imagesc(x_ACS,z_ACS,w, [0 1])
colormap(t3,parula)
c = colorbar;
axis image
title(['Weights, order=',num2str(gamma)])

saveas(gcf,[figDir,'\weightsT',num2str(iAcq),'.png'])

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
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

for ii = 1:size(BR,3)
    t2 = nexttile; 
    imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
    colormap(t2,turbo)
    axis equal tight
    title(['RSLD, \mu=',num2str(mu(ii),2)])
end
c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

saveas(gcf,[figDir,'\weightedT',num2str(iAcq),'.png'])

%%
