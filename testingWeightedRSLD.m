clear,clc
close all
addpath('./functions_v7');
%%
targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID316V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% FOR LOOPING
% Attenuation values from T1 to T8 and background, in that order
groundTruthTargets = [0.52,0.55,0.74,0.81,0.75,0.97,0.95,0.95,0.55];
c1x = 1.9; c1z = 1.93; roiL = 0.9; roiD = 0.6;
%%
for iAcq = 1:8
%iAcq = 7;
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

%% Au = b

b = (log(Sp) - log(Sd)) - (diffraction_compensation);

A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
A = [A1 A2];

%% Standard SLD
% [u,~] = cgs(A'*A,A'*b(:));
% 
% % BS: Beta. Attenuation coefficient slopes of blocks.
% % CS: Constants of blocks.
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
% title('SLD')
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';


%% Regularized SLD
% Regularization: Au = b
tol = 1e-3;

clear mask
mask = ones(m,n,p);
% mu = 1e4*[0.4,1.2,3.6];
%mu = logspace(2,4,5);
%mu = 1.2e4;
mu = 3.2e3;

BR = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1/10;
    [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),mu1,mu2,m,n,tol,mask(:));

    BR(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end

% Getting metrics
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BR(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
RSLD(iAcq) = r;


% Plotting
% figure('Units','centimeters', 'Position',[5 5 30 8]);
% tiledlayout(1,size(BR,3)+1);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% for ii = 1:size(BR,3)
%     t2 = nexttile; 
%     imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
%     colormap(t2,turbo)
%     axis equal tight
%     title(['RSLD, \mu=',num2str(mu(ii),2)])
% end
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';

%saveas(gcf,[figDir,'\T',num2str(iAcq),'.png'])



%  ======================================================
%  WEIGHTS FROM LOCAL SNR 
% ======================================================
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
colorbar;
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
BRW = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1/10;
    [Bn,Cn] = AlterOpti_ADMM(A1w,A2w,bw,mu1,mu2,m,n,tol,mask(:));

    BRW(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


% % Plotting
% figure('Units','centimeters', 'Position',[5 5 30 8]);
% tiledlayout(1,4);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% for ii = 1:size(BR,3)
%     t2 = nexttile; 
%     imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
%     colormap(t2,turbo)
%     axis equal tight
%     title(['RSLD, \mu=',num2str(mu(ii),2)])
% end
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';


% Getting metrics
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRW(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
WRSLD1(iAcq) = r;


%% NEW WEIGHTS

envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp-(nw-1)/2:zp+(nw-1)/2,xw:xw+nx-1);
        sub_block_d = envelope(zd-(nw-1)/2:zd+(nw-1)/2,xw:xw+nx-1);
        
        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
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
colorbar;
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
BRW2 = zeros(m,n,length(mu));
CR = zeros(m,n,length(mu));
for mm = 1:length(mu)
    mu1 = mu(mm);
    mu2 = mu1/10;
    [Bn,Cn] = AlterOpti_ADMM(A1w,A2w,bw,mu1,mu2,m,n,tol,mask(:));

    BRW2(:,:,mm) = (reshape(Bn*8.686,m,n));
    CR(:,:,mm) = (reshape(Cn,m,n));
end


% % Plotting
% figure('Units','centimeters', 'Position',[5 5 30 8]);
% tiledlayout(1,4);
% t1 = nexttile;
% imagesc(x,z,Bmode,dynRange)
% axis image
% colormap(t1,gray)
% colorbar(t1,'westoutside')
% title('Bmode')
% 
% for ii = 1:size(BR,3)
%     t2 = nexttile; 
%     imagesc(x_ACS,z_ACS,BR(:,:,ii), attRange)
%     colormap(t2,turbo)
%     axis equal tight
%     title(['RSLD, \mu=',num2str(mu(ii),2)])
% end
% c = colorbar;
% c.Label.String = 'Att. [db/cm/MHz]';


% Getting metrics
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD);

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
AttInterp = interp2(X,Z,BRW2(:,:,1),Xq,Zq);

r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(inc) - groundTruthTargets(9)) /...
    groundTruthTargets(9),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
WRSLD2(iAcq) = r;

%% Plotting three results
figure('Units','centimeters', 'Position',[5 5 30 8]);
tiledlayout(1,4);
t1 = nexttile;
imagesc(x,z,Bmode,dynRange)
axis image
colormap(t1,gray)
colorbar(t1,'westoutside')
title('Bmode')

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BR(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['RSLD, \mu=',num2str(mu(1),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRW(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['Weighted RSLD, \mu=',num2str(mu(1),2)])

t2 = nexttile; 
imagesc(x_ACS,z_ACS,BRW2(:,:,1), attRange)
colormap(t2,turbo)
axis equal tight
title(['Weighted RSLD v2, \mu=',num2str(mu(1),2)])

c = colorbar;
c.Label.String = 'Att. [db/cm/MHz]';

targetDir = fullfile(figDir,['T',num2str(iAcq)]);
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);
close all

end

%% Results
resultsRSLD = struct2table(RSLD);
resultsWRSLD1 = struct2table(WRSLD1);
resultsWRSLD2 = struct2table(WRSLD2);

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(resultsRSLD.meanInc, resultsRSLD.stdInc,'o', 'LineWidth',2)
hold on,
errorbar(resultsWRSLD1.meanInc, resultsWRSLD1.stdInc, 'x', 'LineWidth',2)
errorbar(resultsWRSLD2.meanInc, resultsWRSLD2.stdInc, 's', 'LineWidth',2)
plot(groundTruthTargets(1:8), "_k"	, 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 1.5])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'})
title('Results inclusion')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

figure('Units','centimeters', 'Position',[5 5 12 8])
errorbar(resultsRSLD.meanBack, resultsRSLD.stdBack,'o', 'LineWidth',2)
hold on,
errorbar(resultsWRSLD1.meanBack, resultsWRSLD1.stdBack, 'x', 'LineWidth',2)
errorbar(resultsWRSLD2.meanBack, resultsWRSLD2.stdBack, 's', 'LineWidth',2)
yline(groundTruthTargets(9), 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 1.5])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'})
title('Results background')
xlabel('Target')
ylabel('Attenuation [db/cm/MHz]')

%% CV
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.stdInc./resultsRSLD.meanInc*100,'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.stdInc./resultsWRSLD1.meanInc*100, 'x', 'LineWidth',2)
plot(resultsWRSLD2.stdInc./resultsWRSLD2.meanInc*100, 's', 'LineWidth',2)
plot(groundTruthTargets(1:8), "_k"	, 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 60])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'})
title('Results inclusion')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.stdBack./resultsRSLD.meanBack*100,'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.stdBack./resultsWRSLD1.meanBack*100, 'x', 'LineWidth',2)
plot(resultsWRSLD2.stdBack./resultsWRSLD2.meanBack*100, 's', 'LineWidth',2)
yline(groundTruthTargets(9), 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([0 60])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'})
title('Results background')
xlabel('Target')
ylabel('Coefficient of Variation [%]')

%% Mean Percentage error
figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.MPEInc, 'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.MPEInc, 'x', 'LineWidth',2)
plot(resultsWRSLD2.MPEInc, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([-30 100])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'}, 'Location','northwest')
title('Results inclusion')
xlabel('Target')
ylabel('MPE %')

figure('Units','centimeters', 'Position',[5 5 12 8])
plot(resultsRSLD.MPEBack, 'o', 'LineWidth',2)
hold on,
plot(resultsWRSLD1.MPEBack, 'x', 'LineWidth',2)
plot(resultsWRSLD2.MPEBack, 's', 'LineWidth',2)
yline(0, 'k--', 'LineWidth',2)
hold off
xlim([0.5 8.5])
ylim([-30 100])
grid on
legend({'RSLD','WRSLD v1','WRSLD v2'}, 'Location','northwest')
title('Results background')
xlabel('Target')
ylabel('MPE %')

%% ROI
figure('Units','centimeters', 'Position',[5 5 12 8])
[~,~,hC] = imoverlay2(Bmode,AttInterp,dynRange,attRange,0.5,x,z,back|inc,x,z);
xlabel('x [cm]')
ylabel('y [cm]')
hC.Label.String = 'Attenuation [db/cm/MHz]';

targetDir = fullfile(figDir,'Overall');
if(~exist(targetDir,"dir")), mkdir(targetDir); end
save_all_figures_to_directory(targetDir);


