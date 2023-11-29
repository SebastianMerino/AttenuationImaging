%% PHANTOMSSS
clear,
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
%      '\ID316V2\06-08-2023-Generic'];

croppedDir = [baseDir,'\cropped'];
croppedFiles = dir([croppedDir,'\*.mat']);
NptodB = log10(exp(1))*20;

figDir = 'C:\Users\sebas\Pictures\ISBI2024\23-11';
mkdir(figDir)
%% For looping each phantom

for iAcq = 2:3
fprintf("Phantom no. %i, %s\n",iAcq,croppedFiles(iAcq+5).name);
load(fullfile(croppedDir,croppedFiles(iAcq+5).name));

%% Plotting B-mode
dynRange = [-50,0];
attRange = [0.4,1.1];
bsRange = [-2 2];

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x,z,Bmode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
axis image
colormap(gray)
title('Bmode')
c = colorbar;
c.Label.String = 'dB';
fontsize(gcf,8,'points')

%% Setup
% ROI
c1x = 1.95; c1z = 1.93;
roiL = 1; roiD = 0.6;
roiLz = 1.5;
%roiL = 1.2; roiD = 0.6;
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(x,z);
rI = 0.6; rB = 1.2; % Both
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
[back,inc] = getRegionMasks(x,z,c1x,c1z,roiL,roiD,roiLz);

%%
groundTruthTargets = [0.97,0.95,0.95,0.55];

% System of eq
b = (log(Sp) - log(Sd)) - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
tol = 1e-3;
mask = ones(m,n,p);
%% RSLD-TV

muB = 10.^(3:0.5:4);
muC = 10.^(0.5:0.5:2.5);
% switch iAcq
%     case 1
%         muB = 10^3; muC = 10^2.5;
%     case 2
%         muB = 10^3.5; muC = 10^2;
%     case 3
%         muB = 10^3; muC = 10^0.5;
% end
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        % tic
        [Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        % tic
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

fprintf("TV: muB = 10^%.1f, muC= 10^%.1f\n",log10(muBopt),log10(muCopt))
BRTV = BRopt;

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x_ACS,z_ACS,BRTV, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-TV')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off


%% British Columbia Approach

envelope = abs(hilbert(sam1));

SNR = zeros(m,n);
for jj=1:n
    for ii=1:m
        xw = x0(jj) ;   % x window
        zp = z0p(ii);
        zd = z0d(ii);

        sub_block_p = envelope(zp:zp+nz/2-1,xw:xw+nx-1);
        sub_block_d = envelope(zd:zd+nz/2-1,xw:xw+nx-1);

        temp = [sub_block_p(:) sub_block_d(:)];
        SNR(ii,jj) = mean(temp)/std(temp);
    end
end

SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
aSNR = 1; bSNR = 0.1;
desvMin = 15;
w = aSNR./(1 + exp(bSNR.*(desvSNR - desvMin)));


%muB = 10.^(1.5:0.5:4.5);
%muC = 10.^(1:0.5:3.5);
muB = 10.^(2:0.5:3.5);
muC = 10.^(1.5:0.5:2.5);
% switch iAcq
%     case 1
%         muB = 10^3; muC = 10^3;
%     case 2
%         muB = 10^3; muC = 10^2.5;
%     case 3
%         muB = 10^2.5; muC = 10^1;
% end
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        %tic
        [Bn,Cn] = AlterOptiAdmmAnisWeighted(A1,A2,b(:),muB(mmB),muC(mmC),...
        m,n,tol,mask(:),w);
        %toc
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        
        % figure('Units','centimeters', 'Position',[5 5 5 4]);
        % imagesc(x_ACS,z_ACS,BR, attRange)
        % xlabel('Lateral [cm]'), ylabel('Axial [cm]')
        % colormap(turbo)
        % axis image
        % title('RSLD-SWTV')
        % subtitle(['muB=',num2str(muB(mmB),2),', muC=',num2str(muC(mmC),2)])
        % c = colorbar;
        % c.Label.String = 'ACS [dB/cm/MHz]';
        % fontsize(gcf,8,'points')
        
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end
%%
% SWTV: muB = 10^2.5, muC= 10^2.0
fprintf("SWTV: muB = 10^%.1f, muC= 10^%.1f\n",log10(muBopt),log10(muCopt))
BRBC = BRopt;

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x_ACS,z_ACS,BRBC, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-SWTV')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

%% Minimizing BS log ratio
% muB = 10.^(3:0.5:4);
% muC = 10.^(0.5:0.5:2.5);
switch iAcq
    case 1
        muB = 10^3.5; muC = 10^2.5;
    case 2
        muB = 10^3; muC = 10^0.5;
    case 3
        muB = 10^3; muC = 10^0.5;
end
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        % tic
        [Bn,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(mmB),muC(mmC),m,n,tol,mask(:));
        % tic
        BR = reshape(Bn*NptodB,m,n);
        CR = reshape(Cn*NptodB,m,n);

        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end

fprintf("TVL1: muB = 10^%.1f, muC= 10^%.1f\n",log10(muBopt),log10(muCopt))
BRTVL1 = BRopt;

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x_ACS,z_ACS,BRTVL1, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-TVL1')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

%% Weighted approach
% Initial estimation
muB = 10^3; muC = 10^0.5;
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muB(1),muC(1),m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight function
ratioCutOff = 6;
order = 5;
reject = 0.1;
extension = 3;
w = (1-reject)*(1./((bscMap/ratioCutOff).^(2*order) + 1))+reject;
w = movmin(w,extension);

W = repmat(w,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

muB = 10.^(3:0.5:5);
muC = 10.^(0.5:0.5:2.5);
minRMSE = 100;
for mmB = 1:length(muB)
    for mmC = 1:length(muC)
        %tic
        [Bn,Cn] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muB(mmB),muC(mmC),m,n,tol,mask(:),w);
        %toc

        BR = reshape(Bn*NptodB,m,n);
        CR = (reshape(Cn*NptodB,m,n));
        AttInterp = interp2(X,Z,BR,Xq,Zq);
        RmseInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
            "omitnan") ;
        RmseBack = mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
            "omitnan");
        RMSE = sqrt((RmseInc + RmseBack)/2);
        if RMSE<minRMSE
            minRMSE = RMSE;
            muBopt = muB(mmB);
            muCopt = muC(mmC);
            BRopt = BR;
            CRopt = CR;
        end
    end
end
%%
fprintf("WFR: muB = 10^%.1f, muC= 10^%.1f\n",log10(muBopt),log10(muCopt))
BRWTik = BRopt;

figure('Units','centimeters', 'Position',[5 5 5 3.8]);
imagesc(x_ACS,z_ACS,BRWTik, attRange)
xlabel('Lateral [cm]'), ylabel('Axial [cm]')
colormap(turbo)
axis image
title('RSLD-WFR')
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
fontsize(gcf,8,'points')
hold on 
rectangle('Position',[x0mask z0mask roiL roiLz], 'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask-roiD-roiL/2 z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
rectangle('Position',[x0mask+roiL+roiD z0mask roiL/2 roiLz],...
    'LineStyle','--', 'LineWidth',1)
hold off

%%
groundTruthTargets = [0.97,0.95,0.95,0.55];
AttInterp = interp2(X,Z,BRTV,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsTV(iAcq) = r;


AttInterp = interp2(X,Z,BRBC,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsSWTV(iAcq) = r;

AttInterp = interp2(X,Z,BRTVL1,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsTVL1(iAcq) = r;

AttInterp = interp2(X,Z,BRWTik,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.MPEInc = mean( (AttInterp(inc) - groundTruthTargets(iAcq)) /...
    groundTruthTargets(iAcq),"omitnan") * 100;
r.MPEBack = mean( (AttInterp(back) - groundTruthTargets(end)) /...
    groundTruthTargets(end),"omitnan") * 100;
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
MetricsWFR(iAcq) = r;

end

%%
results1 = struct2table(MetricsTV);
results2 = struct2table(MetricsSWTV);
results3 = struct2table(MetricsTVL1);
results4 = struct2table(MetricsWFR);

disp('Bias Inc')
disp(results1.MPEInc)
disp(results2.MPEInc)
disp(results3.MPEInc)
disp(results4.MPEInc)

disp('Bias Back')
disp(results1.MPEBack)
disp(results2.MPEBack)
disp(results3.MPEBack)
disp(results4.MPEBack)

disp('CNR')
disp(results1.cnr)
disp(results2.cnr)
disp(results3.cnr)
disp(results4.cnr)

%%
save_all_figures_to_directory(figDir,'figure')
close all
