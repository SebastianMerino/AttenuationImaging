setup

dataDir = 'C:\Users\sebas\Documents\Data\Attenuation\Simulation\liver_RED';
refDir = fullfile(dataDir,'ref');
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = fullfile(dataDir,'test');
targetFiles = dir(fullfile(dataDir,"rf_qus_livernew202410_AC_test.mat"));

%%
blockParams.xInf = -2;
blockParams.xSup = 2;
blockParams.zInf = 1.7;
blockParams.zSup = 5;
blockParams.blocksize = [15 15];
blockParams.freqL = 3.5e6;
blockParams.freqH = 6.5e6;
blockParams.overlap = 0.8;

NptodB = log10(exp(1))*20;
alpha0Ref = 0.6; gammaRef = 1;

tol = 1e-3;

muBtv = 10^2.5; muCtv = 10^2.5;

iAcq = 1;
% for iAcq = 1:length(targetFiles)
%% Loading sample
fprintf("Simulation no. %i, %s\n",iAcq,targetFiles(iAcq).name);
out = load(fullfile(dataDir,targetFiles(iAcq).name));
xBm = out.x*1e2; % [cm]
zBm = out.z'*1e2; % [cm]
c0 = 1540;
sam1 = out.rf;
fs = out.fs;
kgrid = out.kgrid;
attenuation_map = out.attenuation_map;

% Plot region of interest B-mode image
bMode = db(hilbert(sam1));
bMode = bMode - max(bMode(:));

% SLD
[Sp,Sd,x_ACS,z_ACS,f] = getSld(sam1,xBm,zBm,fs,c0,blockParams);
sld = log(Sp) - log(Sd);
[m,n,p] = size(sld);

L = (z_ACS(2) - z_ACS(1))/(1 - blockParams.overlap)/2;   % (cm)
%% Generating Diffraction compensation
% Generating references
att_ref = alpha0Ref*f.^gammaRef/NptodB;
att_ref_map = zeros(m,n,p);
for jj=1:n
    for ii=1:m
        att_ref_map(ii,jj,:) = att_ref;
    end
end

for ff = 1:length(refFiles)
    out = load(fullfile(refDir,refFiles(ff).name));
    rfRef(:,:,ff) = out.rf(:,:,end);
end

sldRef = getSld(rfRef,xBm,zBm,fs,c0,blockParams);
compensation = sldRef - 4*L*att_ref_map;

%% RSLD-TV
b = sld - (compensation);
A1 = kron( 4*L*f , speye(m*n) );
A2 = kron( ones(size(f)) , speye(m*n) );
mask = ones(m,n,p);

tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muBtv,muCtv,m,n,tol,mask(:));
toc
BRSLD = (reshape(Bn*NptodB,m,n));
CRSLD = (reshape(Cn*NptodB,m,n));


%% Figures
% Creating masks and ideal map
[Xq,Zq] = meshgrid(xBm,zBm);
attIdeal = interp2(kgrid.y*100, kgrid.x*100 - kgrid.x(1)*100, attenuation_map, Xq, Zq);
roi = ones(size(bMode));

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.8];
bsRange = [-15,15];

figure('Units','centimeters', 'Position',[5 5 24 8]);
tiledlayout(1,4, "Padding","tight", 'TileSpacing','compact');

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
axis image
c = colorbar(t1, 'westoutside');
c.Label.String = 'dB';
title('B-mode')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

t2 = nexttile;
imagesc(xBm,zBm,attIdeal,attRange)
xlabel('Lateral [cm]'), % ylabel('Axial [cm]')
% axis equal
% xlim([x_ACS(1) x_ACS(end)]),
% ylim([z_ACS(1) z_ACS(end)]),
axis image
title('Ideal')

t3 = nexttile;
[~,hB,hColor] = imOverlayInterp(bMode,BRSLD,dynRange,attRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

t4 = nexttile;
[~,hB,hColor] = imOverlayInterp(bMode,CRSLD,dynRange,bsRange,1,...
    x_ACS,z_ACS,roi,xBm,zBm);
title('RSLD')
% colorbar off
xlabel('Lateral [cm]')

colormap(t1,gray)
colormap(t2,turbo)
colormap(t3,turbo)
colormap(t4,parula)

%%
save_all_figures_to_directory(resultsDir,char("sam"+iAcq+"fig"))
close all

%end


