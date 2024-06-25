clear,clc
close all

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Liver'];
targetDir = fullfile(baseDir,'sample');
refsDir = fullfile(baseDir,'ref');
resultsDir = fullfile(baseDir,'results','liverHetero');

% convertRfFilesToMat(targetDir)
% convertRfFilesToMat(refDir)

if (~exist(resultsDir,"dir")), mkdir(resultsDir); end
targetFiles = dir([targetDir,'\*.mat']);

%%
blocksize = 8;     % Block size in wavelengths
overlap_pc      = 0.8;
ratio_zx        = 12/8;
% ratio_zx        = 1;

% rect = [0.5, 3.3, 2.8, 3.7]; % Just liver, large ROI
% rect = [0.8, 4.2, 2.3, 1.7]; % Just liver, small ROI
% rect = [0.5, 3.3, 2.8, 2.5]; % Just liver, medium ROI
rect = [0.5, 1.8, 2.8, 4]; % liver & muscle
% rect = [];

% Bandwidth
fixedBW = true;
ratio = db2mag(-20);
freq_L = 1.4e6; freq_H = 5.3e6;

% Weight parameters
ratioCutOff = 10;
order = 5;
reject = 0.1;
extension = 3;

% SWTV
aSNR = 5; bSNR = 0.09;
desvMin = 15;

% reg FINAL VERSION
muBtv = 10^3.5; muCtv = 10^3.5;
muBswtv = 10^3; muCswtv = 10^2.5;
muBtvl1 = 10^3.5; muCtvl1 = 10^2;
muBwfr = 10^3.5; muCwfr = 10^2;

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.5];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%%
dataCols = zeros(length(targetFiles),16);
iAcq = 1;

fprintf("Patient no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(targetDir,targetFiles(iAcq).name));

dx = x(2)-x(1);
dz = z(2)-z(1);
sam1 = RF(:,:,1);

xFull = x*1e2; % [cm]
zFull = z*1e2; % [cm]

BmodeFull = db(hilbert(sam1));
BmodeFull = BmodeFull - max(BmodeFull(:));
%%

maskReady = false;

if ~maskReady
    figure('Units','centimeters', 'Position',[3 5 35 15]),
    tiledlayout(1,2)
    nexttile,
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 7])
    
    nexttile,
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 7])
    
    confirmation = '';
    while ~strcmp(confirmation,'Yes')
        h = drawfreehand('Multiclick',true);
        confirmation = questdlg('Sure?');
        if strcmp(confirmation,'Cancel')
            break
        elseif strcmp(confirmation,'No')
            delete(h)
        end
    end
    regionMask = createMask(h);
else    
    load(fullfile(baseDir,'results','homogeneous',[patient,'.mat']))
    figure('Units','centimeters', 'Position',[3 5 20 15]),
    imagesc(xFull,zFull,BmodeFull,dynRange); axis image; 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 min([3.5,max(zFull)])])
end


%%
maskLiver = regionMask;
save(fullfile(baseDir,'liverMask.mat'),'maskLiver')