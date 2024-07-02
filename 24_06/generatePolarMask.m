clear,clc
close all

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Liver_24_06_28\set1'];
targetDir = fullfile(baseDir,'sample');
refsDir = fullfile(baseDir,'ref');
resultsDir = fullfile(baseDir,'results','initial');

% convertRfFilesToMat(targetDir)
% convertRfFilesToMat(refDir)

if (~exist(resultsDir,"dir")), mkdir(resultsDir); end
targetFiles = dir([targetDir,'\*.mat']);

%%

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.5];
bsRange = [-15 15];
NptodB = log10(exp(1))*20;

%%
iAcq = 5;

fprintf("Patient no. %i, %s\n",iAcq,targetFiles(iAcq).name);
load(fullfile(targetDir,targetFiles(iAcq).name));

xFull = th; % [deg]
zFull = (r-r(1))*1e2; % [cm]

BmodeFull = db(hilbert(rf));
BmodeFull = BmodeFull - max(BmodeFull(:));
%%

maskReady = false;

if ~maskReady
    figure('Units','centimeters', 'Position',[3 5 30 15]),
    tiledlayout(1,2)
    nexttile,
    imagesc(xFull,zFull,BmodeFull,dynRange); 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 10])
    
    nexttile,
    imagesc(xFull,zFull,BmodeFull,dynRange); 
    colormap gray; clim(dynRange);
    hb2=colorbar; ylabel(hb2,'dB')
    xlabel('\bfLateral distance (cm)'); ylabel('\bfAxial distance (cm)');
    ylim([0.1 10])
    
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
    ylim([0.1 10])
end


%%
maskLiver = regionMask;
[~,~,~] = mkdir(fullfile(baseDir,'masks'));
save(fullfile(baseDir,'masks',targetFiles(iAcq).name),'maskLiver')