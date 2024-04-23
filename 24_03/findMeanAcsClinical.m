% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
clear,clc
close all

%% Clinical case
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Thyroid_Data_PUCP_UTD'];
refsDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\REFERENCES'];
resultsDir = fullfile(baseDir,'results','all_cases');
T = readtable('params_new.xlsx');
[~,~,~] = mkdir(resultsDir);

%%
muB0 = 1e3; muC0 = 10^0;
ratioCutOff     = 10;
order = 5;
reject = 0.3;
extension = 3; % 1 or 3

% muBtv = 10^4; muCtv = 10^1;
% muBwfr = 10^4; muCwfr = 10^0;
muBtv = 10^3.5; muCtv = 10^1;
muBwfr = 10^3.5; muCwfr = 10^0;

% Plotting constants
dynRange = [-50,0];
attRange = [0,1.5];
bsRange = [-20 20];
NptodB = log10(exp(1))*20;

%%
% colTv = [];
% adTv = [];
% colWfr = [];
% adWfr = [];
% clases = [];
T.acsTv = nan(height(T),1);
T.acsWfr = nan(height(T),1);
T.numBlocks = nan(height(T),1);
sldLines = zeros(70,height(T));
%% Loading case FULL VERSION
for iAcq = 1:height(T)
    %%
    patient = num2str(T.patient(iAcq));
    load(fullfile(resultsDir,[patient,'.mat']));
    %%
    % clase = T.clase{iAcq};
    % if clase =='C'
    %     colTv = [colTv; sum(BR.*regionMaskAcs,'all')./sum(regionMaskAcs,'all')];
    %     colWfr = [colWfr; sum(BRWFR.*regionMaskAcs,'all')./sum(regionMaskAcs,'all')];
    %     clases = [clases, {clase}];
    % elseif clase == 'A'
    %     adTv = [adTv; sum(BR.*regionMaskAcs,'all')./sum(regionMaskAcs,'all')];
    %     adWfr = [adWfr; sum(BRWFR.*regionMaskAcs,'all')./sum(regionMaskAcs,'all')];
    %     clases = [clases, {clase}];
    % end
    
    % dz = 0.6352 mm, dx = 0.4471 mm
    % area = 0.2840 mm2
    T.numBlocks(iAcq) =  sum(regionMaskAcs,'all');
    if     T.numBlocks(iAcq) >= 150
    T.acsTv(iAcq) = mean(BR(regionMaskAcs));
    % T.acsWfr(iAcq) =  sum(BRWFR.*w,'all')./sum(w,'all');
    T.acsWfr(iAcq) = mean(BRWFR(regionMaskAcs));
    end
    % sldLines(:,iAcq) = sldLine;
end
%%

% figure,
% boxchart(categorical(T.clase), T.numBlocks, 'MarkerStyle','o');

figure('Units','centimeters', 'Position', [5 5 20 10]),
tiledlayout(1,2)
nexttile,
boxchart(categorical(T.clase), T.acsTv, 'MarkerStyle','o');
ylim([-1,2])
grid on
title('TV')

nexttile,
boxchart(categorical(T.clase), T.acsWfr, 'MarkerStyle','o');
ylim([-1,2])
grid on
title('WFR')

%%
binEdges = -0.5:0.2:2;
figure('Units','centimeters', 'Position', [5 5 10 10]),
tiledlayout(2,1)
nexttile,
histogram(T.acsTv, 'BinEdges',binEdges)
title('TV')

nexttile,
histogram(T.acsWfr, 'BinEdges',binEdges)
title('WFR')
