% ====================================================================== %
% Script to create arbitrary masks for homogeneous ROIs in clinical data. 
% Created on March 25, 2024
% ====================================================================== %
setup,
warning('off'); %% Important to turn warns in loop off

baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
'Attenuation\TEST'];
resultsDir = fullfile(baseDir,'results');

resultFiles = dir(fullfile(resultsDir,'*_rect.mat'));
nFiles = length(resultFiles);
%%

varNames = {'acs','method','patient','sample'};
varTypes = {'double','string','int16','int16'};
T = table('Size',[nFiles*3,4], 'VariableTypes',varTypes,'VariableNames',varNames);


resultFiles = dir(fullfile(resultsDir,'*_results.mat'));
nFiles = length(resultFiles);

for ii = 0:nFiles-1
    fileName = resultFiles(ii+1).name;
    patient = str2double(fileName(1:3));
    sample = str2double(fileName(5:6));
    load(fullfile(resultsDir,fileName));
    T(ii*3 + 1,:) = {mean(BR,'all'),"RSLD",patient,sample};
    T(ii*3 + 2,:) = {mean(BSWTV,'all'),"SWTV",patient,sample};
    T(ii*3 + 3,:) = {mean(BSWIFT,'all'),"SWIFT",patient,sample};
end


%%

figure('Units','centimeters', 'Position', [5 5 20 10]),
boxchart(categorical(T.method, {'RSLD','SWTV','SWIFT'}), ...
    T.acs, 'MarkerStyle','o');
ylim([-0.4,1.6])
grid on
title('ACS')
% legend
