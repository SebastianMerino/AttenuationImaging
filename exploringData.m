clear,clc
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\', ...
    'Attenuation\DataQUS_4_Merino'];rfDir = [baseDir,'\Hashimoto'];
rfFiles = dir([rfDir,'\*.mat']);

refDir = [baseDir,'\References\P4-CUELLO-2'];
refFiles = dir([refDir,'\*.mat']);
%load([refDir,'\',refFiles(1).name])
%%
load([rfDir,'\',rfFiles(1).name]);
dynRange = [-60,-10];

imagesc(x,z,Bmode,dynRange)
colormap gray
colorbar
axis equal tight
