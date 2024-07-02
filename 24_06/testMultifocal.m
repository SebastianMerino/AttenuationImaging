clear,clc

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Attenuation\Simulation\24_06_27_layered'];
% refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Attenuation\Simulation\24_06_26_ref'];
resultsDir = fullfile(targetDir,'results');

[~,~] = mkdir(resultsDir);
targetFiles = dir([targetDir,'\rf*.mat']);
% refFiles = dir([refDir,'\rf*.mat']);
% tableName = 'simuInc.xlsx';

%%
iAcq = 3;
load(fullfile(targetDir,targetFiles(iAcq).name));
x = x*100; z = z*100;

Bmode1 = getBmode(rf(:,:,1));
Bmode2 = getBmode(rf(:,:,2));
Bmode3 = getBmode(rf(:,:,3));
dynRange = [-50 0];

%%
[~,Z] = meshgrid(x,z);
concatRf = rf(:,:,1).*(Z<1.5) + rf(:,:,2).*(Z>=1.5 & Z<2.5) + ...
    rf(:,:,3).*(Z>=2.5);
meanRf = mean(rf,3);

meanBm = getBmode(meanRf);
concatBm = getBmode(concatRf);

figure,tiledlayout(2,3)
nexttile,
imagesc(x,z,Bmode1, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,Bmode2, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,Bmode3, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,meanBm, dynRange)
colormap gray
axis image

nexttile,
imagesc(x,z,concatBm, dynRange)
colormap gray
axis image

%%
% figure,
% plot(z,mean(concatBm(:,46:50),2))
% hold on
% plot(z,mean(meanBm(:,46:50),2))
% grid on
% axis tight

figure,
plot(z,mean(concatBm,2))
hold on
plot(z,mean(meanBm,2))
grid on
axis tight
%% Lines
ix = 15;
figure('Units','centimeters', 'Position',[5 5 30 15]),
tiledlayout(3,1)
nexttile,
plot(z,rf(:,ix,1))
grid on
axis tight
xlim([z(1), z(end)])
ylim([-1,1]*1e4)

nexttile,
plot(z,rf(:,ix,2))
grid on
axis tight
xlim([z(1), z(end)])
ylim([-1,1]*1e4)

nexttile,
plot(z,rf(:,ix,3))
grid on
xlim([z(1), z(end)])
ylim([-1,1]*1e4)
%%
figure('Units','centimeters', 'Position',[5 5 30 15]),

plot(z,concatRf(:,ix))
% grid on
% axis tight
% 
% nexttile,
hold on
plot(z,meanRf(:,ix))
grid on
axis tight

%%
function bm = getBmode(rf)
    bm = db(hilbert(rf));
    bm = bm - max(bm(:));
end

