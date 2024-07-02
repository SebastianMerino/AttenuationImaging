clear,clc

targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
    'Attenuation\simulations_processed\24_06_26_inc'];

targetFiles = dir([targetDir,'\rf*.mat']);
% refFiles = dir([refDir,'\rf*.mat']);
% tableName = 'simuInc.xlsx';

%%
iAcq = 2;
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

%% Attempt 1
transition = 20;
foc1 = 1./(1+exp(transition*(Z-1.5)));
foc2 = 1./(1+exp(-transition*(Z-1.5))) + 1./(1+exp(transition*(Z-2.5))) - 1;
foc3 = 1./(1+exp(-transition*(Z-2.5)));
figure, hold on
plot(z,foc1(:,1))
plot(z,foc2(:,1))
plot(z,foc3(:,1))
plot(z,foc1(:,1) + foc2(:,1) + foc3(:,1))
hold off
grid on 
axis tight

newRf = rf(:,:,1).*foc1 + rf(:,:,2).*foc2 + ...
    rf(:,:,3).*foc3;

%% Attempt 2
transition = 20;
foc1 = abs(Z-1).^-2./(abs(Z-1).^-2 + abs(Z-2).^-2 + abs(Z-3).^-2);
foc2 = abs(Z-2).^-2./(abs(Z-1).^-2 + abs(Z-2).^-2 + abs(Z-3).^-2);
foc3 = abs(Z-3).^-2./(abs(Z-1).^-2 + abs(Z-2).^-2 + abs(Z-3).^-2);
figure, hold on
plot(z,foc1(:,1))
plot(z,foc2(:,1))
plot(z,foc3(:,1))
plot(z,foc1(:,1) + foc2(:,1) + foc3(:,1))
hold off
grid on 
axis tight

newRf = rf(:,:,1).*foc1 + rf(:,:,2).*foc2 + ...
    rf(:,:,3).*foc3;
%%
meanBm = getBmode(meanRf);
concatBm = getBmode(concatRf);
newBm = getBmode(newRf);

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

nexttile,
imagesc(x,z,newBm, dynRange)
colormap gray
axis image


%%
figure,
plot(z,mean(concatBm(:,46:50),2))
hold on
plot(z,mean(meanBm(:,46:50),2))
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

