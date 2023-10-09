clear,clc
close all
addpath('./functions_att');

targetDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID316V2\06-08-2023-Generic'];
% targetDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

refDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Attenuation' ...
    '\ID544V2\06-08-2023-Generic'];
% refDir = ['C:\Users\smerino.C084288\Documents\MATLAB\Datasets\' ...
%     'Attenuation\ID544V2\06-08-2023-Generic'];

croppedDir = [targetDir,'\cropped'];
figDir = [targetDir,'\fig'];
if (~exist(figDir,"dir")), mkdir(figDir); end

%% Loading sample
iAcq = 4;
load([croppedDir,'\T',num2str(iAcq),'.mat'])
load([refDir,'\compensation.mat']);

%% SNR
envelope = abs(hilbert(sam1));

Nwz = floor(nw/3);
Nwx = floor(nx/6);
iz = 1:Nwz:size(sam1,1);
ix = 1:Nwx:size(sam1,2);
SNR = zeros(length(iz)-2,length(ix)-2);
for jj=1:length(ix)-2
    for ii=1:length(iz)-2
        xw = ix(jj) ;   % x window
        zw = iz(ii);    % z window
        subBlock = envelope(zw:zw+2*Nwz-1,xw:xw+2*Nwx-1);        
        SNR(ii,jj) = mean(subBlock(:))./std(subBlock(:));
    end
end

xSNR = x(ix(1:end-2)+Nwx);
zSNR = z(iz(1:end-2)+Nwz);

figure, tiledlayout(1,2)
nexttile,
imagesc(xSNR,zSNR,db(SNR))
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')

[Xq,Zq] = meshgrid(x,z);
[X,Z] = meshgrid(xSNR,zSNR);
SNRinterp = interp2(X,Z,SNR,Xq,Zq,"cubic",1);

nexttile,
imagesc(x,z,db(SNRinterp))
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')

%%
tiledlayout(1,2)
nexttile,
imagesc(x,z,db(SNRinterp))
c = colorbar;
ylabel(c,'dB')
axis image
title('SNR')

nexttile,
L = watershed(-db(SNRinterp));
imagesc(x,z,L)
c = colorbar;
axis image
title('SNR')


%%
figure('Units','centimeters', 'Position',[5 5 30 8]), 
tiledlayout(1,3)
t1 = nexttile;
imagesc(x,z,Bmode)
colormap(t1,gray)
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Bmode')

L = watershed(Bmode);
t2 = nexttile;
imagesc(x,z,L)
colormap(t1,parula)
axis image
xlabel('x [cm]'), ylabel('z [cm]')
title('Bmode')
