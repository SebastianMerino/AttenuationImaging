function [] = rf2movie()
filename = 'phantom.rf';
%filename = 'whiteA.rf';
%clear;
close all; clc;

addpath(genpath(pwd))

[rf,header] = RPread(filename);

rf2env = @(x) abs(hilbert(x));   % x is a matrix
normrf = @(x) 20*log10(x) - 20*log10(max(max(x)));

[L1, L2, L3] = size(rf);

fs = header.sf;
x = linspace(0,0.038,L2)*100;   % [cm]
c = 1540;
z = (1:L1)*1/fs*c/2*100;   % [cm]
im1 = nan(size(rf));

Nfig = 1;

figure(Nfig);
hfig = gcf;
%title(['Phantom: ',material(1:8),' ',num2str(fmin),'MHz to ',num2str(fmax),'MHz'])
%set(gcf, 'Position', get(0, 'Screensize'));
set(hfig,'units','normalized','outerposition',[0 0 1 1])
font = 55;
set(gca,'fontsize',font)

for ii = 1:L3
    Im_db = normrf(rf2env(rf(:,:,ii)));
    imagesc(x,z,Im_db); 
    
    if ii == 1
        axis image;
        colormap gray;
                
        hb1 = colorbar;
        ylabel(hb1,'dB','FontSize', font)
        %title('B-mode image');
        xlabel('\bfLateral distance (cm)');
        ylabel('\bfAxial distance (cm)');
        set(gca,'FontSize',font);
    end
    caxis([-60 0]);
    title(ii)
    pause(0.0005);
end
end


