clear,clc

baseDir = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets\Attenuation\ID316V2\06-08-2023-Generic';
rawDir = [baseDir,'\06-08-2023-Generic'];

%%



%% B-mode image
for iFile = 1:length(rfFiles)
    load([rfDir,'\',rfFiles(iFile).name]);
    
    % Frequency response rf - 3 to 9 MHz
    %freqResponse = median(abs(fft(RF)),2);
    %f = (1:size(RF,1)) /size(RF,1) * fs;
    %plot(f/1e6,freqResponse)
    
    % Plotting Bmode
    dx = x(2)-x(1);
    dz = z(2)-z(1);
    x = x*1e2; % [cm]
    z = z*1e2; % [cm]
    dynRange = [-60,-10];
    
    figure, imagesc(x,z,Bmode,dynRange)
    colormap gray
    colorbar
    axis equal tight
end
