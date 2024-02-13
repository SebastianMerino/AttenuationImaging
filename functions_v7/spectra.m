%% Spectra
function [spect,psnr_Sp]=spectra(block,windowing,saran_layer,~,NFFT)
% Computes the average of the parallel echoes spectra in the ROI
%
% Inputs:
%   block           Data matrix, size nw x nx
%   windowing       Window vector, size nw x nx (the same for all cols)
%   saran_layer     Vector containing spectrum correction for the saran
%                   layer, of size NFFT. 0 is case there is no correction.
%   nw              Axial length of block (deprecated)
%   NFFT            Number of FFT points
%
% Outputs:
%   spect           Average of the parallel echoes spcetrin the ROI
%   psnr_Sp         Log spectrum
%       

block = block - mean(block);
block = block.*windowing;

spect = abs(fft(block,NFFT,1));     % Fourier transform proximal window
spect = spect.^2;                   % Sp is Intensity Now 

spect = mean(spect,2);   % Sp is the averaga of the parallel echoes in the ROI


% Saran-wrap correction factor for phantoms
if all(saran_layer)
    spect = spect./saran_layer;
end

% psnr_Sp = 10*log10(max(spect)./spect(end/2,:)); 
psnr_Sp = mean(log(spect),2);
end