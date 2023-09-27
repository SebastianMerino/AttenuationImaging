%% Spectra
function [spect,psnr_Sp]=spectra(block,windowing,saran_layer,nw,NFFT)
% Returns the power spectrum on the block. 
% It calculates the spectrum of each column and averages all columns.
% 
% Inputs: 
%       block           matrix of RF signal in block
%       windowing       window function as a vector
%       saran_layer     idk
%       nw              number of points in the window
%       NFFT            number of points for FFFT
%  
% Outputs:
%       spect           column vector containing the power spectrum
%       psnr_Sp         ratio from the max intensity to the nyquist
%                       intensity in db.

block = block - ones(nw,1)*mean(block);
block = block.*windowing;

spect = abs(fft(block,NFFT,1));   % Fourier transform proximal window
spect = spect.^2;        % Sp is Intensity Now 
spect = mean(spect,2);   % Sp is the average of the parallel echoes in the ROI
%spect = spect./sum(spect); % NORMALIZATION

% Swaran-wrap correction factor for phantoms
if saran_layer == '1'
T = saran_wrap(band);
spect = spect./T;
end

psnr_Sp = 10*log10(max(spect)./spect(end/2,:)); 
%psnr_Sp = min(psnr_Sp); 

end