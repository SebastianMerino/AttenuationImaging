%% Spectra
function spect = spectra2(block,Nw,Noverlap,f,fs)

block = block - mean(block,1);
%windowing = tukeywin(Nw,0.25);
windowing = hamming(Nw);

spect = zeros(length(f),1);
for ix = 1:size(block,2)
    spectLine = spectrogram(block(:,ix),windowing,Noverlap,f,fs);
    spect = spect + mean(abs(spectLine).^2,2);        % Sp is Intensity Now 
end
spect = spect/size(block,2);

% amp normalization (DOES NOT WORK)
% spect = spect/sum(spect);
end