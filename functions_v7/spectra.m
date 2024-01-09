%% Spectra
function [spect,psnr_Sp]=spectra(block,windowing,saran_layer,nw,NFFT,T)

block = block - ones(nw,1)*mean(block);
block = block.*windowing;

% figure(10800); 
% subplot(211)
% plot(p(:,floor((d+1)/2)),'b'); hold on;

%p=p-ones(size(p,1),1)*mean(p);

% plot(p(:,floor((d+1)/2)),'r'); hold off; xlim([0 650]); ylim([-0.05 0.05]);

spect = abs(fft(block,NFFT,1));   % Fourier transform proximal window
spect = spect.^2;        % Sp is Intensity Now 

% GABRIEL IF_IA_QEA 11-12-2015
% Sp = if_ia_qea(Sp);
% save data_gabriel.mat Sp
% figure(8888);
% plot(1:4096,Sp(:,1));
% Sp2=if_ia_qea(Sp(:,1))
% hold on
% plot(1:4096,Sp2(:,1),'r');

%psnr_Sp = 10*log10(max(Sp)./Sp(end/2,:)); 
%psnr_Sp = min(psnr_Sp); 

spect = mean(spect,2);   % Sp is the averaga of the parallel echoes in the ROI

%psnr_Sp = 10*log10(ma|x(Sp)./Sp(end/2,:)); 

% Swaran-wrap correction factor for phantoms
if saran_layer == 1
%T = saran_wrap(band);
spect = spect./T;
%disp("Done")
end

%psnr_Sp= 10*log10(max(Sp(:))/Sp(end/2)); 

% subplot(212)
% plot(band((1:NFFT/2+1))*1e-6,10*log10(Sp((1:NFFT/2+1))/max(Sp((1:NFFT/2+1)))),'k');   %
% legend('Spectrum of Window'); 
% title('Spectrum vs Frequency'); xlabel('Frequency (MHz)'); ylabel('Spectrum');
% axis([0 125 -80 0]);% ylim([-80 0]);

% r = z(z);
% diff = 1;
% if diff == 1
% Ds = diffraction_V2(r,band);
% spect = spect./Ds';
% end

psnr_Sp = 10*log10(max(spect)./spect(end/2,:)); 
%psnr_Sp = min(psnr_Sp); 

end