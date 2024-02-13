% ======================================================================
% Script to demonstrate effect of hyperchoic vs hypoechoic phantoms
% on borders.
% Created on Jan 6, 2024
% ======================================================================
clear, clc
%%
fs = 1.3320e+08;
nz = 556/2;
t = (0:nz-1)/fs;
fc1 = 5e6;
fc2 = 6e6;

t0 = 0.55e-6;
a2 = 1;
a1 = 1;
rf = a1*-cos(2*pi*fc1.*t).*(t<t0) + a2*cos(2*pi*fc2.*t).*(t>=t0) ;

NFFT = 1024;
f = (0:NFFT-1)/NFFT*fs;

spec = db(fft(rf,NFFT));
spec = spec - max(spec);

figure, tiledlayout(2,1);
nexttile, plot(t,rf)
xline(t0,'--')
axis tight
xlabel('t [s]')

nexttile, plot(f,spec)
xlim([0 10e6])
ylim([-30 0])
grid on
xlabel('Freq [Hz]')

%%
a1 = 1;
a2 = 0.1;
rf = -a1*cos(2*pi*fc1.*t).*(t<t0) + a2*cos(2*pi*fc2.*t).*(t>=t0) ;

spec = db(fft(rf,NFFT));
spec = spec - max(spec);

figure, tiledlayout(2,1);
nexttile, plot(t,rf)
xline(t0,'--')
axis tight
xlabel('t [s]')

nexttile, plot(f,spec)
xlim([0 10e6])
ylim([-30 0])
grid on
xlabel('Freq [Hz]')

% %%
% a2 = 1;
% a1 = 0.1;
% rf = -a1*cos(2*pi*fc1.*t).*(t<t0) + a2*cos(2*pi*fc2.*t).*(t>=t0) ;
% figure, plot(t,rf)
% 
% spec = db(fft(rf,NFFT));
% spec = spec - max(spec);
% figure, plot(f,spec)
% xlim([0 10e6])
% ylim([-30 0])
% grid on
% 
%%
nw = round(nz/2);
s = spectrogram(rf,tukeywin(nw,0.25),round(nw*0.5),NFFT,'twosided',fs);
spec = mean(db(s),2);
% spec = db(mean(abs(s),2));
spec = spec - max(spec);
figure, tiledlayout(2,1);
nexttile, plot(t,rf)
xline(t0,'--')
axis tight
xlabel('t [s]')

nexttile, plot(f,spec)
xlim([0 10e6])
ylim([-30 0])
grid on
xlabel('Freq [Hz]')

save_all_figures_to_directory('./24_02','border')