% ====================================================================== %
% Script to plot weight maps
% Created on Jan 23, 2024
% ====================================================================== %
clear,clc
close all
addpath('./functions_v7');
addpath('./AttUtils');
%%
bscMap = -25:0.1:25;

ratioCutOff     = 15;
order = 5;
reject = 0.3;
extension = 3; % 1 or 3

w = (1-reject)* (1./((bscMap/ratioCutOff).^(2*order) + 1)) + reject;
w = movmin(w,extension);

figure('Units','centimeters', 'Position',[5 5 10 5])
plot(bscMap,w),
grid on
xlabel('BSC change [dB/cm]')
    ylabel('Weights')
xlim([-25 25])
ylim([0 1])

%%
dBgain = 0.3;
w = db2mag(-dBgain*abs(bscMap));


figure('Units','centimeters', 'Position',[5 5 10 5])
plot(bscMap,db(w)),
grid on
xlabel('BSC change [dB/cm]')
ylabel('Weights [dB]')

xlim([-25 25])
ylim([-15 0])

