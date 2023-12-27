% ====================================================================== %
% Script to store and show ACS measurements in clinical data. 
% Created on Dec 05, 2023
% ====================================================================== %

%% Nodulos adenomatosos
TV = [0.44;0.55;0.79;0.87];
SWTV = [0.44;0.55;0.79;0.87];
TVL1 = [0.47;0.66;0.75;0.82];
WFR = [0.40;0.64;0.75;0.93];

boxplot([TV,SWTV,TVL1,WFR],{'TV','SWTV','TVL1','WFR'})
ylim([0.2 1.6])
xlabel('Method')
ylabel('ACS [dB/cm/MHz]')
grid on

%% Carcinomas
TV = [0.93;1.11;0.78;0.65];
SWTV = [0.93;1.11;0.78;0.65];
TVL1 = [0.90;1.07;0.81;0.81];
WFR = [0.89;0.98;0.78;0.76];
boxplot([TV,SWTV,TVL1,WFR],{'TV','SWTV','TVL1','WFR'})
ylim([0.2 1.6])
xlabel('Method')
ylabel('ACS [dB/cm/MHz]')
grid on

%% Nodulos coloides
TV = [-0.65;0.39;0.17];
SWTV = [-0.65;0.39;0.17];
TVL1 = [0.06;0.46;0.28];
WFR = [0.13;0.47;0.34];
boxplot([TV,SWTV,TVL1,WFR],{'TV','SWTV','TVL1','WFR'})
ylim([-0.1 1.4])
xlabel('Method')
ylabel('ACS [dB/cm/MHz]')
grid on
