function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE, SLOPE)
% function [ATT] = attenuation_phantoms_Np(FREQ, CHOICE, SLOPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION: Define attenuation of References Phantoms for the RPM Method
% INPUT:
%        - FREQ: Vector of frequencies in [MHz]
%        - CHOICE: Type of phantom
%                 - 1 :  phantom "low attenuation"
%                 - 2 :  phantom "high attenuation" 
%                 - 3 :  phantom 3 - "Agar" Referenced in Journal Rouyer TUFFC*
%                 - 4 :  phantom 4 - "Agar & Milk"
%                 - 5 :  high attenuation phantom using Madsen values - spline interpolation
%                 - 6 :  new low att phantom using Madsen values
%                 - 111 : linear attenuation dependence, add SLOPE (usually k-wave simulation)
%                 - 222 : power attenuation dependence, COILA EXVIVO (usually k-wave simulation)
%        - SLOPE: Slope for linear dependency CHOICE = 111 [dB/cm-MHz]
% OUTPUT: 
%     - ATT: Attenuation as function of frequency in [Np/cm]
%
% AUTHOR: Edmundo Arom Miranda based on LIM-Repository
% EXTRA: including "famous" phantom P4-cuello referenced in paper (CHOICE=3):
%
% [1] Rouyer, J., Cueva, T., Yamamoto, T., Portal, A., & Lavarello, R. J. (2016). 
% In vivo estimation of attenuation and backscatter coefficients from human thyroids. IEEE transactions on ultrasonics, 
% ferroelectrics, and frequency control, 63(9), 1253-1261.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch CHOICE
    case 1
        % phantom "low attenuation"
        ATT  = (0.0141*FREQ.^2 + 0.0838*FREQ - 0.0662)/8.6858;      % [Np/cm]
    case 2
        % phantom "high attenuation" 
        ATT = (0.0085*FREQ.^2 +0.5866*FREQ -0.3725)/8.6858;         % [Np/cm]
    case 3
        % phantom 3 - "Agar" Referenced in Journal Rouyer TUFFC*
        ATT =  (0.0076*FREQ.^2 + 0.1189*FREQ -0.0319)/8.6858;       % [Np/cm]
    case 4
        % phantom 4 - "Agar & Milk"
        ATT =  (0.0057*FREQ.^2 + 0.4432*FREQ -0.1000)/8.6858;       % [Np/cm]
    case 5 %high attenuation phantom using Madsen values - spline interpolation
        madsenHighAtt=[1.2000,2.7000,4.4600,6.3600,7.9800,10.2800];
        madsenFreqs=[2.5,5,7.5,10,12,15];
        ATT=interp1(madsenFreqs,madsenHighAtt./8.6858,FREQ,'spline');% [Np/cm]
    case 6 %new low att phantom using Madsen values
        madsenLowAtt=[0.076,0.327,0.729,1.296,2,3.345,5.411];
        madsenFreqs=[1,2.5,5,7.5,10,12,15];
        ATT=interp1(madsenFreqs,madsenLowAtt./8.6858,FREQ,'spline');% [Np/cm]
    
    case 111
        % phantom "linear dependency"
        ATT  = (0.0*FREQ.^2 + SLOPE*FREQ - 0.0)/8.6858;             % [Np/cm]
    case 222
        % phantom "power dependency" (COILA EXVIVO)
        a = 0.28; m = 1.34; % 6 BACKGROUND
        ATT  = (a*FREQ.^m)/8.6858;             % [Np/cm]

end