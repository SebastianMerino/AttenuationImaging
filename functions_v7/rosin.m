
% rosin 1 normal, rosin 2 para L curve
function [threshold,ind2]=rosin(X,Y)

%clear all; close all; clc; home;

% rosin 1
%[Y1, ind] = max(Y);
%Y2 = Y(end);
%X1 = X(ind);
%X2 = X(end);

% rosin 2
ind = 1;
Y1 = Y(ind);
Y2 = Y(end);
X1 = X(ind);
X2 = X(end);

% Line Equation
% y = a*x +b
% y -a*x - b = 0
a = (Y2-Y1)/(X2-X1);
b = Y1 - a*X1;
% Distance from (x0,y0)
% abs(y0 - a*x0 -b)/sqrt(1+a^2)

Y=Y(ind:end);
X=X(ind:end);

distance_rosin = abs(Y - a*X - b)/sqrt(1+a^2);
[~, ind2] = max(distance_rosin);
threshold = X(ind2);

end
