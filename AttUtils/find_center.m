%% Region type: One circular inclusion

% Find center
% input: coordinates of three points
% output: center

function [c_row,c_col,radius]=find_center(x,y)
% Line 1 - Perpendicular to the line (x1,y1),(x2,y2)
% a1*x+b1*y=c1
mp1=[(x(1)+x(2))/2 (y(1)+y(2))/2];   % mid-point
a1=x(2)-x(1);
b1=y(2)-y(1);
c1=a1*mp1(1)+b1*mp1(2);

% Line 2 - Perpendicular to the line (x2,y2),(x3,y3)
% a2*x+b2*y=c2
mp2=[(x(2)+x(3))/2 (y(2)+y(3))/2];   % mid-point
a2=x(3)-x(2);
b2=y(3)-y(2);
c2=a2*mp2(1)+b2*mp2(2);

% Ax=b
A= [a1 b1; a2 b2]; b= [c1; c2];
[center]=cgs(A,b,1e-6,20);
c_col=center(1);
c_row=center(2);

radiusX(1,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radiusX(2,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radiusX(3,1) = sqrt((c_col - x(1))^2 + (c_row-y(1))^2);
radius = mean(radiusX);

end

