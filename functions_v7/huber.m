function [tau] = huber(v,eps)

tau = 2*sqrt(v.^2 + eps^2);
%c = (abs(v) - eps^2) > 0;
%tau = c.*(v) + (1-c).*(eps^(-1)*v.^2);

end