
% TV Andres Leonel Coila
function [TV] = TVcalc_isotropic(B,M,N,mask)

mask(isnan(mask)) = 0;
mask = mask(:);

X = reshape(B,M,N);
Dh = diff(X,[],1);
Dh = [Dh;zeros(1,N)];
Dv = diff(X,[],2);
Dv = [Dv zeros(M,1)];

P = Dh.^2 + Dv.^2;
%eps = 0.1;
%eps = 0.05;
%P = sqrt(P+eps^2);
P = sqrt(P);
TV = norm(P(:).*mask,1);

end


