%% Utility functions
function x = optimTikhonovReg(A,B,Params)
% Choice of regularization function (Tikhonov matrix)
switch (Params.operator)
    case 'I'        % Penalizing greater SWS values
        L=speye(size(A,2));
    case 'G'        % Penalizing greater gradients
        L=speye(size(A,2))-circshift(speye(size(A,2)),[0 1]);
        L(end,:)=0;
    case 'L'
        L = Params.L;
end

x0 = zeros(size(A,2),1);
err = 1;
while err > Params.tolerance
    Lx = L*x0;
    W = spdiags( Params.k/2*( abs(Lx.^2+Params.beta).^(Params.k/2 - 1) ),...
        0, length(x0), length(x0));
    % x = ((A'*A+Params.alpha2 *L'*W*L)\A') *B;
    [x,~] = pcg( A'*A+Params.alpha2 *L'*W*L , A'*B, 1e-6 , 20,[],[],x0);
    err = norm(x-x0)^2/norm(x)^2; 
    x0 = x;
end

end