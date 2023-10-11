
%% Differences for attenuations B and constants K
% m: number of rows
% n: number of columns (it must include the columns for B and C)
% Wx: weight-horizontal
% Wy: weight-vertical
function [B, Bt, BtB] = DiffOper(m,n,Wx,Wy)
D1 = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n+1); % col
D1(:,1) = [];
D1(1,1) = 0;
D2 = spdiags([-ones(m,1) ones(m,1)], [0 1], m,m+1); % row
D2(:,1) = [];
D2(1,1) = 0;
%B = [ kron(speye(m),D) ; kron(D,speye(m)) ];
B = [ Wx*kron(D1,speye(m)) ; Wy*kron(speye(n),D2)];

%B = [ kron(speye(m),D) ]; %Vertical
% B = [ kron(D,speye(m)) ]; % Horizontal
%B = [B B];
%B = [B sparse(size(B,1),size(B,2))];
Bt = B';
BtB = Bt*B;
end
