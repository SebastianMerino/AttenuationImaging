%% Sample matrix X
M = 10;
N = 6;
%M = 200;
%N = 120;

X = ones(M,N);
X(1:floor(M/2),1:floor(N/2)) = 0;
% X = randn(M,N);
figure, imagesc(X)
colorbar
axis equal
axis tight

%% Difference operators
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dx = kron(D,speye(M));

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dy = kron(speye(N),D);

figure, tiledlayout(1,2),
nexttile,
imagesc(Dx)
colorbar
axis image

nexttile,
imagesc(Dy)
colorbar
axis image

%% Gradients
Dh = Dx*X(:);
Dv = Dy*X(:);

figure, tiledlayout(1,2),
nexttile,
imagesc(reshape(Dh,[M N]))
colorbar
axis image

nexttile,
imagesc(reshape(Dv,[M N]))
colorbar
axis image

%% Sample matrix X
M = 200;
N = 120;

X = randn(M,N);
figure, imagesc(X)
colorbar
axis equal
axis tight

%% Gradients
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,end) = [];
D(N,N) = 0;
Dx = kron(D,speye(M));

D = spdiags([-ones(M,1) ones(M,1)], [0 1], M,M+1);
D(:,end) = [];
D(M,M) = 0;
Dy = kron(speye(N),D);
    
Dh = Dx*X(:);
Dv = Dy*X(:);

figure, tiledlayout(1,2),
nexttile,
imagesc(reshape(Dh,[M N]))
colorbar
axis image

nexttile,
imagesc(reshape(Dv,[M N]))
colorbar
axis image

%% Weights
w = ones(M,N);
w(1:M/2,1:N/2) = 0.5;
W = spdiags(w(:),0,M*N,M*N);

Dh = W*Dx*X(:);
Dv = W*Dy*X(:);

figure, tiledlayout(1,3),
nexttile,
imagesc(w)
colorbar
axis image

nexttile,
imagesc(reshape(Dh,[M N]))
colorbar
axis image

nexttile,
imagesc(reshape(Dv,[M N]))
colorbar
axis image



