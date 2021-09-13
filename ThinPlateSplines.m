function g = ThinPlateSplines(data, lambda)
    %Extract X, Y matrices
X = data(:, 1:end - 1);
Y = data(:, end);

%Create matrix K
K = zeros(size(X, 1), size(X, 1));

for i = 1:size(X, 1)
    for j = 1:size(X, 1)
        K(i, j) = norm(X(i, :) - X(j, :));
        K(i, j) = Eta(K(i, j));
    end
end

%Create P matrix
P = ones(size(X, 2) + 1, size(X, 1));
P(2:end, :) = X';


%Concatenate matrices and obtain the parameter estimates
K = K + lambda * eye(size(K, 1));
M = [K, P'];
N = [P, zeros(size(P, 1), size(X, 2) + 1)];
M = [M;N];
yExt = [Y;zeros(size(P, 1), 1)];
parEst = M\yExt;
alphaEst = parEst(1:size(X, 1));
betaEst = parEst(size(X, 1) + 1:end);


%Plot the results
%Define the x-y grid
x = -1.5:0.1:1.5;
y = x;
[xMesh, yMesh] = meshgrid(x);
%Create minimizing function
F = betaEst(1) + betaEst(2) * xMesh + betaEst(3) * yMesh;
for i = 1:size(F, 1)
    for j = 1:size(F, 2)
        for k = 1:length(alphaEst)
            meshVec = [x(i), y(j)];
            F(i, j) = F(i, j) + alphaEst(k) * Eta(norm(meshVec - X(k, :)));
        end
    end
end
g = F;
end