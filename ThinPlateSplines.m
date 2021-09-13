function [g, df] = ThinPlateSplines(data, lambda, meshStart, meshStep, meshEnd)
    %This function takes a dataset, value of the regularization parameter
    %lambda as well as parameters meshStart, meshStep and meshEnd to create
    %a square mesh plot and creates an optimizing function using Thin Plate 
    %Splines regression, and generates the function on the created square mesh
    
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


    %Concatenate matrices and obtain the parameter estimates including
    %the number of degrees of freedom
    KP = [K P'];
    K = K + lambda * eye(size(K, 1));
    M = [K, P'];
    N = [P, zeros(size(P, 1), size(X, 2) + 1)];
    M = [M;N];
    yExt = [Y;zeros(size(P, 1), 1)];
    parEst = M\yExt;
    alphaEst = parEst(1:size(X, 1));
    betaEst = parEst(size(X, 1) + 1:end);
    dfMat = KP / M;
    df = trace(dfMat(:, 1:size(X, 1)));


    %Define the x-y grid and create a function over the defined grid
    x = meshStart:meshStep:meshEnd;
    y = x;
    [xMesh, yMesh] = meshgrid(x);
    %Create minimizing function
    F = betaEst(1) + betaEst(2) * xMesh + betaEst(3) * yMesh;
    for i = 1:size(F, 1)
        for j = 1:size(F, 2)
            for k = 1:length(alphaEst)
                meshVec = [x(j), y(i)];
                F(i, j) = F(i, j) + alphaEst(k) * Eta(norm(meshVec - X(k, :)));
            end
        end
    end
    g = F;
end