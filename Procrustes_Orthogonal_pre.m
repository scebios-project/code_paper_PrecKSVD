function [U, V, D] = Procrustes_Orthogonal_pre(A, B)
% Solve orthogonal, but not orthonormal Procustes problem
% U = argmin_U ||A - U*B||_F
% where U = V*D, with V = orthonormal matrix and D = diagonal matrix
%
% Reference: "Orthogonal, but not Orthonormal, Procrustes Problems",
% Richard Everson, 1997
%
% Method used: Tandem algorithm (2.1)
% 

% Initialize d = 1 1 1 ....
d = ones(1, size(B,1));          % initial D = 1 1 1 ...
D = diag(d);

% Iterate: alternatively find best V, best D
for i = 1:10
    
    % 1. Obtain new V by minimizing || A - V*(D*B)||^2 with classical
    % orthonormal Procrustes
    % For Frobenius norm and trace of matrix, see See: https://math.stackexchange.com/a/2128704
    % For derivating trace of matrix, see: https://www.ics.uci.edu/~welling/teaching/KernelsICS273B/MatrixCookBook.pdf
    V = Procrustes_Orthonormal(A, D*B);
    V_hist{i} = V;
    
    % 2. Obtain D values by minimizing ||A - V*D*B||^2
    % Below not working because V'V is not I
    %numerator   = B*A'*V;      
    %denominator = sum(B.^2, 2);
    %denominator = sum(B.^2, 2) .* sum(V.^2, 1);  % divide by squared norm of B(k,:) * squared norm of V(:,k)
    %for k = 1:size(V,2), d(k) = numerator(k,k)/denominator(k); end
    %
    % See: https://mathoverflow.net/questions/75051/how-do-i-optimize-over-or-take-derivative-wrt-a-square-diagonal-matrix
    %d = optimizeDiag(A, V, B');
    %D = diag(d);
    % Takes very long
    %
    % Gradient descent
    %mu = 0.1;
    for iii = 1:20
        grad = diag(diag( (V'*V)*D*(B*B') - V'*A*B' ));
        mu = 0.05 / max(abs(diag(grad)));
        D = D - mu*grad;
        d = diag(D);
    end
    d_hist{i} = d;
    
    % Compute residual
    Residual(i) = norm(A - V*diag(d)*B, 'fro');
    %if (i > 1) && (Residual(i) > Residual(i-1))
    %    i = i-1;  % make final i = i-1
    %    break
    %end
end
% Final values
V = V_hist{i};
D = diag(d_hist{i});
U = V*D;
end

function w = optimizeDiag(P, X, Y)
% Solve w = argmin_w || P - X * diag(w) * Y^T ||_F
%
% See: https://mathoverflow.net/questions/75051/how-do-i-optimize-over-or-take-derivative-wrt-a-square-diagonal-matrix

% Build A
A = zeros(size(X,1)*size(Y,1), size(X,2));
for j = 1:size(X,2)
    Hj = X(:,j)*Y(:,j)';
    A(:,j) = Hj(:);
end

% Solve argmin_w ||vec(P) - A w ||_2
w = A \ P(:);
end