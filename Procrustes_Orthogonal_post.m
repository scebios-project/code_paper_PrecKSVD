function [U, V, D] = Procrustes_Orthogonal_post(A, B)
% Solve orthogonal, but not orthonormal Procustes problem
% U = argmin_U ||A - B*U||_F
% where U = V*D, with V = orthonormal matrix and D = diagonal matrix
%
% Reference: "Orthogonal, but not Orthonormal, Procrustes Problems",
% Richard Everson, 1997
%
% Method used: Schwarz tandem
% 


d_hist{1} = ones(1, size(B,2));          % initial D = 1 1 1 ...

% Obtain (classical) V by minimizing ||A - B*V||^2
% i.e. minimizing ||A' - V'*B'||^2  = orthonormal Procrustes
V = Procrustes_Orthonormal(A', B');
V = V';
V_hist{1} = V;                           % save initial V
Residual(1) = norm(A - B*V, 'fro');      % initial residual  

% Iterate: alternatively find best D, best V
for i = 2:10
    
    % 1. Obtain Di+1 values by minimizing ||A - B*Vi*Di+1||^2
    numerator   = A'*B*V;      
    denominator = V'*(B'*B)*V;
    for k = 1:size(V,1), d(k) = numerator(k,k)/denominator(k,k); end
    d_hist{i} = d;
    
    % 2. Obtain new V by minimizing || A - B*V*Di+1 ||^2
    %    using Schwarz tandem
    [U2, ~, V2] = svd(diag(1./d)*A'*B);  % here I think it is a mistake in original paper, 
                                         % it's Di+1^-1*A'*B instead of B*Di+1^-1*A', 
                                         % since this is a transposed orthonormal Procrustes problem
    V = U2*V2';
    V = V';
    V_hist{i} = V;

    % Compute residual
    Residual(i) = norm(A - B*V*diag(d), 'fro');
    if Residual(i) > Residual(i-1)
        i = i-1;  % make final i = i-1
        break
    end
end
finalV = V_hist{i};
finald = d_hist{i};
U = finalV*diag(finald);
V = finalV;
D = diag(finald);