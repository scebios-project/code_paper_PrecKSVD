function [D, residuals, D_hist, coef_hist, R_hist] = K_SVD_GLSP_frameDiag(data, D0, niter, L)
% Dictionary learning for GLSP
%
% 1. Run GLSP (apply preconditioner, run OMP without normalization)
% 2. Find Nearest Tight Frame and Diag, with Procrustes orthogonal but not orthonormal
% 3. Find nearest preconditioner with L2 minimization (general matrix)
% 4. Normalize
%



D_hist = {};
coef_hist = {};
R_hist = {};

% Normalize D0
D0 = D0 / norm(D0, 'fro');

D = D0;

% Run half of iteration (iiter+1) as well, to obtain the coefficients
for iiter=1:(niter+1)
    fprintf('%s --- K-SVD FrameDiag iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
    % Compute preconditioner = S^(-1)*U';
    [U,S,~] = svd(D,'econ');
    Prec = S^(-1)*U';
    
    % Apply preconditioner to data and dictionary
    precdata = Prec*data;  % preconditioned data
    precD = Prec*D;        % preconditioned dictionary

    % Run GLSP = OMP without atom normalization (from Blumensath's 'sparsify' toolbox)
    for i = 1:size(precdata,2)
        coef(:,i) = greed_omp(precdata(:,i), precD, size(precD,2), 'solver', 'qr', 'stopCrit', 'M','stopTol',L);
    end    

    % Compute residual matrix
    R = precdata - precD*coef;
    residuals(iiter) = norm(R,'fro');      
    fprintf(', residual = %g\n', residuals(iiter));    

    % Save dictionary and coefficients in history
    %   D{1} = initial dictionary D0
    %   coef_hist{1} = coefficients for initial dictionary D0
    D_hist{iiter} = D;
    coef_hist{iiter} = coef;
    R_hist{iiter} = sqrt(sum(R.^2));
    
    % Stop after computing the coefficients and residual for the final dictionary
    if iiter == (niter+1)
        break;
    end    
    
    unused_atoms = find(sum(abs(coef),2) == 0);
    fprintf('%s --- --- %d unused atoms\n',  datestr(now, 'yy-mm-dd HH:MM:SS'), numel(unused_atoms));
    
    %========================================
    % Update part 1: nearest tight frame * Diag
    %newTF = Procrustes_Orthonormal(precdata, coef);
    [~, newTF, D] = Procrustes_Orthogonal_pre(precdata, coef);
    coef = D*coef;
    %========================================
    
    %========================================
    % Update part 2: nearest preconditioner
    % minimize ||Prec*data - newTF*coef|| = 
    %   = minimize || newTF*coef - Prec*data||
    %   = minimize ||coef'newTF' - data'*Prec'||
    %   = minimize ||coef'newTF' - data'*V*D ||
    % 
    [VD, V, D] = Procrustes_Orthogonal_post(coef'*newTF', data');
    d = diag(D);
    newPrec = (V*D)';
    %========================================
    
    % New Dictionary = newPreconditioner^(-1)* newTightFrame, normalized
    %D = newPrec^(-1)*newTF;
    dinv = d.^-1 / norm(d.^-1, 2);
    D = V*diag(dinv)*newTF;  % newPrec = D'*V', newPrec^-1 = V*D^-1
    
    % Normalize dictionary, to prevent it getting bigger
    %D = D / norm(D, 'fro') * norm(D0, 'fro');
end 
