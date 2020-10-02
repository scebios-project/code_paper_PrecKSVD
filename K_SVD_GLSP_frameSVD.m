%function [D, residuals, D_hist, coef_hist, GLSP_coef_hist, GLSP_supp_hist] = K_SVD_GLSP(data, D0, niter, L)
function [D, residuals, D_hist, coef_hist, R_hist] = K_SVD_GLSP_frame(data, D0, niter, L)

D_hist = {};
coef_hist = {};
R_hist = {};

% Normalize D0
D0 = D0 / norm(D0, 'fro');

D = D0;

% Run half of iteration (iiter+1) as well, to obtain the coefficients
for iiter=1:(niter+1)
    fprintf('%s --- K-SVD_GLSP iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
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
    
    %========================================
    % Update part 1: nearest tight frame
    newTF = Procrustes_Orthonormal(precdata, coef);
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
    
    % Normalize dictionary atoms
    %D = D / norm(D, 'fro') * norm(D0, 'fro');
end 


function 
