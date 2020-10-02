function [D, residuals, D_hist, coef_hist, R_hist] = K_SVD_GLSP_AnyPrecFrameDiag(data, D0, niter, L,...
                                                         replace_unused)
% Dictionary learning for GLSP
%
% 1. Run GLSP (apply preconditioner, run OMP without normalization)
% 2. Find Nearest Tight Frame and Diag, with Procrustes orthogonal but not orthonormal
% 3. Find nearest preconditioner with Nearest Orthogonal Matrixand Diag, with Procrustes orthogonal but not orthonormal
% 



D_hist = {};
coef_hist = {};
R_hist = {};

ndata = size(data,2);

% Normalize D0
D0 = D0 / norm(D0, 'fro');

D = D0;

% Run half of iteration (iiter+1) as well, to obtain the coefficients
for iiter=1:(niter+1)
    fprintf('%s --- K-SVD AnyPrecFrameDiag iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
    % ============= Run sparse coding ==============
    
    % Compute preconditioner = S^(-1)*U';
    [U,S,~] = svd(D,'econ');
    Prec = pinv(S)*U';    % Use pinv(S) instead of S^(-1) because may be singular
    
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
    % Update part 1: nearest tight frame * coef
    %newTF = Procrustes_Orthonormal(precdata, coef);
    [~, newTF, D] = Procrustes_Orthogonal_pre(precdata, coef);
    coef = D*coef;
    %========================================
    
    %========================================
%     % Update part 2: nearest preconditioner
%     % minimize ||Prec*data - newTF*coef|| = 
%     %   = minimize || newTF*coef - Prec*data||
%     %   = minimize ||coef'newTF' - data'*Prec'||
%     %   = minimize ||coef'newTF' - data'*V*D ||
%     % 
%     [VD, V, D] = Procrustes_Orthogonal_post(coef'*newTF', data');
%     d = diag(D);
%     newPrec = (V*D)';
%     
%     % ========================================
%     % New Dictionary = newPreconditioner^(-1)* newTightFrame, normalized
%     %
%     D = newPrec^(-1)*newTF;
%     dinv = d.^-1 / norm(d.^-1, 2);
%     D = V*diag(dinv)*newTF;  % newPrec = D'*V', newPrec^-1 = V*D^-1

    
    
    
    %========================================
    % Y = coef'*newTF'
    % A = data'
    % Y = A X => X = A \ Y
    newPrecGeneral = (data' \ (coef'*newTF'))';
    %========================================
    
    % New Dictionary = newPreconditioner^(-1)* newTightFrame, normalized
    D = pinv(newPrecGeneral)*newTF;   % use pinv() instead of inverse because may be singular
    [U,S,V] = svd(D, 'econ');    
    S = diag(diag(S) / norm(diag(S), 2));
    D = U*S*V';  % newPrec = D'*V', newPrec^-1 = V*D^-1

    %========================================
    % Replace unused atoms
    if replace_unused
        unused_atoms = find(sum(abs(coef),2) == 0);
        used_atoms = find(sum(abs(coef),2) ~= 0);
        median_norm = median(sqrt(sum(D(:,used_atoms).^2)));
        fprintf('%s --- --- Replacing %d unused atoms\n',  datestr(now, 'yy-mm-dd HH:MM:SS'), numel(unused_atoms));
        for i_atom = 1:numel(unused_atoms)
            i = unused_atoms(i_atom);

            % Pick a new atom
            newindex = randi(ndata);
            newatom = data(:, newindex);
            % Replace in dict
            D(:,i) = newatom / norm(newatom) * median_norm;
        end    
    end
    
    % Unit norm
    %for i = 1:size(D,2), D(:,i) = D(:,i) / norm(D(:,i)); end
    
    
    % Normalize dictionary, to prevent it getting bigger
    %D = D / norm(D, 'fro') * norm(D0, 'fro');
end 
