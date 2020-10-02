function [D, residuals, D_hist, coef_hist, R_hist] = K_SVD_GLSP_KsvdPrec(data, D0, niter, L)
% K-SVD in preconditioned space
%
% 1. A apply preconditioner to dictionary and data
% 2. Run K-SVD (OMP without normalization, followed by normal K-SVD update)
%     - note that apply preconditioner + OMP = GLSP
% 3. Undo preconditioner
% 4. Normalize the atoms


D_hist = {};
coef_hist = {};
R_hist = {};

ndata = size(data,2);

% Normalize D0
%D0 = D0 / norm(D0, 'fro');
for i = 1:size(D0,2), D0(:,i) = D0(:,i) / norm(D0(:,i)); end

D = D0;

% Run half of iteration (iiter+1) as well, to obtain the coefficients
for iiter=1:(niter+1)
    fprintf('%s --- K-SVD Prec iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
    % ============= Run sparse coding ==============
    
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
    
    % Coordinate descent - update each atom one by one
    for i = 1:size(precD,2)

        % Keep contribution of current atom
        Rcurr = R + precD(:,i)*coef(i,:); % 
        
        % Restrict R to the data having the current atom in the support
        Rcurr = Rcurr(:, coef(i,:) ~= 0);
        
        if (size(Rcurr,2) > 0)
            % SVD
            [U,~,~] = svd(Rcurr,'econ');
            precD(:,i) = U(:,1);
        end
    end    
    
    D = Prec^(-1) * precD;
    
  
    % Unit norm
    for i = 1:size(D,2), D(:,i) = D(:,i) / norm(D(:,i)); end
    
    
    % Normalize dictionary, to prevent it getting bigger
    %D = D / norm(D, 'fro') * norm(D0, 'fro');
end 
