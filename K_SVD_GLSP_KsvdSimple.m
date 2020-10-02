function [D, residuals, D_hist, coef_hist, R_hist] = K_SVD_GLSP_KsvdSimple(data, D0, niter, L)
% K-SVD in preconditioned space
%
% 1. Run GLSP (apply preconditioner, run OMP without normalization)
% 2. K-SVD update for atoms, based on GLSP coefficients
%     - note that apply preconditioner + OMP = GLSP
% 
% Difference from KsvdPrec: dictionary update in original space, not in
% preconditioned space

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
    fprintf('%s --- K-SVD Simple iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
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
    
    % Simple KSVD: no preconditioned data for update
    R = data - D*coef;
    
    % Coordinate descent - update each atom one by one
    for i = 1:size(D,2)

        % Keep contribution of current atom
        Rcurr = R + D(:,i)*coef(i,:); % 
        
        % Restrict R to the data having the current atom in the support
        Rcurr = Rcurr(:, coef(i,:) ~= 0);
        
        if (size(Rcurr,2) > 0)
            % SVD
            [U,~,~] = svd(Rcurr,'econ');
            D(:,i) = U(:,1);
        end
    end    
end 
