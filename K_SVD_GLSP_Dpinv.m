%function [D, residuals, D_hist, coef_hist, GLSP_coef_hist, GLSP_supp_hist] = K_SVD_GLSP(data, D0, niter, L)
function [D, residuals, D_hist, coef_hist] = K_SVD_GLSP(data, D0, niter, L)

D = D0;
OMPresiduals = [];
GLSPresiduals = [];
D_hist = {};
coef_hist = {};
GLSP_coef_hist = {};
GLSP_supp_hist = {};


for iiter=1:niter
    fprintf('%s --- K-SVD_GLSP iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    
    % OMP Sparse Coding preconditioned
    %[coef, R] = GLSP_via_OMP(data, D, L);
    
    Dpinv = pinv(D);
    newdata = Dpinv*data;   % preconditioned data
    [Ui,Si,Vi] = svd(D,'econ');
    newD = Vi*Vi'; %          % preconditioned dictionary, equivalent to pinv(D)*D

    % OMP without atom normalization (from Blumensath's 'sparsify' toolbox
    for i = 1:size(newdata,2)
        coef(:,i) = greed_omp(newdata(:,i), newD, size(newD,2), 'solver', 'qr', 'stopCrit', 'M','stopTol',L);
    end    
    
    % Use Python GLSP --- SLOW
%     np_result = cell(py.GLS.greedy_least_squares(mat2np(data), mat2np(D), L));
%     GLSP_coef = np2mat(np_result{1});
%     GLSP_supp = np2mat(np_result{2});
%     %GLSP_debias = np2mat(np_result{3});
%     %GLSP_recerr = np2mat(np_result{4});
%     %coef = GLSP_coef(logical(GLSP_supp));
%     coef = GLSP_coef.*GLSP_supp;
    %fprintf('%s --- --- Residual = %g\n', datestr(now, 'yy-mm-dd -HH:MM:SS'), norm(GLSP_coef.*(1-GLSP_supp),'fro'));

    % Compute residual matrix
    R = newdata - newD*coef;
    
    %residuals(iiter) = norm(GLSP_coef.*(1-GLSP_supp),'fro');
    residuals(iiter) = norm(R,'fro');      
    fprintf(', residual = %g\n', residuals(iiter));
    
    %OMPcorr = newD' * R;    
    %GLSPresiduals(iiter) = norm(OMPcorr,'fro');      
    %fprintf(', GLSP residual = %g\n', GLSPresiduals(iiter));
    
    % Coordinate descent - update each atom one by one
    for i = 1:size(newD,2)

        % Keep contribution of current atom
        Rcurr = R + newD(:,i)*coef(i,:); % 
        
        % Restrict R to the data having the current atom in the support
        Rcurr = Rcurr(:, coef(i,:) ~= 0);
        
        if (size(Rcurr,2) > 0)
            % SVD
            [U,~,~] = svd(Rcurr,'econ');
            newD(:,i) = U(:,1);
        end
    end
    
    % Undo preconditioning and find real % check in row space: D norm(newD - Dpinv*D*newD, 'fro')
    %D = V'*newD;
    D = D*newD;
    
    % Normalize D?
    for i = 1:size(D,2), D(:,i) = D(:,i)/norm(D(:,i)); end
    
    % Find nearest SPD matrix
%     A = nearestSPD(newD);
%     [UA, SA, VA] = svd(A, 'econ');
%     assert( norm(UA(:,1:32) - VA(:,1:32), 'fro') < 1e-12);
%     D = Ui * Si * (UA(:,1:32)');
    
    % Save history
    D_hist{iiter} = D;
    coef_hist{iiter} = coef;
    %GLSP_coef_hist{iiter} = GLSP_coef;
    %GLSP_supp_hist{iiter} = GLSP_supp;
end 