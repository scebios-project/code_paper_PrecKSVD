function [D, residuals, D_hist, coef_hist, GLSP_coef_hist, GLSP_supp_hist] = K_SVD_2(useGLSP, data, D0, niter, L)

% Use least-SVD of residuals!!

D = D0;
residuals = [];
D_hist = {};
coef_hist = {};
GLSP_coef_hist = {};
GLSP_supp_hist = {};


for iiter=1:niter
    fprintf('%s --- K-SVD iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);
    if useGLSP
        np_result = cell(py.GLS.greedy_least_squares(mat2np(data), mat2np(D), L));
        GLSP_coef = np2mat(np_result{1});
        GLSP_supp = np2mat(np_result{2});
        %GLSP_debias = np2mat(np_result{3});
        %GLSP_recerr = np2mat(np_result{4});
        
        %coef = GLSP_coef(logical(GLSP_supp));
        coef = GLSP_coef.*GLSP_supp;
        %fprintf('%s --- --- Residual = %g\n', datestr(now, 'yy-mm-dd -HH:MM:SS'), norm(GLSP_coef.*(1-GLSP_supp),'fro'));
        residuals(iiter) = norm(GLSP_coef.*(1-GLSP_supp),'fro');
    else
        param.L = L;
        coef = mexOMP(data, D, param);
        residuals(iiter) = norm(data - D*coef,'fro');      
    end
    
    fprintf(', residual = %g\n', residuals(iiter));
    
    % Compute residual matrix
    R = data - D*coef;
    
    if useGLSP
        % Least-SVD iterations for K-SVD
        % Restrict R to the data having the current atom in the support
        Rcurr = Rcurr(:, coef(i,:) ~= 0);
        
    else
        % Normal K-SVD iteration
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
    
    % Save history
    D_hist{iiter} = D;
    coef_hist{iiter} = coef;
    if useGLSP
        GLSP_coef_hist{iiter} = GLSP_coef;
        GLSP_supp_hist{iiter} = GLSP_supp;
    end
end 