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
    %BAT = precdata * coef';
    %[U_BAT, S_BAT, V_BAT] = svd(BAT, 'econ');
    %newTF = U_BAT * V_BAT';
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
%     % Initial V1.  V1 = procrustes(coef'newTF' - data'*V) = procrustes(newTF*coef - V'*data)^T
%     d_hist{1} = ones(1, size(data,1));
%     [U,S,V] = svd(newTF*coef*data', 'econ');
%     V1 = (U*V')';
%     V1_hist{1} = V1;
%     Resi(1) = norm(coef'*newTF' - data'*V1, 'fro');
%     for i = 2:10
%         % Iterate
%         % 1. obtain D values
%         %for k = 1:size(V1,1)
%         %    num = (newTF(k,:)*coef)*data'*V1(:,k);
%         %    den = V1(:,k)'*data*data'*V1(:,k);
%         %    d(k) = num/den;
%         %end
%         numer = (newTF*coef)*data'*V1;
%         denom = V1'*(data*data')*V1;
%         for k = 1:size(V1,1), d(k) = numer(k,k)/denom(k,k); end
%         d_hist{i} = d;
%         % 2. Obtain new V = procrustes(coef'newTF'D^-1 - data'*V) = procrustes(D^(-1)*newTF*coef - V'*data)^T
%         %[U,S,V] = svd(data'*diag(1./d)*newTF*coef, 'econ');
%         [U, S, V] = svd(diag(1./d)*newTF*coef*data');
%         V1 = (U*V')';
%         V1_hist{i} = V1;
%         
%         % Compute residual
%         Resi(i) = norm(coef'*newTF' - data'*V1*diag(d), 'fro');
%         if Resi(i) > Resi(i-1)
%             %d = d_hist{i-1};
%             %V1 = V1_hist{i-1};
%             %newPrecT = V1*diag(d);
%             %newPrec = newPrecT';
%             i = i-1;  % make final i = i-1
%             break
%         end
%     end
%     finalV1 = V1_hist{i};
%     finald = d_hist{i};
%     newPrecT = finalV1*diag(finald);
%     newPrec = newPrecT';
    %========================================
    
    % New Dictionary = newPreconditioner^(-1)* newTightFrame, normalized
    %D = newPrec^(-1)*newTF;
    dinv = d.^-1 / norm(d.^-1, 2);
    D = V*diag(dinv)*newTF;  % newPrec = D'*V', newPrec^-1 = V*D^-1
    
    
    % Normalize dictionary, to prevent it getting bigger
    %D = D / norm(D, 'fro') * norm(D0, 'fro');
end 
