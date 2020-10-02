function [D, residuals, D_hist, coef_hist] = K_SVD(data, D0, niter, L)

D = D0;
residuals = [];
D_hist = {};
coef_hist = {};

for iiter=1:niter
    fprintf('%s --- K-SVD iter %2d',  datestr(now, 'yy-mm-dd HH:MM:SS'), iiter);

    param.L = L;
    coef = mexOMP(data, D, param);
    residuals(iiter) = norm(data - D*coef,'fro');      
    
    fprintf(', residual = %g\n', residuals(iiter));
    
    % Compute residual matrix
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
    
    % Save history
    D_hist{iiter} = D;
    coef_hist{iiter} = coef;

    %========================================
    % Replace unused atoms
    unused_atoms = find(sum(abs(coef),2) == 0);
    used_atoms = find(sum(abs(coef),2) ~= 0);
    %median_norm = median(sqrt(sum(D(:,used_atoms).^2)));
    fprintf('%s --- --- Replacing %d unused atoms\n',  datestr(now, 'yy-mm-dd HH:MM:SS'), numel(unused_atoms));
    for i_atom = 1:numel(unused_atoms)
        i = unused_atoms(i_atom);
        
        % Pick a new atom
        newindex = randi(ndata);
        newatom = data(:, newindex);
        % Replace in dict
        %D(:,i) = newatom / norm(newatom) * median_norm;
        D(:,i) = newatom / norm(newatom);
    end    
    
end 