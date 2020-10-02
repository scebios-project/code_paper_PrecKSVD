function [coef, resid] = GLSP_via_OMP(data, D, L)
Dpinv = pinv(D);
newdata = Dpinv*data;   % preconditioned data
[~,~,Vi] = svd(D,'econ');
newD = Vi*Vi'; %          % preconditioned dictionary, equivalent to pinv(D)*D

% OMP without atom normalization (from Blumensath's 'sparsify' toolbox
for i = 1:size(newdata,2)
    coef(:,i) = greed_omp(newdata(:,i), newD, size(newD,2), 'solver', 'qr', 'stopCrit', 'M','stopTol',L);
end    

resid = newdata - newD*coef;

