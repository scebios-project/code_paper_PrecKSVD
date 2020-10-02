function Omega = Procrustes_Orthonormal(B, A)
% Solve the orthonormal Procrustes problem
%  Omega = argmin ||B - Omega A||
 
[U,~,V] = svd(B*A', 'econ');
Omega = U*V';

end
