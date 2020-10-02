% Test equivalence of GLSP and preconditioned-OMP, using S^-1UT as
% preconditioner
clear all
close all

load('preproc_mitdb109.mat')
%load('Dicts_iter50_L4.mat')
load('Dicts_iter50_L8.mat')

% Pick a dictionary
DGLSP = DGLSP_hist{12};

% Sparse Coding Parameters
L = 5;

%% GLSP
np_result = cell(py.GLS.greedy_least_squares(mat2np(Cwin), mat2np(DGLSP), L));
GLSPcoef = np2mat(np_result{1});

%% OMP
param.L  = L;
[U,S,V] = svd(DGLSP,'econ');
newdata = S^(-1)*U'*Cwin;
newD = V'; % equivalent to S^(-1)*U'*D

% YES:
% Use OMP without atom normalization!
for i = 1:size(newdata,2)
	OMPcoef(:,i) = greed_omp(newdata(:,i), newD, size(newD,2), 'solver', 'qr', 'stopCrit', 'M','stopTol',L);
end

OMPres = newdata - newD*OMPcoef;
OMPcorr = newD' * OMPres;

% Check equality -> difference should be very small
norm(GLSPcoef - (OMPcoef + OMPcorr), 'fro')

% Check that OMPcorr is in the row space of D (orthogonal to null space) -> difference should be very small
norm(OMPcorr - pinv(newD)*newD*OMPcorr, 'fro')
