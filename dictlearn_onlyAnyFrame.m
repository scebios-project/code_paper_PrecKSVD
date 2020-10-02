clear all
close all

%load('preproc_mitdb109.mat');

%% Select data files
%signal_number = 214;
%signal_number = 106;
%signal_number = 109;
signal_number = 116;
segments = 7;
for seg_num=1:segments
    data_files{seg_num} = sprintf('data/preproc/preproc_mitdb%d_seg%d.mat', signal_number, seg_num);
end

for ifile = 1:numel(data_files)
    %% Load segment data
    load(data_files{ifile});  % load 'Cwin' matrix
    
    fprintf('Loaded %s\n', data_files{ifile});
    
    %% Dictionary learning
    dltrainset = Cwin;
    K = 64;
    iter = 40;
    L = 4;

    % Initialize dictionary with random linear combinations, normalized
    D0 = dltrainset*randn(size(dltrainset,2), K);
    for i = 1:size(D0,2), D0(:,i) = D0(:,i)/norm(D0(:,i)); end

    % K-SVD variants
   
    % 3A. Dictionary learning with Procrustes (TF + Diag) + L2 minimization (preconditioner)
    [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_AnyPrecFrameDiag(dltrainset, D0, iter, L, false);
    savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsAnyPrecFrameDiag_NoReplace_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
    save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');

end