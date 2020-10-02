function dictlearn_multi(signal_number, segments, K, iter, L)

% clear all
% close all

%load('preproc_mitdb109.mat');

%% Select data files
%signal_number = 214;
%segments = 7;
for seg_num=1:segments
    data_files{seg_num} = sprintf('data/preproc/preproc_mitdb%d_multi_seg%d.mat', signal_number, seg_num);
end

parfor ifile = 1:numel(data_files)
    %% Load segment data
    ld = load(data_files{ifile});  % load 'Cwin' matrix
    Cwin = ld.Cwin;
    
    fprintf('Loaded %s\n', data_files{ifile});
    
    %% Dictionary learning
    %K = 64;
    %iter = 40;
    %L = 4;
    
    % Arrange data in different formats
    set_cell = Cwin;                                                                    % cell array of matrices
    set_concat = cell2mat(set_cell);                                                    % horizontal concatenated (more data)
    set_tensor = reshape(set_concat, size(set_cell{1},1), size(set_cell{1},2), []);     % multichannel tensor
    set_double = cell2mat(set_cell');                                                   % vertical concatenated (double-size data)
    
    % Initialize dictionary with random linear combinations, normalized
    %   Normal-size, from all data    
    D0_concat = set_concat*randn(size(set_concat,2), K);        % consider all data
    for i = 1:size(D0_concat,2), D0_concat(:,i) = D0_concat(:,i)/norm(D0_concat(:,i)); end
    %   Double-size
    D0_double = set_double*randn(size(set_double,2), K);        % consider all data
    for i = 1:size(D0_double,2), D0_double(:,i) = D0_double(:,i)/norm(D0_double(:,i)); end
    %   Separate
    D0_separ = {};
    for ich = 1:numel(set_cell)
        D0_separ{ich} = set_cell{ich}*randn(size(set_cell{ich},2), K);        % consider all data
        for i = 1:size(D0_separ{ich},2), D0_separ{ich}(:,i) = D0_separ{ich}(:,i)/norm(D0_separ{ich}(:,i)); end
    end
    
    
    % K-SVD variants
    %
%     % 1. K-SVD Prec: K-SVD in preconditioned space (GLSP + K-SVD update in preconditioned space)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_KsvdPrec(dltrainset, D0, iter, L);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdPrec_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');
% 
%     % 2. K-SVD Simple: GLSP + nomal K-SVD in original space (GLSP + K-SVD update in original space)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_KsvdSimple(dltrainset, D0, iter, L);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdSimple_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');
%     
%     % 3A. Dictionary learning with Procrustes (TF + Diag) + L2 minimization (preconditioner)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_AnyPrecFrameDiag(dltrainset, D0, iter, L, true);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsAnyPrecFrameDiag_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');
% 
%     % 3B. Dictionary learning with Procrustes (TF + Diag) + L2 minimization (preconditioner)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_AnyPrecFrameDiag(dltrainset, D0, iter, L, false);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsAnyPrecFrameDiag_NoReplace_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');    
%     
%     % 4. Dictionary learning with Procrustes (TF + Diag) + Procrustes (Orthon + Diag) (preconditioner)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_frameDiag(dltrainset, D0, iter, L);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsFrameDiag_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');

%     % 1. K-SVD Prec: K-SVD in preconditioned space (GLSP + K-SVD update in preconditioned space)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_KsvdPrec(dltrainset, D0, iter, L);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdPrec_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');


    % Simultaneous K-SVD with OMP
    [DOMP, residOMP, DOMP_hist, coefOMP_hist] = K_SVD_simult(set_tensor, D0_concat, iter, L);
    savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdOMP_simult_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
    %save(savename, 'DOMP','DOMP_hist','residOMP', 'iter','L');
    parsave(savename, DOMP, DOMP_hist, residOMP, iter, L);
    
    % Double-sized normal K-SVD with OMP
    [DOMP, residOMP, DOMP_hist, coefOMP_hist] = K_SVD(set_double, D0_double, iter, L);
    savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdOMP_double_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
    %save(savename, 'DOMP','DOMP_hist','residOMP', 'iter','L');
    parsave(savename, DOMP, DOMP_hist, residOMP, iter, L);
    
    % Concatenated normal K-SVD with OMP
    [DOMP, residOMP, DOMP_hist, coefOMP_hist] = K_SVD(set_concat, D0_concat, iter, L);
    savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdOMP_concat_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
    %save(savename, 'DOMP','DOMP_hist','residOMP', 'iter','L');
    parsave(savename, DOMP, DOMP_hist, residOMP, iter, L);    
    
    % Separate normal K-SVD with OMP for each channel
    DOMP={}; residOMP={};DOMP_hist={};coefOMP_hist={};
    for ich=1:numel(set_cell)
        [DOMP{ich}, residOMP{ich}, DOMP_hist{ich}, coefOMP_hist{ich}] = K_SVD(set_cell{ich}, D0_separ{ich}, iter, L);
    end
    savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdOMP_separ_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
    parsave(savename, DOMP, DOMP_hist, residOMP, iter, L);
    

    
%     % 1. K-SVD Prec: K-SVD in preconditioned space (GLSP + K-SVD update in preconditioned space)
%     [DGLSP, residGLSP, DGLSP_hist, coefGLSP_hist, R_hist] = K_SVD_GLSP_KsvdPrec(dltrainset, D0, iter, L);
%     savename = ['data/dicts/mitdb' num2str(signal_number) '_seg' num2str(ifile) '_DictsKsvdPrec_N' num2str(K) '_iter' num2str(iter) '_L' num2str(L) '.mat'];
%     save(savename, 'DGLSP', 'DGLSP_hist', 'residGLSP', 'iter','L');
    
    
end