function classif(signal_number, segments, N, detclass)
% Run classification with all dictionaries and algorithms, and save outputs

%

% clear all
% close all

% foldername = 'D:\CS\ECG_Classif_MitBih\data\figs\';
% filebase = '2020-08-10_16-13-10';
% 
% filename = [foldername filebase '.mat'];
% load(filename)

% %% Select data files
% signal_number = 109;
% segments = 7;
% N = 64;

%% Classify
X = {}; Y = {}; T = {}; Targets_all = {}; Outputs_all = {};
% [X{1}, Y{1}, T{1}, Targets_all{1}, Outputs_all{1}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat', 'OMP', 'DOMP', true, 'OMP + KSVD');
% [X{2}, Y{2}, T{2}, Targets_all{2}, Outputs_all{2}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', true, 'LSP + KSVD Preconditioned');
% [X{3}, Y{3}, T{3}, Targets_all{3}, Outputs_all{3}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsKsvdSimple_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', true, 'LSP + KSVD Simple');
% [X{4}, Y{4}, T{4}, Targets_all{4}, Outputs_all{4}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsFrameDiag_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', true, 'LSP + Frame Diag');
% [X{5}, Y{5}, T{5}, Targets_all{5}, Outputs_all{5}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', true, 'LSP + Any Prec Frame Diag');
% [X{6}, Y{6}, T{6}, Targets_all{6}, Outputs_all{6}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat', 'GLSP', 'DOMP', true, 'LSP + KSVD');
% [X{7}, Y{7}, T{7}, Targets_all{7}, Outputs_all{7}] = ...
%     classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat', 'OMP', 'DGLSP', true, 'OMP + KSVD Prec');

% For sprintfc(), see http://undocumentedmatlab.com/articles/sprintfc-undocumented-helper-function
format_values = [signal_number*ones(segments, 1)  (1:segments)'];  % [106 1; 106 2; 106 3; ...]
data_files = sprintfc('data/preproc/preproc_mitdb%d_seg%d.mat', format_values);

format_values = [signal_number*ones(segments, 1)  (1:segments)'  N*ones(segments,1)];  % [106 1 64; 106 2 64; 106 3 64; ...]
param_strings = ...
{'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat',         'OMP',  'DOMP',  'OMP + KSVD';
'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat',         'GLSP', 'DGLSP', 'LSP + KSVD Preconditioned';...
'data/dicts/mitdb%d_seg%d_DictsKsvdSimple_N%d_iter40_L4.mat',       'GLSP', 'DGLSP', 'LSP + KSVD Simple';
'data/dicts/mitdb%d_seg%d_DictsFrameDiag_N%d_iter40_L4.mat',        'GLSP', 'DGLSP', 'LSP + Frame Diag';
'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_NoReplace_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', 'LSP + Any Prec Frame Diag';
'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat',          'GLSP', 'DOMP',  'LSP + KSVD';
'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat',         'OMP',  'DGLSP', 'OMP + KSVD Prec';
};
for i = 1:size(param_strings, 1)
    dict_files = sprintfc(param_strings{i}, format_values);
    [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = classifGenerateRoc(data_files, dict_files, detclass, param_strings{i,2}, param_strings{i,3}, true, param_strings{i,4});
end

% Save
datetime = datestr(now, 'YYYY-mm-dd_HH-MM-SS');
savename = sprintf('data/figs/mitdb%d_%dseg_N%d_%s.mat', signal_number, segments, N, datetime);
save(savename);

end