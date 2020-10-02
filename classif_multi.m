function classif_multi(signal_number, segments, N, detclass)
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

%% Classify with multi-channel strategies
X = {}; Y = {}; T = {}; Targets_all = {}; Outputs_all = {};

% For sprintfc(), see http://undocumentedmatlab.com/articles/sprintfc-undocumented-helper-function
format_values = [signal_number*ones(segments, 1)  (1:segments)'];  % [106 1; 106 2; 106 3; ...]
data_files = sprintfc('data/preproc/preproc_mitdb%d_multi_seg%d.mat', format_values);

format_values = [signal_number*ones(segments, 1)  (1:segments)'  N*ones(segments,1)];  % [106 1 64; 106 2 64; 106 3 64; ...]
param_strings = ...
{
'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_simult_N%d_iter40_L4.mat',        'sOMP',  'DOMP',  'sOMP + sKSVD';
'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_double_N%d_iter40_L4.mat',        'dOMP',  'DOMP',  'dOMP + KSVD';
'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_concat_N%d_iter40_L4.mat',        'OMPconcat',  'DOMP',  'OMP + KSVD concat';
'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_separ_N%d_iter40_L4.mat',         'OMPsep',  'DOMP',  'OMP + KSVD separate';
% 'data/dicts/mitdb%d_seg%d_DictsKsvdSimple_N%d_iter40_L4.mat',       'GLSP', 'DGLSP', 'LSP + KSVD Simple';
% 'data/dicts/mitdb%d_seg%d_DictsFrameDiag_N%d_iter40_L4.mat',        'GLSP', 'DGLSP', 'LSP + Frame Diag';
% 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_NoReplace_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', 'LSP + Any Prec Frame Diag';
% 'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat',          'GLSP', 'DOMP',  'LSP + KSVD';
% 'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat',         'OMP',  'DGLSP', 'OMP + KSVD Prec';
};
for i = 1:size(param_strings, 1)
    dict_files = sprintfc(param_strings{i,1}, format_values);
    [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = classifGenerateRoc(data_files, dict_files, detclass, param_strings{i,2}, param_strings{i,3}, true, param_strings{i,4});
end

%% Classify with single-channel strategies

% For sprintfc(), see http://undocumentedmatlab.com/articles/sprintfc-undocumented-helper-function
format_values = [signal_number*ones(segments, 1)  (1:segments)'];  % [106 1; 106 2; 106 3; ...]
data_files = sprintfc('data/preproc/preproc_mitdb%d_seg%d.mat', format_values);

format_values = [signal_number*ones(segments, 1)  (1:segments)'  N*ones(segments,1)];  % [106 1 64; 106 2 64; 106 3 64; ...]
param_strings = ...
{'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat',         'OMP',  'DOMP',  'OMP + KSVD';
'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat',         'GLSP', 'DGLSP', 'LSP + KSVD Preconditioned';...
%'data/dicts/mitdb%d_seg%d_DictsKsvdSimple_N%d_iter40_L4.mat',       'GLSP', 'DGLSP', 'LSP + KSVD Simple';
%'data/dicts/mitdb%d_seg%d_DictsFrameDiag_N%d_iter40_L4.mat',        'GLSP', 'DGLSP', 'LSP + Frame Diag';
%'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_NoReplace_N%d_iter40_L4.mat', 'GLSP', 'DGLSP', 'LSP + Any Prec Frame Diag';
%'data/dicts/mitdb%d_seg%d_DictsKsvdOMP_N%d_iter40_L4.mat',          'GLSP', 'DOMP',  'LSP + KSVD';
%'data/dicts/mitdb%d_seg%d_DictsKsvdPrec_N%d_iter40_L4.mat',         'OMP',  'DGLSP', 'OMP + KSVD Prec';
};
startval = numel(X)+1;
for i = startval:startval+size(param_strings,1)-1
    i1 = i-startval+1;
    dict_files = sprintfc(param_strings{i1,1}, format_values);
    [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = classifGenerateRoc(data_files, dict_files, detclass, param_strings{i1,2}, param_strings{i1,3}, true, param_strings{i1,4});
end


% Save
datetime = datestr(now, 'YYYY-mm-dd_HH-MM-SS');
savename = sprintf('data/figs/mitdb%d_%dseg_multi_N%d_%s.mat', signal_number, segments, N, datetime);
save(savename);

end