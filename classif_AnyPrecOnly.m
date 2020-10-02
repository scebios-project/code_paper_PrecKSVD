clear all
close all

% foldername = 'D:\CS\ECG_Classif_MitBih\data\figs\';
% filebase = '2020-08-10_16-13-10';
% 
% filename = [foldername filebase '.mat'];
% load(filename)


%% Select data files
%signal_number = 214;
signal_number = 106;
segments = 7;
N = 64;

%% Classify
X = {}; Y = {}; T = {}; Targets_all = {}; Outputs_all = {};

i=1; [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = ...[X{5}, Y{5}, T{5}, Targets_all{5}, Outputs_all{5}] = ...
    classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_N%d_iter40_L4.mat', 'GLSP', 'LSP + Any Prec Frame Diag', false);
i=2; [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = ...
    classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_N%d_iter40_L4.mat', 'GLSP', 'LSP + Any Prec Frame Diag', true);
i=3; [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = ...
    classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_NoReplace_N%d_iter40_L4.mat', 'GLSP', 'LSP + Any Prec Frame Diag', false);
i=4; [X{i}, Y{i}, T{i}, Targets_all{i}, Outputs_all{i}] = ...
    classifGenerateRoc(signal_number, segments, N, 'data/dicts/mitdb%d_seg%d_DictsAnyPrecFrameDiag_NoReplace_N%d_iter40_L4.mat', 'GLSP', 'LSP + Any Prec Frame Diag', true);


% datetime = datestr(now, 'YYYY-mm-dd_HH-MM-SS');
% savename = sprintf('data/figs/mitdb%d_%dseg_N%d_AnyFrameOnly_%s.mat', signal_number, segments, N, datetime);
% save(savename);