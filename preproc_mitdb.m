function preproc_mitdb(signal_number)

% Sparse Coding with Anomaly Detection
% Amir Adler · Michael Elad · Yacov Hel-Or ·Ehud Rivlin

% clear all
% close all

load('data/mitdb/MITDB.mat');

% Select MitDB signal
%signal_number = 214;

signal_data     = data{signal_numbers == signal_number}(:,1);   % signal sample data
signal_Fs       = Fs{signal_numbers == signal_number};          % sampling frequency
signal_ann      = ann{signal_numbers == signal_number};         % vector of annotation positions
signal_type     = type{signal_numbers == signal_number};        % vector of annotation values
signal_type_str = [repmat('  ', size(signal_type,1),1) signal_type];  % prepend " " to annotation values for printing


% Signal preprocessing: BP filtering
[b,a] = butter(5,[1 100]/signal_Fs*2, 'bandpass');
x = filtfilt(b, a, signal_data);
plot([signal_data, x])

plot(signal_data);hold on;grid on
plot(signal_ann, signal_data(signal_ann), 'ro');
text(signal_ann, signal_data(signal_ann), signal_type_str);


% Split signal into segments of 5 minutes
seg_seconds = 300;       % 300 seconds = 5 minutes
seg_len = signal_Fs*300;  % number of samples in the segment

% Split sample indices into segments, then use them to segment the data and
% the annotations
signal_ind = 1:numel(signal_data);
signal_ind_matrix = buffer(signal_ind, seg_len);  % Segment with 'buffer()'; Segments = columns

% Take each segment
%seg_num = 3;              % which segment to take
%seg_idx = (seg_num-1)*seg_len+1:seg_num*seg_len;  % indices in x
for seg_num = 1:size(signal_ind_matrix,2)
    seg_idx = signal_ind_matrix(:, seg_num);  % each segment = each column 

    % Ignore padding with 0 at the end of last segment
    seg_idx = seg_idx(seg_idx~=0);
    
    % Segment data
    seg_x   = x(seg_idx);
    
    % Annotation indices in the segment, indexed from start of segment
    annidx = signal_ann>seg_idx(1)& signal_ann <= seg_idx(end);
    seg_ann = signal_ann(annidx);
    seg_ann = seg_ann - seg_idx(1);
    
    % Annotation values in the segment
    seg_type     = signal_type(annidx);
    seg_type_str = signal_type_str(annidx);  % prepend space ' '
    %
    plot(seg_x);hold on;grid on
    %plot(seg_ann, seg_x(seg_ann), 'ro');
    %text(seg_ann, seg_x(seg_ann),seg_type_str);
    %
    seg_ann_P      = seg_ann(seg_type~='L');
    seg_type_str_P = seg_type_str(seg_type~='L');
    plot(seg_ann_P, seg_x(seg_ann_P), 'ro');
    text(seg_ann_P, seg_x(seg_ann_P),seg_type_str_P);
    hold off

    % Get windows of segment
    winlen = 256;
    % wincount = length(seg_x) - winlen + 1;
    % win = zeros(winlen, wincount);
    % for i = 1:wincount
    %     win(:,i) = seg_x(i : i+winlen-1);
    % end
    win = buffer(seg_x, winlen, winlen-1, 'nodelay');  % buffer with overlap len-1
    P = pca(win');
    PCAsize = 32;
    Pr = P(:,1:PCAsize);

    %plot(Pr*Pr'*win(:,1:10))
    Cwin = Pr'*win;

    filename = sprintf('data/preproc/preproc_mitdb%d_seg%d.mat', signal_number, seg_num);
    fprintf('Saving %s\n', filename);
    save(filename, 'Cwin', 'Pr', 'win', 'seg_type', 'seg_ann');
end