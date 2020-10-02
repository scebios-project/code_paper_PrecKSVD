function [X, Y, T, Targets_all, outputs_all] = classifGenerateRoc_simult(data_files, dict_files, detclass, SPalg, Dname, pick_best, titlestr)

% function [X, Y, T, Targets_all, outputs_all] = classifGenerateRoc(signal_number, segments, N, DictFilePattern, ...
%                                                                     SPalg, Dname, pick_best, titlestr)
% 
% for seg_num=1:segments
%     data_files{seg_num} = sprintf('data/preproc/preproc_mitdb%d_seg%d.mat', signal_number, seg_num);
%     %dict_files{seg_num} = sprintf('mitdb%d_seg%d_DictsAnyPrecFrameDiag_iter40_L4.mat', signal_number, seg_num);
%     %dict_files{seg_num} = sprintf('mitdb%d_seg%d_DictsAnyPrecFrameDiag_N%d_iter40_L4.mat', signal_number, seg_num, N);
%     dict_files{seg_num} = sprintf(DictFilePattern, signal_number, seg_num, N);
% end

Targets_all  = [];
outputs_all  = [];

fprintf('%s --- Running classification algorithm %s (simult)\n',  datestr(now, 'yy-mm-dd HH:MM:SS'), SPalg);

for ifile = 1:numel(data_files)

    fprintf('%s --- Loading file %s\n',  datestr(now, 'yy-mm-dd HH:MM:SS'), dict_files{ifile});
    
    load(data_files{ifile})
    load(dict_files{ifile})
    
    % Make data into a multi-dimensional tensot (3rd dimension = channels)
    Cwin_tensor = reshape(cell2mat(Cwin), size(Cwin{1},1), size(Cwin{1},2), numel(Cwin));

    % Sparse coding parameters
    param.L  = 4;

    % Select dictionary specified
    switch Dname
        case 'DOMP'
            D = DOMP;
        case 'DGLSP'
            if pick_best
                % Pick best dictionary
                [minval, minpos] = min(residGLSP);
                D = DGLSP_hist{minpos};
            else
                D = DGLSP;
            end
    end
    
    % Sparse coding wuth specified algorithm
    if strcmp(SPalg, 'OMP')
        %% Classification with OMP
        %OMPcoef = mexOMP(Cwin, D, param);
        for i = 1:size(Cwin_tensor, 2)  % 3rd dimension  = channels
           coef1(:,:,i) = greed_somp_qr(squeeze(Cwin_tensor(:,i,:)), D, size(D,2), 'stopCrit','M','stopTol',param.L);
        end
        OMPcoef = permute(coef1, [1 3 2]);

        % Find reconstruction error
        %winRec = Pr*D*OMPcoef;
        %recerror = (sum((win - winRec).^2));
        recerror = sum(sum((Cwin - mult3(D,OMPcoef)).^2, 3), 1);  % sum along channels and columns
    
    elseif strcmp(SPalg, 'GLSP')
        %% Classification with GLSP
        [GLSPcoef, GLSPresid] = GLSP_via_OMP(Cwin, D, param.L);
        recerror = sum(GLSPresid.^2);

    end

    %% Classify and plot ROC curve
    err = zeros(size(seg_type));
    Targets = (seg_type == detclass)';  % We want 1 when V is detected
    recerror_len = numel(recerror);
    for i = 1:numel(seg_ann)
        segment_i = seg_ann(i) - 127;
        left  = max(segment_i - 100, 1);
        right = min(segment_i + 100, recerror_len);
        if left < right  % it may happen than right < 0 so, skip this 
            err(i)  = max(recerror(left:right));
        else
            err(i) = 0;  % ignore
        end
    end
    outputs = err' / max(err);

    Targets_all  = [Targets_all Targets];
    outputs_all  = [outputs_all outputs];
   
end

figure
plotroc(Targets_all, outputs_all);
title(titlestr);

[X, Y, T] = perfcurve(Targets_all, outputs_all, 1);


function C = mult3(A,B)
% 3D multiplication A*B
for i=1:size(B,3)
    C(:,:,i) = A * B(:,:,i);
end