clear all
close all

foldername = 'D:\CS\ECG_Classif_MitBih\data\figs\';
filebases = {...
% 'mitdb106_7seg_N64_2020-09-06_15-16-03'
% 'mitdb109_7seg_N64_2020-09-06_15-30-14'
% 'mitdb116_7seg_N64_2020-09-06_15-44-45'
% 'mitdb214_7seg_N64_2020-09-06_16-14-43'
'mitdb228_7seg_N64_2020-09-06_21-40-25'
};

for fi = 1:numel(filebases)
    filebase = filebases{fi};

    filename = [foldername filebase '.mat'];
    load(filename)

    figure
    h = plot(X{1}, Y{1}, X{2}, Y{2}, X{3}, Y{3}, X{4}, Y{4}, X{5}, Y{5}, X{6}, Y{6}, X{7}, Y{7});
    set(h,{'LineWidth'},{1;1;1;1;1;1;1})
    set(h,{'LineStyle'},{'-';'-';'-';'--';'-';'-'     ;'-'})
    set(h,{'Marker'},   {'+';'o';'*';'.' ;'x';'square';'diamond'});
    set(h,{'MarkerIndices'},{1:100:length(Y{1})});

    titles = {  'KSVD + OMP',...
                'KSVD Prec + LSP',...
                '(LSP+SVD)+LSP',...
                'Frame Diag + LSP',...
                'Procrustes + LSP',...
                'KSVD + LSP',...
                'KSVD Prec + OMP'...
                };
    % h = plot(X{1}, Y{1}, X{6}, Y{6}, X{3}, Y{3}, X{2}, Y{2});
    % set(h,{'LineWidth'},{1;1;1;2})
    % titles = {  'KSVD + OMP',...
    %             'KSVD + LSP',...
    %             '(LSP+SVD) + LSP',...
    %             'Prec KSVD + LSP',...
    %             };
    axis([0,1,0.5,1]);
    legend(titles, 'Location', 'southeast')
    title('Receiver Operating Characteristics');


    % Save .fig
    figname = [foldername filebase '.fig'];
    savefig(figname)

    % Save PDF
    pdfname = [foldername 'fig_ROC_' filebase '.pdf'];
    ps = get(gcf, 'Position');
    ratio = (ps(3)-ps(1)) / (ps(4)-ps(2));
    paperWidth = 10;
    paperHeight = paperWidth*ratio;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
    print(gcf, '-dpdf', pdfname);
end