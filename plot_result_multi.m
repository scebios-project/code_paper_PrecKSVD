clear all
close all

foldername = 'D:\CS\ECG_Classif_MitBih\data\figs\';
filebases = {...
'mitdb106_7seg_multi_N64_2020-09-09_11-05-20'
'mitdb109_7seg_multi_N64_2020-09-09_11-15-37'
'mitdb116_7seg_multi_N64_2020-09-09_11-21-43'
'mitdb214_7seg_multi_N64_2020-09-09_11-27-41'
'mitdb228_7seg_multi_N64_2020-09-09_11-33-36'
};

for fi = 1:numel(filebases)
    filebase = filebases{fi};

    filename = [foldername filebase '.mat'];
    load(filename)

    figure
    h = plot(X{1}, Y{1}, X{2}, Y{2}, X{3}, Y{3}, X{4}, Y{4}, X{5}, Y{5}, X{6}, Y{6});
    set(h,{'LineWidth'},{2;2;2;2;2;2})
    set(h,{'LineStyle'},{'-';'--';':';'-.';'-';'-'})
    %set(h,{'Marker'},   {'+';'o';'*';'X'});
    set(h,{'MarkerIndices'},{1:100:length(Y{1})});

    titles = {  'Multi-channel: Simultaneous',...
                'Multi-channel: Concatenate',...
                'Multi-channel: Reunion',...
                'Multi-channel: Separate',...
                'Single-channel: KSVD+OMP',...
                'Single-channel: KSVD Prec + LSP'
                };
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
    paperWidth = 15;
    paperHeight = paperWidth*ratio;
    set(gcf, 'paperunits', 'centimeters');
    set(gcf, 'papersize', [paperWidth paperHeight]);
    set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
    print(gcf, '-dpdf', pdfname);
end