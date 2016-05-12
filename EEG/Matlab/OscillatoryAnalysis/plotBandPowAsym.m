function plotBandPowAsym( bpasym )
%% FUNCTION plotBandPowAsym
% PURPOSE:  Plot the eeg2bandpwrs from an EEG.
% PARAMS:   'bpasym' - structure holding matrices of asymmetry results
% USAGE:    
% Note:     Requires the creation of the data by the function analysebandpwrs()
%           WIP: handle whatever protocols have been handed in
% CALLS:    handle graphics

    % Get num fields in bpasym
%     sz = sum(~structfun(@isempty, bpasym));
    % Get max y values so all figures are to the same scale
    maxy = max(max(cell2mat(struct2cell(bpasym))));
    miny = min(min(cell2mat(struct2cell(bpasym))));

    lbwh = (get(0,'ScreenSize')./1.2)+50;
    f1 = figure('Units','pixels','Position',lbwh);
    
    if isfield( bpasym, 'vigec' )
        subplot(2,2,1);    
        bar(bpasym.vigec');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EC Vigilance protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    if isfield( bpasym, 'vigeo' )
        subplot(2,2,2);    
        bar(bpasym.vigeo');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EO Vigilance protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    
    if isfield( bpasym, 'allec' )
        subplot(2,2,3);
        bar(bpasym.allec');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EC all protocols');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    if isfield( bpasym, 'alleo' )
        subplot(2,2,4);    
        bar(bpasym.alleo');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EO all protocols');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    
    hold all;
    f2 = figure('Units','pixels','Position',lbwh);
    
    if isfield( bpasym, 'helec' )
        subplot(2,2,1);
        bar(bpasym.helec');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EC Hello protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    if isfield( bpasym, 'heleo' )
        subplot(2,2,2);    
        bar(bpasym.heleo');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EO Hello protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    
    if isfield( bpasym, 'gbyec' )
        subplot(2,2,3);
        bar(bpasym.gbyec');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EC Goodbye protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end

    if isfield( bpasym, 'gbyec' )
        subplot(2,2,4);
        bar(bpasym.gbyeo');
        set( gca, 'YLim', [miny maxy] );
        title('Aggregate EO Goodbye protocol');
        legend({'ln(P3/P4)' 'ln(C3/C4)' 'ln(F3/F4)'});
        newXTickLabels( gca, {'\delta' '\theta' '\alpha' '\beta' '\gamma'}, true, 1, 'Bands' );
        ylabel('Positive=greater Left');
    end
    
    % Save out figures.
    print( f1, '-dpng', fullfile(pwd, 'ROIasym_AllProtocols') );
    print( f2, '-dpng', fullfile(pwd, 'ROIasym_Hello_Gdbye') );
    
    close all;
end