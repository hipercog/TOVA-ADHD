% INIT PARAMS AND DATA (if doing fresh timefreq) USING TOVA_paper1_loadnclean.m

%% IMPORT DATA AND COMBINE ALL CONDITIONS
% Skip this unless doing fresh time-freq calculations
cEEG = pop_mergeset(pop_mergeset(cEEGrsp, 1:numel(cEEGrsp))...
                  , pop_mergeset(cEEGinb, 1:numel(cEEGinb)));
cEEG = rmgrandbl(cEEG, NaN, 0);
aEEG = pop_mergeset(pop_mergeset(aEEGrsp, 1:numel(aEEGrsp))...
                  , pop_mergeset(aEEGinb, 1:numel(aEEGinb)));
aEEG = rmgrandbl(aEEG, NaN, 0);
clear cEEGrsp cEEGinb aEEGrsp aEEGinb


%% ERSPs per ROI
for r = 2%which_r
    roi = ROI{1, r};
    myReport(['Processing ' ROI{2, r} ' trodes: ' biosemi1020(roi)])
    
    svnm = ['TOVA_CvA_ersp-evok-indu-erps_makeig-BL_' ROI{2, r} '_cyc3-11'];
    if isfile(fullfile(oud, 'ersp', [svnm '.mat']))
        load(fullfile(oud, 'ersp', [svnm '.mat']))
    else
        %% get ERSP, Evoked, and Induced power from each ROI electrode
        bl = NaN;
        
        ersp = cell(numel(roi), 1);
        itc = cell(numel(roi), 1);
        blpow = cell(numel(roi), 1);
        eboot = cell(numel(roi), 1);
        iboot = cell(numel(roi), 1);
        alltf = cell(numel(roi), 1);

        evked = cell(numel(roi), 1);
        evktf = cell(numel(roi), 1);

        iducd = cell(numel(roi), 1);
        idutf = cell(numel(roi), 1);

        for t = 1:numel(roi)
            [ersp{t}, itc{t}, blpow{t}, times, freqs, eboot{t}, iboot{t}, alltf{t}]...
             = newtimef(...
                {cEEG.data(roi(t),:,:) aEEG.data(roi(t),:,:)}...
                , cEEG.pnts...
                , [cEEG.times(1) cEEG.times(end)]...
                , cEEG.srate...
                , [3 0.5] ...
                , 'plotersp', 'off'...
                , 'plotitc', 'off'...
                , 'freqs', [4 30] ...
                , 'nfreqs', 54 ...
                , 'freqscale', 'log'...
                , 'baseline', bl ... % 0 removes all negative values
                , 'commonbase', 'off'...
                , 'trialbase', 'full');
            [evked{t}, ~, ~, ~, ~, ~, ~, evktf{t}]...
             = newtimef(...
                {mean(cEEG.data(roi(t),:,:), 3) mean(aEEG.data(roi(t),:,:), 3)}...
                , cEEG.pnts...
                , [cEEG.times(1) cEEG.times(end)]...
                , cEEG.srate...
                , [3 0.5] ...
                , 'plotersp', 'off'...
                , 'plotitc', 'off'...
                , 'freqs', [4 30] ...
                , 'nfreqs', 54 ...
                , 'freqscale', 'log'...
                , 'baseline', bl ... % 0 removes all negative values
                , 'commonbase', 'off'...
                , 'trialbase', 'full');
            [iducd{t}, ~, ~, ~, ~, ~, ~, idutf{t}]...
             = newtimef({cEEG.data(roi(t),:,:) - mean(cEEG.data(roi(t),:,:), 3)...
                aEEG.data(roi(t),:,:) - mean(aEEG.data(roi(t),:,:), 3)}...
                , cEEG.pnts...
                , [cEEG.times(1) cEEG.times(end)]...
                , cEEG.srate...
                , [3 0.5] ...
                , 'plotersp', 'off'...
                , 'plotitc', 'off'...
                , 'freqs', [4 30] ...
                , 'nfreqs', 54 ...
                , 'freqscale', 'log'...
                , 'baseline', bl ... % 0 removes all negative values
                , 'commonbase', 'off'...
                , 'trialbase', 'full');

        end
        % SAVE CALCULATIONS
        save(fullfile(oud, 'ersp', [svnm '.mat']), 'ersp', 'itc' ,'blpow', 'times'...
                , 'freqs', 'eboot', 'iboot', 'evked', 'evktf', 'iducd', 'idutf'...
                , '-v7.3')
    end

    %% Average ERSP & ITC matrices from all ROI 'trodes
    cersp = mean(structarr2mat(ersp, 1), 3);
    aersp = mean(structarr2mat(ersp, 2), 3);
    cevkd = mean(structarr2mat(evked, 1), 3);
    aevkd = mean(structarr2mat(evked, 2), 3);
    ciduc = mean(structarr2mat(iducd, 1), 3);
    aiduc = mean(structarr2mat(iducd, 2), 3);

    Cerp = [mean(cersp); mean(ciduc); mean(cevkd) / 10];
    Aerp = [mean(aersp); mean(aiduc); mean(aevkd) / 10];
    CvAerp = [mean(cersp - aersp); mean(ciduc - aiduc); mean(cevkd - aevkd) / 10];


    %% Plot the ERSP-EVOKED-INDUCED figure
    figh = figure('Position', [1 1 lbwh(4) lbwh(4)], 'Color', 'w');%, 'Visible', false);

    % init params
    rws = 4;
    cls = 3;
    ts = times / 1000;
    xtls = ''; %(-500:100:500) / 1000;
    ytls = '';
    tpts = [15 101 186];
    ttl = Grup;

    % TOP ROW : Control - ADHD
    % get mirrored color map limit values as +-MAX VALUE
    clim = [-1 1] * max(abs([min(cat(3, cersp, aersp), [], 'all') ...
                             max(cat(3, cersp, aersp), [], 'all')]));
    subplot(rws, cls, 1)
    imagesclogy(ts, freqs, imgaussfilt(cersp, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'ydir', 'normal')
    title(ttl{1})

    ylabel(sprintf('ERSPs\n Frequency (Hz)'), 'FontWeight', 'bold')

    subplot(rws, cls, 2)
    imagesclogy(ts, freqs, imgaussfilt(aersp, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')
    title(ttl{2})

    subplot(rws, cls, 3)
    imagesclogy(ts, freqs, imgaussfilt(cersp - aersp, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')
    title([Grup{1} ' - ' Grup{2}])

    % Color bar and adjustment of it's subplot parent
    pltsz1 = get(gca, 'Position');
    cb = colorbar('EastOutside'); title(cb, 'ERSP (dB)', 'FontSize', 8)
    colormap(eval(clrmp{2}))
    pltsz2 = get(gca, 'Position');
    set(gca, 'Position', [pltsz2(1:2) pltsz1(3) pltsz2(4)])


    % MIDDLE ROW: INDUCED
    % get mirrored color map limit values
    subplot(rws, cls, 4)
    imagesclogy(ts, freqs, imgaussfilt(ciduc, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'ydir', 'normal')

    ylabel(sprintf('Induced power\n '), 'FontWeight', 'bold')
    
    subplot(rws, cls, 5)
    imagesclogy(ts, freqs, imgaussfilt(aiduc, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')

    subplot(rws, cls, 6)
    imagesclogy(ts, freqs, imgaussfilt(ciduc - aiduc, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')

    % COLOR BAR: INDUCED
    pltsz1 = get(gca, 'Position');
    cb = colorbar('EastOutside'); title(cb, 'Induced (dB)', 'FontSize', 8)
    colormap(eval(clrmp{2}))
    pltsz2 = get(gca, 'Position'); % adjust subplot parent
    set(gca, 'Position', [pltsz2(1:2) pltsz1(3) pltsz2(4)])


    % MIDDLE ROW 2: EVOKED
    % get mirrored color map limit values as +-3 STDEVs
    clim = [-1 1] * std(cat(3, cevkd, aevkd), 1, 'all') * 3;
    subplot(rws, cls, 7)
    imagesclogy(ts, freqs, imgaussfilt(cevkd, 2), clim); hold on
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'ydir', 'normal')
    
    ylabel(sprintf('Evoked power\n '), 'FontWeight', 'bold')

    subplot(rws, cls, 8)
    imagesclogy(ts, freqs, imgaussfilt(aevkd, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')

    subplot(rws, cls, 9)
    imagesclogy(ts, freqs, imgaussfilt(cevkd - aevkd, 2), clim)
    line([1 1] * ts(tpts(2)), get(gca, 'YLim'), 'LineStyle', '--')
    set(gca, 'xticklabel', xtls, 'yticklabel', ytls, 'ydir', 'normal')

    % COLOR BAR: EVOKED
    pltsz1 = get(gca, 'Position');
    cb = colorbar('EastOutside'); title(cb, 'Evoked (dB)', 'FontSize', 8)
    colormap(eval(clrmp{2}))
    pltsz2 = get(gca, 'Position'); % adjust subplot parent
    set(gca, 'Position', [pltsz2(1:2) pltsz1(3) pltsz2(4)])

    
    % BOTTOM ROW: ERPs
    ylm = [-2 1];
    subplot(rws, cls, 10)
    line([tpts(2) tpts(2)], ylm); line([1 numel(ts)], [0 0]); hold on
    P = plot(Cerp');
    set(P, {'LineStyle'}, {'--', ':', '-'}', 'Color', 'k', 'LineWidth', 1.2)
    set(gca, 'YLim', ylm, 'xtick', tpts, 'xticklabel', {'-.5' '0' '.5'})
%     line([tpts(2) tpts(2)], ylm, 'LineStyle', '--')
    
    ylabel(sprintf('ERPs\nmean power (dB)'), 'FontWeight', 'bold')

    subplot(rws, cls, 11)
    line([tpts(2) tpts(2)], ylm); line([1 numel(ts)], [0 0]); hold on
    P = plot(Aerp');
    set(P, {'LineStyle'}, {'--', ':', '-'}', 'Color', 'k', 'LineWidth', 1.2)
    set(gca, 'YLim', ylm, 'xtick', tpts, 'xticklabel', {'-.5' '0' '.5'}, 'yticklabel', ytls)
%     line([tpts(2) tpts(2)], ylm, 'LineStyle', '--')
    
    text(0, -3, myReport([ROI{2, r} ' ROI: ' biosemi1020(roi)]), 'Color', 'r', 'FontSize', 8)

    subplot(rws, cls, 12)
    line([tpts(2) tpts(2)], ylm); line([1 numel(ts)], [0 0]); hold on
    P = plot(CvAerp');
    set(P, {'LineStyle'}, {'--', ':', '-'}', 'Color', 'k', 'LineWidth', 1.2)
    set(gca, 'YLim', ylm, 'xtick', tpts, 'xticklabel', {'-.5' '0' '.5'}, 'yticklabel', ytls)

    % Legend
    lg = legend(P, {'ERSP' 'Induced', 'Evoked/10'}...
               , 'Location', 'NorthEast'...
               , 'FontSize', 8 ...
               , 'Box', 'off');
    lg.Position = [lg.Position(1)+lg.Position(3) lg.Position(2:4)];
    lg.ItemTokenSize = [10 3];
    
    % Label for Time axes
    xlabel('Time (s)'...
        , 'FontWeight', 'bold'...
        , 'Units', 'normalized')
    
    
    %% PRINT FIGURE
    print(figh, '-dsvg', fullfile(oud, 'ersp', svnm))
    close
    myReport(['Processed ' ROI{2, r} ' trodes: ' biosemi1020(roi)]);
    
    %% CLEAR WORKSPACE FOR NEXT DATASET
    clear t roi *erp *boot clim* de_mask ts xts ttl fh figh svnm alltf *ersp *evk*...
        *idu* freqs times plt* rws cls *tls itc blpow
end %% End of ROI-select loop ----


function eeg = rmgrandbl(eeg, dims, poststim)

    % The Makeig 'local' method, remove each epoch's own baseline mean
    if isnan(dims)
        eeg.data = rmbase(eeg.data, eeg.pnts, 1:(eeg.pnts / 2));
    else
        BL = mean(eeg.data(:, 1:eeg.pnts/2, :), dims);
        poststim = squeeze(mean(eeg.data(:, 1:eeg.pnts/2, :), 2)) * poststim;

        if isscalar(dims) %treat each channel separately
            for e = 1:eeg.trials
                eBL = [BL repmat(poststim(:, e), 1, eeg.pnts/2)];
                eeg.data(:, :, e) = eeg.data(:, :, e) - eBL;
            end
        else %use a BL derived from average of all channels - WHY?
    %         BL = repmat([BL poststim * ones(1, eeg.pnts/2)], eeg.nbchan, 1, eeg.trials);
        end
    end
end