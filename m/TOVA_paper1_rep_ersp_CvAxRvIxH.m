cond = {'response' 'inhibition'};

alpha = [0.05 0.0005];
naccu = [200 2000];
ntest = numel(alpha);
bl = NaN;
figh = figure('Position', get(0, 'ScreenSize'), 'Color', 'w');
% figh = figure('Position', [1 1 lbwh(4) lbwh(3)], 'Color', 'w');%, 'Visible', false);
svnm = 'TOVA_ersp-TEST_makeig-BL_cyc3-11';
rws = 6;
cls = 4;
plt = 0;
clim = [-2 2];


%% ERSPs
%TODO: test other than morlet wavelets?
%TODO - PLAY WITH PARAMS??
%                 , 'detrend', 'on'...
%                 , 'rmerp', 'on'...
for c = which_c
for h = which_h

    s = cellfun(@(x) [svnm '_H' h '_' cond{c} '_' x '.mat'], ROI(2, :), 'Un', 0);
    if ~all(isfile(fullfile(oud, 'ersp', s)))
        eeg = eval(['cEEG' cnd{c} '_H' h]);
        cEEG = pop_mergeset(eeg, 1:numel(eeg));
        cEEG = rmgrandbl(cEEG, NaN, 0);
        eeg = eval(['aEEG' cnd{c} '_H' h]);
        aEEG = pop_mergeset(eeg, 1:numel(eeg));
        aEEG = rmgrandbl(aEEG, NaN, 0);
    end
    
    
    %% Per ROI
    for r = which_r
        roi = ROI{1, r}(4);
        myReport(['Processing ' ROI{2, r} ' trodes: ' biosemi1020(roi)]);
        savename = [svnm '_H' h '_' cond{c} '_' ROI{2, r}];
        if ~isfile(fullfile(oud, 'ersp', [savename '.mat']))
            load(fullfile(oud, 'ersp', [savename '.mat']))
        else
            % get ERSP and ITC matrices from each ROI electrode
            ersp = cell(numel(roi), 1);
            itc = cell(numel(roi), 1);
            blpow = cell(numel(roi), 1);
            eboot = cell(numel(roi), 1);
            iboot = cell(numel(roi), 1);
            alltf = cell(numel(roi), 1);

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

        end


            %% CONDSTAT! - see ERSP v1 if you want to get ITC vals...
            % -----------------
            % get power values
            tfC = mean(structarr2mat(alltf, 1), 4);
            tfCpwr = tfC .* conj(tfC);
            % remove baseline if data was baselined
            blC = (10 .^ (mean(structarr2mat(blpow, 1), 3) / 20))';
            if ~isnan(blC)
                tfCpwr = tfCpwr ./ repmat(blC, [1 size(tfCpwr, 2) size(tfCpwr, 3)]);
            end

            % get power values
            tfA = mean(structarr2mat(alltf, 2), 4);
            tfApwr = tfA .* conj(tfA);
            % remove baseline if data was baselined
            blA = (10 .^ (mean(structarr2mat(blpow, 2), 3) / 20))';
            if ~isnan(blA)
                tfApwr = tfApwr ./ repmat(blA, [1 size(tfApwr, 2) size(tfApwr, 3)]);
            end

            % Get condition differences and their bootstrapped 95% CIs
            CAdiff = cell(1, ntest);
            CA95CIs = cell(1, ntest);
            for a = 1:ntest
                [CAdiff{a}, CA95CIs{a}, ~, ~] = condstat(...
                    'log10(mean(arg1(:,:,X),3))'...
                    , naccu(a)...
                    , alpha(a)...
                    , 'both'...
                    , 'abs'...
                    , {tfCpwr tfApwr});
            end
            clear a tfC tfA blC blA tfCpwr tfApwr


            %% GET FDR CORRECTED MASK! - for the DIFFERENCE data 
            % see ERSP v1 if you want masks for CTRL, ADHD...
            % -----------------
            de_mask = cell(2, ntest);
            for d = 1:ntest
                de_mask{1, d} =...
                    CAdiff{d} < CA95CIs{d}(:,:,1) | CAdiff{d} > CA95CIs{d}(:,:,2);
                if d < ntest
                    de_mask{2, d} = edge(de_mask{1, d}, 'log');
                end
            end
            de_mask{2, ntest} = sum(cat(3, de_mask{2, :}, de_mask{1, ntest}), 3);
            clear d
        
            %% SAVE THE CALCULATIONS
            save(fullfile(oud, 'ersp', [savename '.mat']), 'ersp', 'itc' ,'blpow'...
                , 'times', 'freqs', 'eboot', 'iboot', 'CAdiff', 'CA95CIs', 'de_mask'...
                , '-v7.3')
        end

        
        %% Average ERSP & ITC matrices from all ROI 'trodes
        cersp = mean(structarr2mat(ersp, 1), 3);
        aersp = mean(structarr2mat(ersp, 2), 3);
        
        
        %% Plot the ERSP-EVOKED-INDUCED figure
        pltoff = (r - 1) + 12 ^ (r - 1);
        
        % init time & title
        ts = times / 1000;
        [~, zeropt] = min(abs(ts));
        xts = (-500:100:500) / 1000;
        ttl = Grup;%cellfun(@(x) [x [', H' h cond{c} ', ' ROI{2, r}]], Grup, 'Un', 0);

        % TOP ROW : Control - ADHD
%         % get mirrored color map limit values
%         clim = [-1 1] * max(abs([min(cat(3, c1, a1), [], 'all') ...
%                                  max(cat(3, c1, a1), [], 'all')]));
        subplot(rws, cls, plt + pltoff)
        imagesclogy(ts, freqs, imgaussfilt(cersp, 2), clim)
        if plt + pltoff == 1
            ytls = get(gca, 'YTickLabels');
            ylabel(sprintf('%s\n Frequency (Hz)', ttl{1}), 'FontWeight', 'bold')
            title(['H' h ', ' cond{c}])
        elseif any(plt + pltoff == 2:4)
            title(['H' h ', ' cond{c}])
            ytls = '';
        elseif plt + pltoff == 13
            ytls = get(gca, 'YTickLabels');
            ylabel(sprintf('%s\n', ttl{1}), 'FontWeight', 'bold')
        else
            ytls = '';
        end
        set(gca, 'xtick', xts, 'xticklabel', '', 'yticklabel', ytls, 'ydir', 'normal')
        line([1 1] * ts(zeropt), get(gca, 'YLim'), 'LineStyle', '--')
        
        % Color bar and adjustment of it's subplot parent
        if plt + pltoff == 4
            pltsz1 = get(gca, 'Position');
            cb = colorbar('EastOutside'); title(cb, 'ERSP (dB)', 'FontSize', 9)
            colormap(eval(clrmp{2}))
            pltsz2 = get(gca, 'Position');
            set(gca, 'Position', [pltsz2(1:2) pltsz1(3) pltsz2(4)])
        end

        subplot(rws, cls, plt + pltoff + 4)
        imagesclogy(ts, freqs, imgaussfilt(aersp, 2), clim)
        set(gca, 'xtick', xts, 'xticklabel', '', 'yticklabel', ytls, 'ydir', 'normal')
        line([1 1] * ts(zeropt), get(gca, 'YLim'), 'LineStyle', '--')
        if any(plt + pltoff + 4 == [5 17])
            ylabel(sprintf('-- %s ROI --\n%s\n ', ROI{2, r}, ttl{2}), 'FontWeight', 'bold')
        end

        subplot(rws, cls, plt + pltoff + 8)
        fh = ctap_imagesclogy(ts, freqs, 10 * CAdiff{2} .* de_mask{1, 1}...
            , 'clim', clim);
        set(fh, 'AlphaData', ~(de_mask{1, 1} - de_mask{1, 2}) + 0.2)
        line([1 1] * ts(zeropt), get(gca, 'YLim'), 'LineStyle', '--')
        if any(plt + pltoff + 8 == 21:24)
            set(gca, 'xtick', xts, 'yticklabel', ytls, 'ydir', 'normal')
            xtickangle(45)
        else
            set(gca, 'xtick', xts, 'xticklabel', '', 'yticklabel', ytls, 'ydir', 'normal')
        end
        if any(plt + pltoff + 8 == [9 21])
            ylabel(sprintf('p(%s - %s)\n', Grup{1}, Grup{2}), 'FontWeight', 'bold')
        end

        % Label for Time axes
        if plt + pltoff + 8 == 24
            xlabel('Time (s)'...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [1.1 -0.05]...
                , 'Rotation', 45)
        end

        myReport(['Processed ' savename ' trodes: ' biosemi1020(roi)]);
    end %% End of ROI-select loop ----
    
    plt = plt + 1;
    
end %% End of TOVA-halves loop ----
end %% End of conditions loop ----

print(figh, '-dsvg', fullfile(oud, 'ersp', [svnm '_v4']))

clear t roi CA* clim* *ersp de_mask ts xts ttl fh svnm cls rws blpow *tf itc ...
    ntest naccu alpha *boot plt* ytls times freqs s r h c cb bl


%% SUB FUNCTIONS
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