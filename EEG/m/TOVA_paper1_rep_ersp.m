ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extended_anal';
ROI = {[7 15 19 21 23 28 36] [7 19 21 36] [15 23 28] [81 83 85 87] [68 85 100];
        'parieto-occipital' 'parietal' 'occipital' 'frontocentral' 'frontal'};
group = {'Control' 'ADHD'};
cond = {'rsp' 'inb'; 'response' 'inhibition'};
lbwh = get(0,'ScreenSize');
clrmp = {'whatjet(''what'', [0.9 0.9 0.9], ''stops'', [0 0.15 0.25 0.1 0.1 0.25 0.15])'...
         'whatjet(''what'', [0.95 0.95 0.95])'...
         'viridis';
         'greytight' 'whiter' 'viridis'};


%% ERPIM x pre-stim alpha phase
cEEG = pop_mergeset( ...
    pop_loadset('filepath', ind, 'filename', 'MERGED_H1_CORRECT_RANDOMISED_CONTROL.set')...
  , pop_loadset('filepath', ind, 'filename', 'MERGED_H2_CORRECT_RANDOMISED_CONTROL.set'));
ix = cellfun(@(x) any(ismember(x, 'cor_rsp')), {cEEG.epoch.eventtype});
cEEGrsp = pop_select(cEEG, 'trial', find(ix));
ix = cellfun(@(x) any(ismember(x, 'cor_inb')), {cEEG.epoch.eventtype});
cEEGinb = pop_select(cEEG, 'trial', find(ix));

aEEG = pop_mergeset( ...
    pop_loadset('filepath', ind, 'filename', 'merged_h1_correct_randomised_adhd.set')...
  , pop_loadset('filepath', ind, 'filename', 'merged_h2_correct_randomised_adhd.set'));
ix = cellfun(@(x) any(ismember(x, 'cor_rsp')), {aEEG.epoch.eventtype});
aEEGrsp = pop_select(aEEG, 'trial', find(ix));
ix = cellfun(@(x) any(ismember(x, 'cor_inb')), {aEEG.epoch.eventtype});
aEEGinb = pop_select(aEEG, 'trial', find(ix));


%% Topoplot the neg+pos peaks
cax = {[-2 0] 'maxmin' 'absmax'};
% itv = 565:616;
m = 1; % using greytight since it works for ERSPs
roi = [ROI{1, 1} ROI{1, 5}];

for s = 1%:numel(cax) %using only matched scales since they seem best
    t0 = 513;
    tinc = 150;
    
    for t = tinc:tinc:500 % five time slices
        fh = figure('Position', [1 1 lbwh(4) lbwh(3)]);
        plt = 1;
        t1 = find(cEEG.times > t, 1);
        sgt = sprintf('%d-%dms', t-tinc, t);
        
        for c = size(cond, 2):-1:1 %conditions across subplot rows
            cnd = ['EEG' cond{1, c}];
            if s == 1
                dat = zeros(2, 128);
                for g = 1:2
                    grp = lower(group{g}(1));
                    dat(g, :) = ...
                        mean(mean(eval([grp cnd]).data(:, t0:t1, :), 3), 2);
                end
                scl = [min(dat, [], 'all') max(dat, [], 'all')];
            else
                scl = cax{s};
            end
            for g = 1:2 %groups across columns
                grp = lower(group{g}(1));

                subplot(2, 2, plt)
                dat = mean(mean(eval([grp cnd]).data(:, t0:t1, :), 3), 2);
                chlocs = eval([grp cnd]).chanlocs;
                topoplot(dat, chlocs...
                    , 'electrodes', 'on'...
                    , 'headrad', 0 ...
                    , 'plotrad', 0.6 ...
                    , 'emarker', {'.', 'k', 2, 1}...
                    , 'emarker2', {roi, 'o', 'k', 7, 1}...
                    , 'maplimits', scl...
                    , 'colormap', eval(clrmp{1, m}))
                hold on
                topoplot([], chlocs...
                    , 'style', 'blank'...
                    , 'headrad', 0 ...
                    , 'plotrad', 0.6 ...
                    , 'plotchans', roi(1:7)...
                    , 'electrodes', 'on'...
                    , 'emarker', {'o', 'm', 5, 2}...
                    , 'colormap', eval(clrmp{1, m}))
                topoplot([], chlocs...
                    , 'style', 'blank'...
                    , 'headrad', 0.5 ...
                    , 'plotrad', 0.6 ...
                    , 'plotchans', roi(8:10)...
                    , 'electrodes', 'on'...
                    , 'emarker', {'x', 'm', 6, 2}...
                    , 'colormap', eval(clrmp{1, m}))
                title(sprintf('%s-%s:%s', group{g}, cond{2, c}, sgt))
                cb = colorbar('EastOutside'); title(cb, '\muV', 'FontSize', 9)

                plt = plt + 1;
            end
%             if s == 1, scl = [-3 0.5]; end
        end
        t0 = t1 + 1;
        
        set(fh, 'Color', 'w')
        sgtitle(['Scalp maps at ' sgt ' post-stimulus'], 'Color', 'r')
        if s == 1, scl = 'match'; end
        svnm = ['TOVA_topoplot-' sgt '_' scl 'Scales_' clrmp{2, m}];
%         print(fh, '-dpng', fullfile(oud, 'topos', svnm))
        export_fig(fullfile(oud, 'topos', svnm), '-png')
        close
    end
end
clear fh grp cnd cax plt itv s m c g dat cb

    % inc = num2str(round((1000 / cEEGrt.srate) * inc - 2));
    % text(-3, -1, ['Scalp maps at ' inc 'ms intervals around RT'], ...
    %   'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')


%% ERSPs
alpha = [0.05 0.0005];
naccu = [200 2000];
ntest = numel(alpha);
cyc = {0 [3 0.5]};
meth = {'FFT' 'wavelet3-05'};
inc = 2:0.5:4;
% 17 channels spread around the scalp
chs = [94 72 115 54 1 100 68 81 85 23 7 36 10 39 19 119 58];
biosemi = {'C30' 'C8' 'D19' 'B22' 'A1' 'D4' 'C4' 'C17' 'C21'...
            'A23' 'A7' 'B4' 'A10' 'B7' 'A19' 'D23' 'B26'};
ten20 =   {'AF7' 'AF8' 'C3' 'C4' 'Cz' 'F3' 'F4' 'Fpz' 'Fz'...
            'Oz' 'P3' 'P4' 'PO7' 'PO8' 'Pz' 'T7' 'T8'};


%% ERSP - Pz
% chs = 19; % debug 
for ch = 1:numel(chs)
%     for l = 1:numel(cyc)

        for c = 1:size(cond, 2)
            for i = 1:numel(inc)
                cnd = ['EEG' cond{1, c}];
%                 [ersp,itc,powbase,times,freqs,erspboot,itcboot]  =  ...
                newtimef(...
                    {eval(['c' cnd]).data(chs(ch),:,:) eval(['a' cnd]).data(chs(ch),:,:)}...
                    , cEEG.pnts...
                    , [cEEG.times(1) cEEG.times(end)]...
                    , cEEG.srate...
                    , cyc{l} ...
                    , 'freqs', [inc(i) 30] ...
                    , 'nfreqs', 27 ...
                    , 'freqscale', 'log'...
                    , 'alpha', 0.05...
                    , 'baseline', [-1000 -1]...
                    , 'commonbase', 'on'...
                    , 'trialbase', 'full'...
                    , 'title', cellfun(@(x) ...
                        [x [' - ' cond{2, c} ' - ' ten20{ch}]], group, 'Un', 0));

                %TODO - PLAY WITH PARAMS??
                % 'detrend', 'on'
                % 'rmerp', 'on'

                fh = gcf;
                fh.Position = lbwh;
                if ~isfolder(fullfile(oud, 'ersp', ten20{ch}))
                    mkdir(fullfile(oud, 'ersp', ten20{ch}))
                end
                svnm = ['TOVA_ersp_' cond{2, c} '_' ten20{ch} '_' meth{l}...
                        '_logfreq' strrep(num2str(inc(i)), '.', 'pt') '-30'];
                print(fh, '-dsvg', fullfile(oud, 'ersp', ten20{ch}, svnm))
                close
            end
        end
%     end
end
clear ch chs biosemi ten20 c i cnd cyc inc svnm fh


%% ERSP - ROIs


%% FFT
% for c = 1:size(cond, 1)
%     % get ERSP and ITC matrices from each TRD electrode
%     for t = 1:numel(ROI)
%         cnd = ['EEG' cond{1, c}];
%         [ersp{t}, itc{t}, blpow{t}, times, freqs, eboot{t}, iboot{t}, alltf{t}] = ...
%         newtimef(...
%             {eval(['c' cnd]).data(ROI(t),:,:) eval(['a' cnd]).data(ROI(t),:,:)}...
%             , cEEG.pnts...
%             , [cEEG.times(1) cEEG.times(end)]...
%             , cEEG.srate...
%             , 0 ...
%             , 'plotersp', 'off'...
%             , 'plotitc', 'off'...
%             , 'freqs', [2 30] ...
%             , 'nfreqs', 27 ...
%             , 'freqscale', 'log'...
%             , 'alpha', 0.05...
%             , 'baseline', [-1000 -1]...
%             , 'commonbase', 'on'...
%             , 'trialbase', 'full'...
%             , 'title'...
%             , cellfun(@(x) [x [' - ' cond{2, c} ' - ROI']], group, 'Un', 0));
% 
%         %TODO - PLAY WITH PARAMS??
%         % 'detrend', 'on'
%         % 'rmerp', 'on'
%     end
%     % Average ERSP & ITC matrices from all TRD 'trodes, composite as fig
% 
%     fh = gcf;
%     fh.Position = lbwh;
%     print(fh, '-dsvg', fullfile(oud, 'ersp', ['TOVA_ersp_' cond{2, c}...
%                 '_ROI_FFT_logfreq2-30']))
%     close
% end

%% morlet wavelets

% c = 1; i = 5; % debug

for r = 1:size(ROI, 2)
    roi = ROI{1, r};
    % get ERSP and ITC matrices from each ROI electrode
    ersp = cell(numel(roi), 1);
    itc = cell(numel(roi), 1);
    blpow = cell(numel(roi), 1);
    eboot = cell(numel(roi), 1);
    iboot = cell(numel(roi), 1);
    alltf = cell(numel(roi), 1);
    
    for c = 1:size(cond, 2)
        for t = 1:numel(roi)
            cnd = ['EEG' cond{1, c}];
            [ersp{t}, itc{t}, blpow{t}, times, freqs, eboot{t}, iboot{t}, alltf{t}] = ...
            newtimef(...
                {eval(['c' cnd]).data(roi(t),:,:) eval(['a' cnd]).data(roi(t),:,:)}...
                , cEEG.pnts...
                , [cEEG.times(1) cEEG.times(end)]...
                , cEEG.srate...
                , [3 0.5] ...
                , 'plotersp', 'off'...
                , 'plotitc', 'off'...
                , 'freqs', [4 30] ...
                , 'nfreqs', 54 ...
                , 'freqscale', 'log'...
                , 'baseline', NaN ... % 0 removes all negative values
                , 'commonbase', 'on'...
                , 'trialbase', 'full');

            %TODO - PLAY WITH PARAMS??
%                 , 'alpha', 0.05...
%                 , 'detrend', 'on'...
%                 , 'rmerp', 'on'...
%                 , 'pcontour', 'on'...
%                 , 'mcorrect', 'fdr'...
%                 , 'plotersp', 'on'...
%                 , 'plotitc', 'on'...

        end

        
%% Average ERSP & ITC matrices from all ROI 'trodes
        ceroi = mean(structarr2mat(ersp, 1), 3);
        aeroi = mean(structarr2mat(ersp, 2), 3);
        
        % try to replicate the ERSP-process in newtimef but outcome is bad
%         ceroi = mean(structarr2mat(alltf, 1), 4);
%         ceroi = ceroi .* conj(ceroi);
%         ceroi = log10(mean(ceroi, 3));
%         
%         aeroi = mean(structarr2mat(alltf, 2), 4);
%         aeroi = aeroi .* conj(aeroi);
%         aeroi = log10(mean(aeroi, 3));
        
        % get color map limit values
        erspmax = max(cat(3, ceroi, aeroi), [], 'all');
        erspmin = min(cat(3, ceroi, aeroi), [], 'all');
        clim = [-1 1] * max(abs([erspmin erspmax])); %mirrored colormap limits
        climnmx = [erspmin erspmax]; %colormap limits from range of data
        clear erspmin erspmax
        
        
        %% CONDSTAT!
        % get power values
        tfC = mean(structarr2mat(alltf, 1), 4);
        tfCpwr = tfC .* conj(tfC);
        % remove baseline if data was baselined
        blC = (10 .^ (mean(structarr2mat(blpow, 1), 3) / 20))';
        if ~isnan(blC)
            tfCpwr = tfCpwr ./ repmat(blC, [1 size(tfCpwr, 2) size(tfCpwr, 3)]);
        end
        % get ITC if you want that
%         tfC = tfC ./ repmat(bl1 / 2, [1 size(tfC, 2) size(tfC, 3)]);
%         tfCnorm = tfC ./ sqrt(tfC .* conj(tfC));

        % get power values
        tfA = mean(structarr2mat(alltf, 2), 4);
        tfApwr = tfA .* conj(tfA);
        % remove baseline if data was baselined
        blA = (10 .^ (mean(structarr2mat(blpow, 2), 3) / 20))';
        if ~isnan(blA)
            tfApwr = tfApwr ./ repmat(blA, [1 size(tfApwr, 2) size(tfApwr, 3)]);
        end
        % get ITC if you want that
%         tfA = tfA ./ repmat(bl2 / 2, [1 size(tfA, 2) size(tfA, 3)]);
%         tfAnorm = tfA ./ sqrt(tfA .* conj(tfA));
        
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
    %             , {tfCnorm tfAnorm});
        end
% {'log10(mean(arg1(:,:,X),3))' 'mean(arg2(:,:,X),3)'}... % original formula
        clear tfC tfA blC blA tfCpwr tfApwr


        %% GET FDR CORRECTED MASK!
        % -----------------
        % for the CONTROL data
%         [PbootC, ~, PboottrialsC] = bootstat( tfC ...
%             , 'mean(arg1,3);'...
%             , 'boottype', 'shuffle' ...
%             , 'label', 'ERSP'...
%             , 'bootside', 'both'...
%             , 'naccu', naccu ...
%             , 'basevect', find(times < 0)...
%             , 'alpha', alpha...
%             , 'dimaccu', 2 );
%         
%         % average the power & compute the mask
%         if ismatrix(PboottrialsC), PboottrialsC = PboottrialsC'; end
%         c_exactp_ersp = compute_pvals(ceroi, PboottrialsC);
%         c_alfdr = fdr(c_exactp_ersp, alpha);
%         fprintf('ERSP multiplicity correction by FDR, alpha = %3.6f\n', c_alfdr);
%         c_maskersp = c_exactp_ersp <= c_alfdr;
%         
%         % for the PATIENT data
%         [PbootA, ~, PboottrialsA] = bootstat( tfA ...
%             , 'mean(arg1,3);'...
%             , 'boottype', 'shuffle' ...
%             , 'label', 'ERSP'...
%             , 'bootside', 'both'...
%             , 'naccu', naccu ...
%             , 'basevect', find(times < 0)...
%             , 'alpha', alpha...
%             , 'dimaccu', 2 );
%         
%         % average the power & compute the mask
%         if ismatrix(PboottrialsA), PboottrialsA = PboottrialsA'; end
%         a_exactp_ersp = compute_pvals(ceroi, PboottrialsA);
%         a_alfdr = fdr(a_exactp_ersp, alpha);
%         fprintf('ERSP multiplicity correction by FDR, alpha = %3.6f\n', a_alfdr);
%         a_maskersp = a_exactp_ersp <= a_alfdr;

        % for the DIFFERENCE data
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


        %% Plot the figure
        ts = times / 1000;
        xts = (-500:100:500) / 1000;
        ttl = cellfun(@(x) [x [', ' cond{2, c} ', ROI']], group, 'Un', 0);
        
        for m = 1:size(clrmp, 2)
            if m == 3, clim = climnmx; end
            
            fh = figure('Position', [1 1 lbwh(3) lbwh(4)], 'Color', 'w');
            
            % Top row : actual data
            subplot(2, 3, 1)
            imagesclogy(ts, freqs, imgaussfilt(ceroi, 2), clim)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title(ttl{1})
            xtickangle(45)
            
            ylabel('Frequency (Hz)', 'FontWeight', 'bold')
            
            subplot(2, 3, 2)
            imagesclogy(ts, freqs, imgaussfilt(aeroi, 2), clim)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title(ttl{2})
            xtickangle(45)
            
            subplot(2, 3, 3)
            imagesclogy(ts, freqs, imgaussfilt(ceroi - aeroi, 2), clim)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title([group{1} ' - ' group{2}])
            xtickangle(45)
            
            % Color bar and adjustment of it's subplot parent
            pltsz1 = get(gca, 'Position');
            cb = colorbar('EastOutside'); title(cb, 'ERSP (dB)', 'FontSize', 9)
            colormap(eval(clrmp{1, m}))            
            pltsz2 = get(gca, 'Position');
            set(gca, 'Position', [pltsz2(1:2) pltsz1(3) pltsz2(4)])
            
            % Bottom row : significances
            subplot(2, 3, 4)
            imagesclogy(ts, freqs, 10 * CAdiff{2}, clim)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title(['log(' group{1} ' - ' group{2} ')'])
            xtickangle(45)
            
            subplot(2, 3, 5)
            imagesclogy(ts, freqs, 10 * CAdiff{2} .* de_mask{2, ntest}, clim)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title(['p(' group{1} ' - ' group{2} ')'])
            xtickangle(45)
            
            subplot(2, 3, 6)
            h = my_imagesclogy(ts, freqs, 10 * CAdiff{2} .* de_mask{1, 1}...
                , 'clim', clim);
            set(h, 'AlphaData', ~(de_mask{1, 1} - de_mask{1, 2}) + 0.2)
            set(gca, 'xtick', xts, 'ydir', 'normal')
            title(['p(' group{1} ' - ' group{2} ')'])
            xtickangle(45)
            
            xlabel('Time (s)'...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [1.1 -0.05]...
                , 'Rotation', 45)

            svnm = ['TOVA_ersp_nobl_' cond{2, c} '_ROI-' ROI{2, r}...
                    '_wavelet3-05_logfreq4-30_' clrmp{2, m}];
            print(fh, '-dsvg', fullfile(oud, 'ersp', svnm))
            close
        end
        svnm = ['TOVA_ersp_calcs_nobl_' cond{2, c} '_ROI-' ROI{2, r} '.mat'];
        save(fullfile(oud, 'ersp', svnm), 'ersp', 'itc' ,'blpow', 'times'...
            , 'freqs', 'eboot', 'iboot', 'CAdiff', 'CA95CIs', '-v7.3')

    end %% End of conditions loop ----
    
    clear t cnd roi
end %% End of ROI-select loop ----
    

%% SUB-FUNCTIONS
% -----------
% function pvals = compute_pvals(dat, srg, tail)
%     
%     if nargin < 3
%         tail = 'both';
%     end
%     
%     if myndims(dat) > 1        
%         if size(dat,2) ~= size(srg, 2) || myndims(srg) == 2
%             if size(dat,1) == size(srg, 1)
%                 srg = repmat(reshape(srg, [size(srg,1) 1 size(srg,2)])...
%                     , [1 size(dat,2) 1]);
%             elseif size(dat,2) == size(srg, 1)
%                 srg = repmat(reshape(srg, [1 size(srg,1) size(srg,2)])...
%                     , [size(dat,1) 1 1]);
%             else
%                 error('Permutation statistics array size error');
%             end
%         end
%     end
% 
%     srg = sort(srg, myndims(srg)); % sort last dimension
%     
%     if myndims(srg) == 1    
%         srg(end+1) = dat;        
%     elseif myndims(srg) == 2
%         srg(:,end+1) = dat;        
%     elseif myndims(srg) == 3
%         srg(:,:,end+1) = dat;
%     else
%         srg(:,:,:,end+1) = dat;
%     end
% 
%     [~, idx] = sort( srg, myndims(srg) );
%     [~, mx]  = max( idx,[], myndims(srg));        
%                 
%     len = size(srg,  myndims(srg) );
%     pvals = 1-(mx-0.5)/len;
%     if strcmpi(tail, 'both')
%         pvals = min(pvals, 1-pvals);
%         pvals = 2*pvals;
%     end
% end
%     
% function val = myndims(a)
%     if iscolumn(a)
%         val = 1;
%     else
%         val = ndims(a);
%     end
% end


%% TO AVERAGE CHANNELS' TIME SERIES AND GET ERSP, WE CAN JUST DO THIS:
%     for i = 1:numel(inc)
%         newtimef(...
%             {mean(cEEGrsp.data(TRD,:,:)) mean(aEEGrsp.data(TRD,:,:))}...
%             , cEEG.pnts...
%             , [cEEG.times(1) cEEG.times(end)]...
%             , cEEG.srate...
%             , cyc{l} ...
%             , 'freqs', [inc(i) 30] ...
%             , 'nfreqs', 27 ...
%             , 'freqscale', 'log'...
%             , 'alpha', 0.05...
%             , 'baseline', [-1000 -1]...
%             , 'commonbase', 'on'...
%             , 'trialbase', 'full'...
%             , 'title', cellfun(@(x) [x [' - ' cond '- ' NM]], group, 'Un', 0))
% 
%         fh = gcf;
%         fh.Position = lbwh;
%         print(fh, '-dsvg', fullfile(oud, ersp, svnm))
%         close
%     end