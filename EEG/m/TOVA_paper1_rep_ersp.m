ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extended_anal';
ROI = [7 15 19 21 23 28 36];
Pz = 19;
group = {'Control' 'ADHD'};
condition = {'inhibition' 'response'};
lbwh = get(0,'ScreenSize');
cmp = {'jet' 'whitejet' 'viridis'};


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
grp = {'c' 'a'};
cnd = {'inb' 'rsp'};

cax = {[-2 0] 'maxmin' 'absmax'};
itv = 565:616;

for s = 1:numel(cax)
    for m = 1:numel(cmp)
        
        fh = figure('Position', [1 1 lbwh(4) lbwh(3)], 'Color', 'w');
        plt = 1;
        scl = cax{s};
        for c = 1:2 %conditions across rows
            for g = 1:2 %groups across columns

                subplot(2, 2, plt)
                dat = mean(mean(eval([grp{g} 'EEG' cnd{c}]).data(:, itv, :), 3), 2);
                topoplot(dat, eval([grp{g} 'EEG' cnd{c}]).chanlocs...
                    , 'electrodes', 'on'...
                    , 'emarker', {'.','k',[],1}...
                    , 'emarker2', {ROI, 'o','m',3,1}...
                    , 'maplimits', scl...
                    , 'colormap', eval(cmp{m}))
                title(sprintf('%s-%s', group{g}, condition{c}))
                cb = colorbar('EastOutside'); title(cb, '\muV', 'FontSize', 9)

                plt = plt + 1;
            end
            if s == 1, scl = [-3 0.5]; end
        end
        % inc = num2str(round((1000 / cEEGrt.srate) * inc - 2));
        sgt = sgtitle('Scalp maps at 100-200ms post-stimulus', 'Color', 'r');
        % text(-3, -1, ['Scalp maps at ' inc 'ms intervals around RT'], ...
        %   'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
        % colormap(whitejet)
        % export_fig(fullfile(oud, 'TOVA-topoplots-ERPxRT.pdf'))
        if s == 1, scl = 'match'; end
        print(fh, '-dsvg', fullfile(oud, 'topos'...
            , ['TOVA-topoplots-ERSP-' scl 'Scales_' cmp{m}]))
        close
    end
end
clear fh grp cnd cax plt itv s m c g dat cb


%% ERSPs
alpha = 0.05;
cyc = {0 [3 0.5]};
meth = {'FFT' 'wavelet3-05'};
cond = {'rsp' 'inb'; 'response' 'inhibition'};
inc = 2:0.5:4;
% 17 channels spread around the scalp
chs = [94 72 115 54 1 100 68 81 85 23 7 36 10 39 19 119 58];
biosemi = {'C30' 'C8' 'D19' 'B22' 'A1' 'D4' 'C4' 'C17' 'C21' 'A23' 'A7' 'B4' 'A10' 'B7' 'A19' 'D23' 'B26'};
ten20 =   {'AF7' 'AF8' 'C3' 'C4' 'Cz' 'F3' 'F4' 'Fpz' 'Fz' 'Oz' 'P3' 'P4' 'PO7' 'PO8' 'Pz' 'T7' 'T8'};


%% ERSP - Pz
for ch = 1:numel(chs)
%     for l = 1:numel(cyc)

        for c = 1:size(cond, 1)
            for i = 1:numel(inc)
                cnd = ['EEG' cond{1, c}];
                [ersp,itc,powbase,times,freqs,erspboot,itcboot]  =  ...
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
                    , 'title'...
                    , cellfun(@(x) [x [' - ' cond{2, c} ' - ' ten20{ch}]], group, 'Un', 0));

                %TODO - PLAY WITH PARAMS??
                % 'detrend', 'on'
                % 'rmerp', 'on'

                %TODO - scalpmap?
                % which('chanlocs128_pist.elp')

                fh = gcf;
                fh.Position = lbwh;
                if ~isfolder(fullfile(oud, 'ersp', ten20{ch}))
                    mkdir(fullfile(oud, 'ersp', ten20{ch}))
                end
                print(fh, '-dsvg', fullfile(oud, 'ersp', ten20{ch}...
                    , ['TOVA_ersp_' cond{2, c} '_' ten20{ch} '_' meth{l}...
                    '_logfreq' strrep(num2str(inc(i)), '.', 'pt') '-30']))
                close
            end
        end
%     end
end


%% ERSP - ROIs
ROI = [15 19 21];
% get ERSP and ITC matrices from each ROI electrode
ersp = cell(numel(ROI), 1);
itc = cell(numel(ROI), 1);
powbase = cell(numel(ROI), 1);
erspboot = cell(numel(ROI), 1);
itcboot = cell(numel(ROI), 1);


%% FFT
for c = 1:size(cond, 1)
    % get ERSP and ITC matrices from each TRD electrode
    for t = 1:numel(ROI)
        cnd = ['EEG' cond{1, c}];
        [ersp{t}, itc{t}, powbase{t}, times{t}, freqs{t}, erspboot{t}, itcboot{t}]  =  ...
        newtimef(...
            {eval(['c' cnd]).data(ROI(t),:,:) eval(['a' cnd]).data(ROI(t),:,:)}...
            , cEEG.pnts...
            , [cEEG.times(1) cEEG.times(end)]...
            , cEEG.srate...
            , 0 ...
            , 'plotersp', 'off'...
            , 'plotitc', 'off'...
            , 'freqs', [2 30] ...
            , 'nfreqs', 27 ...
            , 'freqscale', 'log'...
            , 'alpha', 0.05...
            , 'baseline', [-1000 -1]...
            , 'commonbase', 'on'...
            , 'trialbase', 'full'...
            , 'title'...
            , cellfun(@(x) [x [' - ' cond{2, c} ' - ROI']], group, 'Un', 0));

        %TODO - PLAY WITH PARAMS??
        % 'detrend', 'on'
        % 'rmerp', 'on'

        %TODO - scalpmap?
        % which('chanlocs128_pist.elp')
    end
    % Average ERSP & ITC matrices from all TRD 'trodes, composite as fig

    fh = gcf;
    fh.Position = lbwh;
    print(fh, '-dsvg', fullfile(oud, 'ersp', ['TOVA_ersp_' cond{2, c}...
                '_ROI_FFT_logfreq2-30']))
    close
end

%% morlet wavelets

% for c = 1:size(cond, 1)
%     for i = 1:numel(inc)
        
        for t = 1:numel(ROI)
            cnd = ['EEG' cond{1, c}];
            [ersp{t}, itc{t}, powbase{t}, times, freqs, erspboot{t}, itcboot{t}]  =  ...
            newtimef(...
                {eval(['c' cnd]).data(ROI(t),:,:) eval(['a' cnd]).data(ROI(t),:,:)}...
                , cEEG.pnts...
                , [cEEG.times(1) cEEG.times(end)]...
                , cEEG.srate...
                , [3 0.5] ...
                , 'plotersp', 'off'...
                , 'plotitc', 'off'...
                , 'alpha', 0.05...
                , 'freqs', [inc(i) 30] ...
                , 'nfreqs', 27 ...
                , 'freqscale', 'log'...
                , 'baseline', [-1000 -1]...
                , 'commonbase', 'on'...
                , 'trialbase', 'full');

            %TODO - PLAY WITH PARAMS??
            % 'detrend', 'on'
            % 'rmerp', 'on'

            %TODO - scalpmap?
            % which('chanlocs128_pist.elp')
        end
        
        % Average ERSP & ITC matrices from all ROI 'trodes, composite as fig
        ceroi = sbf_get_avg_mat(ersp, 1);
        aeroi = sbf_get_avg_mat(ersp, 2);
        P = {ceroi aeroi ceroi - aeroi};
        erspmax = max(max(max(cat(3, ceroi, aeroi))));
        erspmin = min(min(min(cat(3, ceroi, aeroi))));
        clim = [erspmin erspmax];
        
%         ciroi = sbf_get_avg_mat(itc, 1);
%         airoi = sbf_get_avg_mat(itc, 2);
%         R = {ciroi airoi ciroi - airoi};


        %% GET FDR CORRECTED MASK!
        % -----------------
        [Pboot, Pboottrialstmp, Pboottrials] = bootstat( P ...
            , 'mean(arg1,3);'...
            , 'boottype', 'shuffle' ...
            , 'label', 'ERSP'...
            , 'bootside', 'both'...
            , 'naccu', 200 ...
            , 'basevect', find(times < 0)...
            , 'alpha', alpha...
            , 'dimaccu', 2 );
        
        % average the power & compute the mask
        if ndims(P) == 4,     Pavg = mean(P, 4);
        elseif ndims(P) == 3, Pavg = mean(P, 3);
        end
        if ismatrix(Pboottrials)
            Pboottrials = Pboottrials'; 
        end
        exactp_ersp = compute_pvals(Pavg, Pboottrials);
        alphafdr = fdr(exactp_ersp, alpha);
        if alphafdr ~= 0
            fprintf('ERSP multiplicity correction by FDR, alpha_fdr = %3.6f\n'...
                , alphafdr);
        else
            fprintf('ERSP multiplicity correction by FDR, not significant\n');
        end
        maskersp = exactp_ersp <= alphafdr;
        

        %% Plot the figure
        fh = figure();
        subplot(1, 4, 1)
        imagesclogy(times, freqs, ceroi, clim)
        set(gca, 'ydir', 'normal')
        subplot(1, 4, 2)
        imagesclogy(times, freqs, aeroi, clim)
        set(gca, 'ydir', 'normal')
        subplot(1, 4, 3)
        imagesclogy(times, freqs, ceroi - aeroi, clim)
        set(gca, 'ydir', 'normal')
        cbar
        subplot(1, 4, 4)
        imagesclogy(times, freqs, maskersp, clim)
        set(gca, 'ydir', 'normal')
        
        colormap(eval(cmp{2}))
        
%         inputdata = alltfX;
%         switch g.type
%             case 'coher'
%                 s = int2str(size(alltfX,3));
%                 f = [ 'sum(arg1,3)./sqrt(sum(arg1.*conj(arg1),3))/ sqrt(' s ');' ];
%             case 'phasecoher'
%                 f =  'mean(arg1,3);'; 
%                 inputdata = alltfX./sqrt(alltfX.*conj(alltfX));
%             case 'phasecoher2'
%                 f = 'sum(arg1,3)./sum(sqrt(arg1.*conj(arg1)),3);';
%         end
%         [Rboot, Rboottmp, Rboottrials] = bootstat(inputdata...
%             , formula...
%             , 'boottype', 'shuffle' ...
%             , 'label', 'ITC'...
%             , 'bootside', 'upper'...
%             , 'naccu', 200 ...
%             , 'basevect', find(times{1} < 0)...
%             , 'alpha', 0.05...
%             , 'dimaccu', 2 );
% 
%         % Pass average matrices to some EEGLAB function to plot ERSP?
% %         title = cellfun(@(x) [x [' - ' cond{2, c} ' - ROI']], group, 'Un', 0);
% 
%         if ndims(P) == 3
%             P = squeeze(P(2,:,:,:));
%             R = squeeze(R(2,:,:,:));
%             mbase = squeeze(mbase(2,:));
%             ERP = mean(squeeze(data(1,:,:)),2);
%         else      
%             ERP = mean(data,2);
%         end
%         if strcmpi(g.plottype, 'image')
%             plottimef(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, maskersp, maskitc, g);
%         else
%             plotallcurves(P, R, Pboot, Rboot, ERP, freqs, timesout, mbase, g);
%         end
%         
%         fh = gcf;
%         fh.Position = lbwh;
%         print(fh, '-dsvg', fullfile(oud, 'ersp', ['TOVA_ersp_' cond{2, c}...
%             '_ROI_wavelet3-05_logfreq' strrep(num2str(inc(i)), '.', 'pt')...
%             '-30_' cmp{3}]))
%         close
%     end
% end

function out = sbf_get_avg_mat(mats, dim)
    cel2mat = cellfun(@(x) x{dim}, mats, 'UniformOutput', false);
    cel2mat = cat(3, cel2mat{:});
    out = mean(cel2mat, 3);

end

% -----------
function pvals = compute_pvals(dat, srg, tail)
    
    if nargin < 3
        tail = 'both';
    end
    
    if myndims(dat) > 1        
        if size(dat,2) ~= size(srg, 2) || myndims(srg) == 2
            if size(dat,1) == size(srg, 1)
                srg = repmat(reshape(srg, [size(srg,1) 1 size(srg,2)])...
                    , [1 size(dat,2) 1]);
            elseif size(dat,2) == size(srg, 1)
                srg = repmat(reshape(srg, [1 size(srg,1) size(srg,2)])...
                    , [size(dat,1) 1 1]);
            else
                error('Permutation statistics array size error');
            end
        end
    end

    srg = sort(srg, myndims(srg)); % sort last dimension
    
    if myndims(srg) == 1    
        srg(end+1) = dat;        
    elseif myndims(srg) == 2
        srg(:,end+1) = dat;        
    elseif myndims(srg) == 3
        srg(:,:,end+1) = dat;
    else
        srg(:,:,:,end+1) = dat;
    end

    [~, idx] = sort( srg, myndims(srg) );
    [~, mx]  = max( idx,[], myndims(srg));        
                
    len = size(srg,  myndims(srg) );
    pvals = 1-(mx-0.5)/len;
    if strcmpi(tail, 'both')
        pvals = min(pvals, 1-pvals);
        pvals = 2*pvals;
    end
end
    
function val = myndims(a)
    if iscolumn(a)
        val = 1;
    else
        val = ndims(a);
    end
end


% TO AVERAGE CHANNELS' TIME SERIES AND GET ERSP, WE CAN JUST DO THIS:
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
%         print(fh, '-dsvg', fullfile(oud...
%         , ['TOVA_ersp_' cond '_' NM '_' meth{l} '_logfreq' strrep(num2str(inc(i)), '.', 'pt') '-30']))
%         close
%     end