%% COMPUTE SPECTRA
timedfilename = ['tova_spectra_' datestr(now, 'THHMMSS')];
topospec = false;
csvout = false;
specdat = struct;
part = {'BL' 'PS' 'all'};
slce = [1 512; 513 1024; 1 1024];
freq = [4 8 10 16];
frlm = [4 16];
pwlm = [0  14];
mplm = 'absmax';
which_p = 1:2;

for g = which_g
for c = which_c
for h = which_h
    eeg = eval([grp{g}(1) 'EEG' cnd{c} '_H' h]);
    eeg = pop_mergeset(eeg, 1:numel(eeg));
    
    for r = which_r %These loops only for generating data
    for p = which_p
        nm = [grp{g} '_' cnd{c} '_H' h '_' part{p} '_' ROI{2, r}];
        slice = slce(p, 1):slce(p, 2);
        specdat.(nm) = ...
            spectopo(eeg.data(:, slice, :), numel(slice), 512 ...
            , 'limits', [frlm pwlm]...
            , 'chanlocs', eeg.chanlocs...
            , 'freqrange', frlm...
            , 'plotchans', ROI{1, r}...
            , 'winsize', 512 ...
            , 'overlap', 384 ...
            , 'plot', 'off'...
            , 'title', [nm ',ws=512,ol=384']);
    end
    end
    if topospec
        nm = [grp{g} '_' cnd{c} '_H' h '_']; %#ok<*UNRCH>
        fh = figure('Position', round(lbwh ./ 1.5), 'Visible', 'off');
        spectopo(eeg.data(:, slice, :), numel(slice), 512 ...
            , 'freq', freq...
            , 'limits', [frlm pwlm]...
            , 'chanlocs', eeg.chanlocs...
            , 'emarker2', {[ROI{1, 1} ROI{1, 2}], 'o', 'k', 3, 1}...%mark ROIs
            , 'freqrange', frlm...
            , 'winsize', 512 ...
            , 'overlap', 384 ...
            , 'maplimits', mplm...
            , 'title', [nm 'topos,ws=512,ol=384']);
    %         , 'plotchans', [ROI{1, 1} ROI{1, 2}]...%plot only ROIs
        colormap(cmap)
        if isfolder(oud)
            print(fh, '-dsvg', fullfile(oud, 'spectra', [nm 'topos']))
        end
        close all;
    end
end
end
end
if ~isempty(fieldnames(specdat)) && isfolder(oud)
    save(fullfile(oud, 'spectra', timedfilename), 'specdat')
end


%% LOAD THE SPECTRAL DATA
if isfile(fullfile(oud, 'spectra', timedfilename))
    load(fullfile(oud, 'spectra', timedfilename))
% elseif isfile(fullfile(oud, 'spectra', 'tova_spectra_v1.mat'))
%     load(fullfile(oud, 'spectra', 'tova_spectra_v1.mat'))
end

%% MAKE A PANEL FIGURE
figh = figure('Position', round(lbwh .* [1 1 0.6 0.5]), 'Color', 'w');
rws = 2;
cls = 4;
plt = 1;
clrs = {[0.5 0 0 0.5] [0 0 0.5 0.5]};
ltyp = {'-' '-.'};
xt = 1:2:frlm(2)-3;
xtl = '';
ytl = pwlm(1):2:pwlm(2);

for r = which_r
for c = which_c
for h = which_h
    subplot(rws, cls, plt)
    if any(plt == [2:4 6:8])
        ytl = '';
    elseif plt == 5
        ytl = pwlm(1):2:pwlm(2);
        xtl = frlm(1):2:frlm(2);
    end
    
    for g = which_g
    for p = which_p
        nm = [grp{g} '_' cnd{c} '_H' h '_' part{p} '_' ROI{2, r}];
        specmean = mean(specdat.(nm));

        lyn = plot(specmean(frlm(1):frlm(2)));
        set(lyn, 'Color', clrs{g} .* p, 'LineStyle', ltyp{p}, 'LineWidth', 1.25)
        set(gca, 'xtick', xt, 'xticklabel', xtl, 'ytick', ytl)
        hold on

    end
    end
    axis([frlm - 3 pwlm])
    if any(plt == 1:4)
        title([Cond{c} ', H' h])
    end
    if plt == 1
        ylabel(sprintf('%s\n10*log_{10}(µV^2/Hz)', upper(ROI{2, r})))
    elseif plt == 5
        ylabel(sprintf('%s\n ', upper(ROI{2, r})))
    elseif plt == 8
        xlabel('Frequency (Hz)')
        lg = legend({'Control BL' 'Control PS' 'ADHD BL' 'ADHD PS'}'...
            , 'box', 'off');
        lg.ItemTokenSize = [15 3];
        lg.Position = lg.Position + 0.005;
    end
    plt = plt + 1;
end
end
end

%% PRINT FIGURE
if isfolder(oud)
    print(figh, '-dsvg', fullfile(oud, 'spectra', timedfilename))
end


%% COMPUTE RT-SPLIT SPECTRA
% for rt = {'lo' 'hi'}
%     for g = which_g
%         eeg = eval([grp{g}(1) 'EEG' rt{:} 'RT']);
% 
%         for p = which_p
%             nm = [grp{g} '_' rt{:} 'RT_' part{p} '_'];
%             slice = slce(p, 1):slce(p, 2);
%             for r = 1:2
%                 fh = figure('Position', round(lbwh ./ 1.5), 'Visible', 'off');
%                 specdat.([nm ROI{2, r}]) = ...
%                     spectopo(eeg.data(:, slice, :), numel(slice), 512 ...
%                     , 'limits', [frlm pwlm]...
%                     , 'chanlocs', eeg.chanlocs...
%                     , 'freqrange', frlm...
%                     , 'plotchans', ROI{1, r}...
%                     , 'winsize', 512 ...
%                     , 'overlap', 384 ...
%                     , 'title', [nm ROI{2, r} ',ws=512,ol=384']);
%                 print(fh, '-dsvg'...
%                         , fullfile(oud, 'spectra', ROI{2, r}, [nm ROI{2, r}]))
%             end
% 
%             fh = figure('Position', round(lbwh ./ 1.5), 'Visible', 'off');
%             spectopo(eeg.data(:, slice, :), numel(slice), 512 ...
%                 , 'freq', freq...
%                 , 'limits', [frlm pwlm]...
%                 , 'chanlocs', eeg.chanlocs...
%                 , 'freqrange', frlm...
%                 , 'plotchans', [ROI{1} ROI{3}]...
%                 , 'winsize', 512 ...
%                 , 'overlap', 384 ...
%                 , 'maplimits', mplm...
%                 , 'title', [nm 'topos,ws=512,ol=384'])
%             colormap(cmap)
%             print(fh, '-dsvg', fullfile(oud, 'spectra', [nm 'topos']))
%             close all;
%         end
%     end
% end


%% OUTPUT SPECTRA VALS
if isfolder(oud) && csvout
    roich = repmat(num2cell(1:8), 1, 4);
    names = cellfun(@(x) repmat({x}, 1, 8), {'4Hz' '8Hz' '10Hz' '16Hz'}, 'Un', 0);
    names = cellfun(@(x, y) [x num2str(y)], [names{:}], roich, 'Un', 0);
    specsub = structfun(@(x) x(:, [4 8 10 16]), specdat, 'Un', 0);
    specflat = structfun(@(x) x(:), specsub, 'Un', 0);
    specTab = array2table(struct2mat(specflat)...
                  , 'VariableNames', fieldnames(specdat)...
                  , 'RowNames', names);
    writetable(specTab, fullfile(oud, 'spectra', [timedfilename '.csv'])...
                  , 'WriteRowNames', true)
end

clear which_* spec* roich names c g r p h fh nm ROI...
    freq frlm mplm pwlm slce slice part
