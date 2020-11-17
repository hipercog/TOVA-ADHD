% INIT VARS
pth = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extanal/PLV';

roi = {[7 19 36 21 15 23 28] [68 85 100]
        'parieto-occipital' 'frontal'};
roi = sort([roi{1, :}]);
tx = biosemi1020(roi);

grp = {'Control' 'ADHD'};

lbwh = get(0,'ScreenSize');
cmap = {'greytight'
'whatjet(''what'', [0.9 0.9 0.9], ''stops'', [0 0.15 0.25 0.1 0.1 0.25 0.15])'};
cmap = eval(cmap{2, 1});


%% Load data
CHLOCS = readlocs(which('chanlocs128_pist.elp'));
fs = dir(fullfile(pth, 'data', '*EEG_PLV*.mat'));
load(fullfile(pth, 'data', fs(1).name))
Awdws = sldngwdws;
load(fullfile(pth, 'data', fs(2).name))
Cwdws = sldngwdws;
Nwdws = sum(~cellfun(@isempty, sldngwdws));

for w = 1:Nwdws
    Cwdws{w}{4} = false(size(Cwdws{w}{1}));
    Awdws{w}{4} = false(size(Awdws{w}{1}));
end

tmp = squeeze(mean(Cwdws{1}{1}));
[conIx(:, 1), conIx(:, 2)] = ind2sub(size(tmp), find(tmp)); %Indexing
clear sldngwdws tmp fs w


%% Stat testing - with bootstrap 95% CIs
for c = 1:size(conIx, 1)
    ttl = sprintf('plv:%s-%s', tx{conIx(c, 1)}, tx{conIx(c, 2)});
    
    Cwdws{Nwdws}{4}(:, conIx(c, 1), conIx(c, 2)) = ...
        Cwdws{Nwdws}{2}(:, conIx(c, 1), conIx(c, 2), 1) >...
        Awdws{Nwdws}{2}(:, conIx(c, 1), conIx(c, 2), 2);
    Awdws{Nwdws}{4}(:, conIx(c, 1), conIx(c, 2)) = ...
        Awdws{Nwdws}{2}(:, conIx(c, 1), conIx(c, 2), 1) >...
        Cwdws{Nwdws}{2}(:, conIx(c, 1), conIx(c, 2), 2);

    for swidx = 1:Nwdws - 1
        Cwdws{swidx}{4}(:, conIx(c, 1), conIx(c, 2)) = ...
            Cwdws{swidx}{2}(:, conIx(c, 1), conIx(c, 2), 1) >...
            Awdws{swidx}{2}(:, conIx(c, 1), conIx(c, 2), 2);
        Awdws{swidx}{4}(:, conIx(c, 1), conIx(c, 2)) = ...
            Awdws{swidx}{2}(:, conIx(c, 1), conIx(c, 2), 1) >...
            Cwdws{swidx}{2}(:, conIx(c, 1), conIx(c, 2), 2);
        
        % write aggregated CI-diff test result for whole trial
        if swidx > 1
            w1ts = Cwdws{swidx - 1}{3};
            w2ts = Cwdws{swidx}{3};
            w1ix = ismember(w1ts, w2ts);
            w2ix = ismember(w2ts, w1ts);
            allx = ismember(Cwdws{Nwdws}{3}, w1ts(w1ix));
            Cwdws{Nwdws}{4}(allx, conIx(c, 1), conIx(c, 2)) = ...
                Cwdws{Nwdws}{4}(allx, conIx(c, 1), conIx(c, 2)) &...
                Cwdws{swidx - 1}{4}(w1ix, conIx(c, 1), conIx(c, 2)) &...
                Cwdws{swidx}{4}(w2ix, conIx(c, 1), conIx(c, 2));
            Awdws{Nwdws}{4}(allx, conIx(c, 1), conIx(c, 2)) = ...
                Awdws{Nwdws}{4}(allx, conIx(c, 1), conIx(c, 2)) &...
                Awdws{swidx - 1}{4}(w1ix, conIx(c, 1), conIx(c, 2)) &...
                Awdws{swidx}{4}(w2ix, conIx(c, 1), conIx(c, 2));
        end
    end
end


%% Visualise bootstrap 95% CIs

swidx = Nwdws;
for c = 1:size(conIx, 1)
    ttl = sprintf('plv:%s-%s', tx{conIx(c, 1)}, tx{conIx(c, 2)});
    CgtA = Cwdws{swidx}{4}(:, conIx(c, 1), conIx(c, 2));
    AgtC = Awdws{swidx}{4}(:, conIx(c, 1), conIx(c, 2));
    
    fh = sbf_plot_sw(conIx(c, 1), conIx(c, 2), Cwdws, Awdws...
        , 'CgtA', CgtA...
        , 'AgtC', AgtC...
        , 'sldg', false...
        , 'smth', 4 ...
        , 'ttl', ttl);
    print(fh, '-dpng', fullfile(pth, sprintf('TOVA-PLV-testCI-%s', ttl)))
    close
end


%% PLV whole duration plotting
figure('Position', lbwh)
hold on

plvs = Cwdws{Nwdws}{1};
plvCI = Cwdws{Nwdws}{2};
plvTs = Cwdws{Nwdws}{3};
lh = plot(plvTs, plvs(:, 1, 10), 'LineWidth', 6);
lh.Color = [lh.Color 0.5];
lh = plot(plvTs, plvCI(:, 1, 10, 1), 'Color', lh.Color);
lh.Color = [lh.Color 0.5];
lh = plot(plvTs, plvCI(:, 1, 10, 2), 'Color', lh.Color);
lh.Color = [lh.Color 0.5];

plvs = Awdws{Nwdws}{1};
plvCI = Awdws{Nwdws}{2};
plvTs = Awdws{Nwdws}{3};
lh = plot(plvTs, plvs(:, 1, 10), ':', 'LineWidth', 6, 'Color', lh.Color);
lh = plot(plvTs, plvCI(:, 1, 10, 1), ':', 'Color', lh.Color);
lh = plot(plvTs, plvCI(:, 1, 10, 2), ':', 'Color', lh.Color);
legend({'Ctrl' '' '' 'ADHD'})


%% PLV sliding windows plotting
% sbf_plot_sw(1, 2, Cwdws, Awdws, 'sldg', false)
sbf_plot_sw(6, 10, Cwdws, Awdws, 'N', 9)


%% PLV time course plot matrix
rows = 9;
cols = 9;
% dsg = find(cEEG.times > tms(1), 1) + 1 : find(cEEG.times <= tms(2), 1, 'last');
% lgnd = '';
pltix = zeros(9, 9);
pltix(1:81) = 1:81;
pltix = pltix';

figh = figure('Position', lbwh, 'Color', 'w');

for p = 1:size(conIx, 1)
    
    subplot(rows, cols, pltix(conIx(p, 1), conIx(p, 2) - 1))
% %     erp = [Cplv(:, conIx(i, 1), conIx(i, 2)) Aplv(:, conIx(i, 1), conIx(i, 2))];
% %     erp = [Cplv(:, conIx(i, 1), conIx(i, 2)) Aplv(:, conIx(i, 1), conIx(i, 2))];
% %     cCIs = [CplvCI(:, conIx(p, 1), conIx(p, 2), 1) CplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
% %     aCIs = [AplvCI(:, conIx(p, 1), conIx(p, 2), 1) AplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
%     cCIs = [CloRTplvCI(:, conIx(p, 1), conIx(p, 2), 1) CloRTplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
%     aCIs = [ChiRTplvCI(:, conIx(p, 1), conIx(p, 2), 1) ChiRTplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
%     sh=plot([cCIs aCIs]);
% %     sh = ctap_plot_basic_erp(erp', find(cEEG.times(dsg) == 0), cEEG.srate...
% %                     , 'timeunit', 'ms'...
% %                     , 'overploterp', CIs'...
% %                     , 'lgnd', lgnd...
% %                     , 'lgndloc', 'southeast');
% % %                     , 'ylimits', [0.15 0.9]...
    
    sbf_plot_sw(conIx(p, 1), conIx(p, 2), Cwdws, Awdws, 'lbwh', NaN, 'sldg', false)
    
%     if p == size(conIx, 1) - 1
% %         lgnd = {'Ctrl PLV' 'ADHD PLV' 'Ctrl CI1' 'Ctrl CI2'};
% %         lgnd = {'Ctrl:rt<md' 'Ctrl:rt>md' 'ADHD:rt<md' 'ADHD:rt>md'};
%     end
%     if p < size(conIx, 1)
%         xlabel(''); xticklabels({});
%     else
%         xlabel('Time (ms)')%; xtickangle(sh, 45)
%         legend({'Ctrl CI1' 'Ctrl CI2' 'ADHD CI1' 'ADHD CI2'})
%     end
%     if p == 1, ylabel('PLV'); else, ylabel(''); end
%     title([tx{conIx(p, 1)} '<>' tx{conIx(p, 2)}])
        % ', p=' num2str(round(pPLV(i), 4), 4)]);
end


%% PLV lead-matrix and scalp-map plotting
mxall = 0.121;
for swidx = 1:Nwdws
    
    Cplv = Cwdws{swidx}{1};
    Aplv = Awdws{swidx}{1};
    plvTs = Cwdws{swidx}{3};

    mnCplv = squeeze(mean(Cplv));
    mnAplv = squeeze(mean(Aplv));
    [conIx(:, 1), conIx(:, 2)] = ind2sub(size(mnCplv), find(mnCplv)); %Indexing
    tstPLV = (mnCplv - mnAplv);
%     mxdff = max(abs(tstPLV), [], 'all');
    tstPLV = (tstPLV + mxall) / (2 * mxall);
    plvlim = [min(cat(3, mnCplv, mnAplv), [], 'all')...
              max(cat(3, mnCplv, mnAplv), [], 'all')];

    %topoplot_connect structures
    Cds.chanPairs = conIx;
    Cds.connectStrength = mnCplv(mnCplv > 0);
    Cds.connectStrengthLimits = [0 1];%plvlim;

    Ads.chanPairs = conIx;
    Ads.connectStrength = mnAplv(mnAplv > 0);
    Ads.connectStrengthLimits = [0 1];%plvlim;

    Tds.chanPairs = conIx;
    Tds.connectStrength = tstPLV(tstPLV ~= 0.5);
    Tds.connectStrengthLimits = [0 1];%[min(pPLV) max(pPLV)];

    % and figure
    fh = figure('Position', lbwh, 'Color', 'w');%[lbwh(1:3) lbwh(4) * 0.5])

    colormap(cmap)
    cmpN = 128;

    % all-to-all matrix plots
    subplot(2, 3, 1)
    im = image(mnCplv .* cmpN + cmpN);
    title(grp{1}); xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    cb.Limits = [cmpN cmpN*2];
    cb.Ticks = cmpN:64:cmpN*2;
    cb.TickLabels = round((cb.Ticks - cmpN) ./ cmpN, 2);
    cb.Title.String = 'PLV';

    subplot(2, 3, 2)
    image(cmpN - mnAplv .* cmpN)
    title(grp{2}); xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    set(cb, 'YDir', 'reverse')
    cb.Limits = [cmpN cmpN*2];
    cb.Ticks = cmpN:64:cmpN*2;
    cb.TickLabels = round((flip(cb.Ticks) - cmpN) ./ cmpN, 2);
    cb.Title.String = 'PLV';

    subplot(2, 3, 3)
    image(tstPLV .* 256)
    title([grp{1} ' - ' grp{2}]); xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    cb.Ticks = [1 64:64:256];
    cb.TickLabels = round(cb.Ticks ./ 256 * (2 * mxall) - mxall, 3);
    cb.Title.String = 'PLV diff';

    % Cartoon head plots
    chlocs = CHLOCS(roi);
    [chlocs.labels] = deal(tx{:});

    subplot(2, 3, 4)
    topoplot_connect(Cds, chlocs, 'colormap', cmap(128:255, :), 'showlabels', 1)

    subplot(2, 3, 5)
    topoplot_connect(Ads, chlocs, 'colormap', flip(cmap(1:128, :)), 'showlabels', 1)

    subplot(2, 3, 6)
    topoplot_connect(Tds, chlocs, 'colormap', cmap, 'showlabels', 1)


    %%
%     tightfig
    if swidx < 10
        nm = sprintf('%s%d', 'TOVA-PLV-cR-wdw', swidx);
    else
        nm = 'TOVA-PLV-corResp';
    end
    svnm = sprintf('%s_%d:%d', nm, round(plvTs(1)), round(plvTs(end)));
    print(fh, '-dpng', fullfile(pth, svnm))
    close

end


%% Subfunctions
function fh = sbf_plot_sw(row, col, C, A, varargin)
    
    P = inputParser;

    P.addRequired('row', @isnumeric)
    P.addRequired('col', @isnumeric)
    P.addRequired('C', @iscell)
    P.addRequired('A', @iscell)

    P.addParameter('N', sum(~cellfun(@isempty, A)), @isnumeric)
    P.addParameter('lbwh', get(0,'ScreenSize'), @isnumeric)
    P.addParameter('smth', 40, @isscalar)
    P.addParameter('lgnd', {'Ctrl', 'ADHD'}, @iscellstr)
    P.addParameter('sldg', true, @islogical)
    P.addParameter('ttl', '', @ischar)
    P.addParameter('CgtA', 0, @islogical)
    P.addParameter('AgtC', 0, @islogical)

    P.parse(row, col, C, A, varargin{:})
    P = P.Results;
    
    if ~isnan(P.lbwh)
        fh = figure('Position', P.lbwh);
    else
        fh = gcf;
    end
    
    xlim([C{10}{3}(1) C{10}{3}(end)])
    hold on
    if P.sldg
        wdws = P.N:-1:1;
    else
        wdws = P.N;
    end
    linemult = ceil(numel(wdws) / 5);
    for i = wdws
        plvTs = C{i}{3};
        
        lh = plot(plvTs, C{i}{1}(:, row, col)...
            , 'LineWidth', 3 * linemult...
            , 'HandleVisibility', 'off');
        lh.Color = [lh.Color 0.5];
        clrl1 = lh.Color;
        
        % CTRL CI filled-patch
        plvCIlo = smooth(C{i}{2}(:, row, col, 1), P.smth, 'loess')';
        plvCIhi = smooth(C{i}{2}(:, row, col, 2), P.smth, 'loess')';
        patch([plvTs fliplr(plvTs)], [plvCIlo fliplr(plvCIhi)], [0.9 0.9 0.9]...
            , 'FaceAlpha', 0.25 ...
            , 'EdgeColor', lh.Color...
            , 'HandleVisibility', 'off')
%         lh = plot(plvTs, plvCIlo, 'LineWidth', 1 * linemult, 'Color', lh.Color);
%         lh.Color = [lh.Color 0.5];
%         lh = plot(plvTs, plvCIhi, 'LineWidth', 1 * linemult, 'Color', lh.Color);
%         lh.Color = [lh.Color 0.5];

        % ADHD CI filled-patch
        plvCIlo = smooth(A{i}{2}(:, row, col, 1), P.smth, 'loess')';
        plvCIhi = smooth(A{i}{2}(:, row, col, 2), P.smth, 'loess')';
        patch([plvTs fliplr(plvTs)], [plvCIlo fliplr(plvCIhi)], [0.9 0.9 0.9]...
            , 'FaceAlpha', 0.25 ...
            , 'LineStyle', ':'...
            , 'EdgeColor', lh.Color...
            , 'HandleVisibility', 'off')
        
        plot(plvTs, C{i}{1}(:, row, col)...
            , 'LineWidth', 3 * linemult, 'Color', clrl1)
        plot(plvTs, A{i}{1}(:, row, col), ':'...
            , 'LineWidth', 3 * linemult, 'Color', lh.Color)
%         plot(plvTs, plvCIlo, ':', 'LineWidth', 1 * linemult, 'Color', lh.Color)
%         plot(plvTs, plvCIhi, ':', 'LineWidth', 1 * linemult, 'Color', lh.Color)
    end
    ylm = ylim;
    nany = plvTs;
    if any(P.CgtA)
        nany(~P.CgtA') = NaN;
        plot(nany, ones(1, numel(plvTs)) .* ylm(1), 'Color', 'r', 'LineWidth', 3)
        CgtApc = round(sum(P.CgtA) / numel(P.CgtA) * 100);
        res = sprintf('C>A %d%%', CgtApc);
        P.lgnd{end + 1} = res;
        P.ttl = [P.ttl ', ' res];
    end
    if any(P.AgtC)
        nany(~P.AgtC') = NaN;
        plot(nany, ones(1, numel(plvTs)) .* ylm(1), 'Color', 'b', 'LineWidth', 3)
        AgtCpc = round(sum(P.AgtC) / numel(P.AgtC) * 100);
        res = sprintf('A>C %d%%', AgtCpc);
        P.lgnd{end + 1} = res;
        P.ttl = [P.ttl ', ' res];
    end
    if ~all(cellfun(@isempty, P.lgnd))
        legend(P.lgnd)
    end
    if ~isempty(P.ttl)
        title(P.ttl)
    end
end