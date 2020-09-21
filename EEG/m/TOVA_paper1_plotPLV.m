% INIT VARS
pth = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extanal/PLV';

ROI = {[7 19 36 21 15 23 28] [68 85 100]
        'parieto-occipital' 'frontal'};
roi = sort([ROI{1, :}]);
tx = biosemi1020(roi);

grp = {'Control' 'ADHD'};
cnd = {'' 'rt'
    'Target-locked' 'Response-locked'};

lbwh = get(0,'ScreenSize');
clrmp = {'greytight'
'whatjet(''what'', [0.9 0.9 0.9], ''stops'', [0 0.15 0.25 0.1 0.1 0.25 0.15])'};
cmap = eval(clrmp{2, 1});
s = 40;


%% Load data
CHLOCS = readlocs(which('chanlocs128_pist.elp'));
load(fullfile(pth, 'data', 'cEEG_PLV_7380433272.mat'))
Cwdws = sldngwdws;
load(fullfile(pth, 'data', 'aEEG_PLV_7380493809.mat'))
Awdws = sldngwdws;
Nwdws = numel(sldngwdws);
clear sldngwdws


%% Stat testing - replace ANOVA with...? bootstrap?

% pPLV = zeros(1, size(conIx, 1));
% for t = 1:size(conIx, 1)
%     tst = [Cplv(:, conIx(t, 1), conIx(t, 2)) Aplv(:, conIx(t, 1), conIx(t, 2))];
%     conNm = [tx{conIx(t, 1)} '_' tx{conIx(t, 2)}];
%     [pPLV(t), testPLV.(conNm), statsPLV.(conNm)] = anova1(tst, grp, 'off');
% end
% pPLV = mafdr(pPLV);

%% PLV whole duration plotting
figure('Position', lbwh)
hold on

plvs = Cwdws{10}{1};
plvCI = Cwdws{10}{2};
plvTs = Cwdws{10}{3};
lh = plot(plvTs, plvs(:, 1, 10), 'LineWidth', 6);
lh.Color = [lh.Color 0.5];
lh = plot(plvTs, plvCI(:, 1, 10, 1), 'Color', lh.Color);
lh.Color = [lh.Color 0.5];
lh = plot(plvTs, plvCI(:, 1, 10, 2), 'Color', lh.Color);
lh.Color = [lh.Color 0.5];

plvs = Awdws{10}{1};
plvCI = Awdws{10}{2};
plvTs = Awdws{10}{3};
lh = plot(plvTs, plvs(:, 1, 10), ':', 'LineWidth', 6, 'Color', lh.Color);
lh = plot(plvTs, plvCI(:, 1, 10, 1), ':', 'Color', lh.Color);
lh = plot(plvTs, plvCI(:, 1, 10, 2), ':', 'Color', lh.Color);
legend({'Ctrl' '' '' 'ADHD'})


%% PLV sliding windows plotting
sbf_plot_sw(1, 10, Cwdws, Awdws, Nwdws, 'lgnd', {'' 'Ctrl' '' '' 'ADHD'})


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
    
    sbf_plot_sw(conIx(p, 1), conIx(p, 2), Cwdws, Awdws, Nwdws, 'lbwh', NaN)
    
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
for swidx = 1:Nwdws - 1
    
    Cplv = Cwdws{swidx}{1};
    Aplv = Awdws{swidx}{1};
    plvTs = Cwdws{swidx}{3};

    mnCplv = squeeze(mean(Cplv));
    mnAplv = squeeze(mean(Aplv));
    [conIx(:, 1), conIx(:, 2)] = ind2sub(size(mnCplv), find(mnCplv)); %Indexing
    tstPLV = (mnCplv - mnAplv);
    mxdff = max(abs(tstPLV), [], 'all');
    tstPLV = (tstPLV + mxdff) / (2 * mxdff);
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

    subplot(2, 3, 1)
    % sh.Colormap = cmap(128:255, :);
    im = image(mnCplv .* cmpN + cmpN);
    title(grp{1}); xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    cb.Title.String = 'remove me!';
    delete(cb)

    subplot(2, 3, 2)
    image(mnAplv .* cmpN + cmpN)
    title(grp{2}); xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    cb.Limits = [cmpN cmpN*2];
    cb.Ticks = cmpN:64:cmpN*2;
    cb.TickLabels = round((cb.Ticks - cmpN) ./ cmpN, 2);
    cb.Title.String = 'PLV';


    subplot(2, 3, 3)
    image(tstPLV .* 256)
    % image((tstPLV + mxdff) / (2 * mxdff) .* 256)
    title([grp{1} ' - ' grp{2}])
    xticks(1:10); xticklabels(tx); yticklabels(tx);

    cb = colorbar('EastOutside');
    cb.Ticks = [1 64:64:256];
    cb.TickLabels = round(cb.Ticks ./ 256 * (2 * mxdff) - mxdff, 3);
    cb.Title.String = 'PLV diff';

    chlocs = CHLOCS(roi);
    [chlocs.labels] = deal(tx{:});

    subplot(2, 3, 4)
    topoplot_connect(Cds, chlocs, 'colormap', cmap(128:255, :), 'showlabels', 1)

    subplot(2, 3, 5)
    topoplot_connect(Ads, chlocs, 'colormap', cmap(128:255, :), 'showlabels', 1)

    subplot(2, 3, 6)
    topoplot_connect(Tds, chlocs, 'colormap', cmap, 'showlabels', 1)


    %%
%     tightfig
    svnm = sprintf('TOVA-PLV-resp_%d:%d', round(plvTs(1)), round(plvTs(end)));
    % svnm = 'TOVA-PLV-response-topo';
    print(fh, '-dpng', fullfile(pth, svnm))
    close

end


%% Subfunctions
function sbf_plot_sw(row, col, C, A, N, varargin)
    
    P = inputParser;

    P.addRequired('row', @isnumeric)
    P.addRequired('col', @isnumeric)
    P.addRequired('C', @iscell)
    P.addRequired('A', @iscell)
    P.addRequired('N', @isnumeric)

    P.addParameter('lbwh', get(0,'ScreenSize'), @isnumeric)
    P.addParameter('s', 40, @isscalar)
    P.addParameter('lgnd', {''}, @iscellstr)

    P.parse(row, col, C, A, N, varargin{:})
    P = P.Results;
    
    if ~isnan(P.lbwh)
        figure('Position', P.lbwh)
    end
    
    plot(C{10}{3}, ones(1, numel(C{10}{3})) * 0.5, 'Color', 'w')
    hold on
    for i = N - 1:-1:1
        plvs = C{i}{1};
        plvCI = C{i}{2};
        plvCIlo = smooth(plvCI(:, row, col, 1), P.s, 'loess');
        plvCIhi = smooth(plvCI(:, row, col, 2), P.s, 'loess');
        plvTs = C{i}{3};
        lh = plot(plvTs, plvs(:, row, col), 'LineWidth', 6);
        lh.Color = [lh.Color 0.5];
        lh = plot(plvTs, plvCIlo, 'LineWidth', 2, 'Color', lh.Color);
        lh.Color = [lh.Color 0.5];
        lh = plot(plvTs, plvCIhi, 'LineWidth', 2, 'Color', lh.Color);
        lh.Color = [lh.Color 0.5];

        plvs = A{i}{1};
        plvCI = A{i}{2};
        plvCIlo = smooth(plvCI(:, row, col, 1), P.s, 'loess');
        plvCIhi = smooth(plvCI(:, row, col, 2), P.s, 'loess');
        plvTs = A{i}{3};
        lh = plot(plvTs, plvs(:, row, col), ':', 'LineWidth', 6, 'Color', lh.Color);
        plot(plvTs, plvCIlo, ':', 'LineWidth', 2, 'Color', lh.Color)
        plot(plvTs, plvCIhi, ':', 'LineWidth', 2, 'Color', lh.Color)
    end
    if ~isempty(P.lgnd{:})
        legend(P.lgnd)
    end
end