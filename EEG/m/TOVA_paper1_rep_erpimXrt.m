ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
% ind = '/wrk/group/hipercog/project_TOVA/ANALYSIS/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extended_anal';
% oud = '/wrk/group/hipercog/project_TOVA/ANALYSIS/paper1/PLV';

ROI = {[7 19 36 21 15 23 28] [68 85 100]
        'parieto-occipital' 'frontal'};
grp = {'Control' 'ADHD'};
cnd = {'' 'rt'
    'Target-locked' 'Response-locked'};

lbwh = get(0,'ScreenSize');
clrmp = {'greytight'
'whatjet(''what'', [0.9 0.9 0.9], ''stops'', [0 0.15 0.25 0.1 0.1 0.25 0.15])'};
cmap = eval(clrmp{2, 1});


%% CORRECT RESPONSE DATA
cEEG = pop_loadset('filepath', ind...
                 , 'filename', 'CONTROL_MERGED_randomised_COR_RSP_TOTAL.set');
aEEG = pop_loadset('filepath', ind...
                 , 'filename', 'ADHD_MERGED_randomised_COR_RSP_TOTAL.set');
%% epoch aligned to RT
ix = [cEEG.event(ismember({cEEG.event.type}, 'cor_rsp')).duration] > 600;
cEEGrt = pop_epoch(pop_select(cEEG, 'notrial', ix), {'RESP'}, [-0.6 0.4]);
ix = [aEEG.event(ismember({aEEG.event.type}, 'cor_rsp')).duration] > 600;
aEEGrt = pop_epoch(pop_select(aEEG, 'notrial', ix), {'RESP'}, [-0.6 0.4]);
%% get structs for group-wise erp-images - NOT NEEDED?
% control_erpimXrt = std_erpimage(ctrlEEG, 'channels', 1:128 ...
%     , 'nlines', 1500, 'smoothing', 100, 'recompute', 'on'...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration'...
%     , 'fileout', fullfile(oud, cttl));
% adhdin_erpimXrt = std_erpimage(adhdEEG, 'channels', 1:128 ...
%     , 'nlines', 4100, 'smoothing', 285, 'recompute', 'on'...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration'...
%     , 'fileout', fullfile(oud, attl));


%% TESTING: test negative and positive peaks around RT
tstwn = [150 250 330 430; -170 -70 0 100];% in ms
meanps = ones(size(tstwn, 2) / 2, size(ROI, 2), size(cnd, 2));
latps = ones(size(tstwn, 2) / 2, size(ROI, 2), size(cnd, 2));

for c = 1:size(cnd, 2)
    
for r = 1:size(ROI, 2)

    CEEG = eval(['cEEG' cnd{1, c}]);
    AEEG = eval(['aEEG' cnd{1, c}]);

    diffNA = NaN(1, AEEG.trials - CEEG.trials);
    wdwlab = {'2' '3'
            'pre' 'post'};

    trodes = ROI{1, r};
    tstwdw = arrayfun(@(x) find(CEEG.times >= x, 1), tstwn(c, :));
    
    if numel(trodes) > 1
        cN = squeeze(mean(CEEG.data(trodes, tstwdw(1):tstwdw(2), :)));
        aN = squeeze(mean(AEEG.data(trodes, tstwdw(1):tstwdw(2), :)));
        cP = squeeze(mean(CEEG.data(trodes, tstwdw(3):tstwdw(4), :)));
        aP = squeeze(mean(AEEG.data(trodes, tstwdw(3):tstwdw(4), :)));
    else
        cN = squeeze(CEEG.data(trodes, tstwdw(1):tstwdw(2), :));
        aN = squeeze(AEEG.data(trodes, tstwdw(1):tstwdw(2), :));
        cP = squeeze(CEEG.data(trodes, tstwdw(3):tstwdw(4), :));
        aP = squeeze(AEEG.data(trodes, tstwdw(3):tstwdw(4), :));
    end

    [pNaov, testTN, statsN] = anova1([[mean(cN) diffNA]; mean(aN)].', grp, 'off'); %#ok<*ASGLU>
    [hN, pNks, ks2statN] = kstest2(mean(cN), mean(aN));
    [pPaov, testTP, statsP] = anova1([[mean(cP) diffNA]; mean(aP)].', grp, 'off');
    [hP, pPks, ks2statP] = kstest2(mean(cP), mean(aP));
    testTN = [testTN {'KS-test'; pNks; hN; ks2statN}]; %#ok<*AGROW>
    testTP = [testTP {'KS-test'; pPks; hP; ks2statP}];

    svnm = ['anova-wdw-mean-N' wdwlab{c, 1} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTN(2:end, :), 'VariableNames', testTN(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))
    svnm = ['anova-wdw-mean-P' wdwlab{c, 2} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTP(2:end, :), 'VariableNames', testTP(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))

    [~, lCN] = min(cN); lCN = CEEG.times(lCN + tstwdw(1) - 1);
    [~, lAN] = min(aN); lAN = AEEG.times(lAN + tstwdw(1) - 1);
    [~, lCP] = min(cP); lCP = CEEG.times(lCP + tstwdw(3) - 1);
    [~, lAP] = min(aP); lAP = AEEG.times(lAP + tstwdw(3) - 1);

    [plNaov, testTlN, statslN] = anova1([[lCN diffNA]; lAN].', grp, 'off');
    [hlN, plNks, ks2statlN] = kstest2(lCN, lAN);
    [plPaov, testTlP, statslP] = anova1([[lCP diffNA]; lAP].', grp, 'off');
    [hlP, plPks, ks2statlP] = kstest2(lCP, lAP);
    testTlN = [testTlN {'KS-test'; plNks; hlN; ks2statlN}];
    testTlP = [testTlP {'KS-test'; plPks; hlP; ks2statlP}];

    svnm = ['anova-wdw-lat-N' wdwlab{c, 1} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTlN(2:end, :), 'VariableNames', testTlN(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))
    svnm = ['anova-wdw-lat-P' wdwlab{c, 2} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTlP(2:end, :), 'VariableNames', testTlP(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))

    meanps(:, r, c) = [pNks pPks];
    latps(:, r, c) = [plNks plPks];
end %roi
end %condition
clear CEEG AEEG

% Storey FDR multi-comparison correction
pfdr = [reshape(meanps, [1 numel(meanps)]) reshape(latps, [1 numel(latps)])];
pfdr = mafdr(pfdr, 'lambda', 10^-8:0.001:0.05);
meanps = reshape(pfdr(1:8), [2 2 2]);


%% PLOTTING: 
% plot sorted by RT
% fh = sbf_erpim2x2(cEEG, aEEG, 100, 285, -100 ...
%                 , ROI{1, 1}, ROI{1, 2}, ctl, atl, ROI{2, 1}, ROI{2, 2}...
%                 , ttl, [-8.1 8.1], {'cor_rsp'}, 'duration');
% print(fh, '-dsvg', fullfile(oud, 'erpimXrt', 'TOVA-ERPIMxRT'))
% % plot sorted by TARGET
% fh = sbf_erpim2x2(cEEGrt, aEEGrt, 100, 285, -100 ...
%                 , ROI{1, 1}, ROI{1, 2}, ctl, atl, ROI{2, 1}, ROI{2, 2}...
%                 , ttl, [-8.1 8.1], {'TRGT'}, 'latency');
% print(fh, '-dsvg', fullfile(oud, 'erpimXrt', 'TOVA-RTlocked-ERPIMxTarget'))


%% PLOT ERPIM x RT
rows = 3;
cols = 4;
cax = [-8 8];
cbr = 'off';
startms = [-200 -600];
sortype = {'cor_rsp' 'TRGT'};
sortfld = {'duration' 'latency'};
smth = [100 285];
ylab = 'Reaction Time (ms)';
erp = cell(numel(grp), size(cnd, 2), size(ROI, 2));

plt = 1;

fh = figure('Units','pixels',...
        'Position', lbwh .* [1 1 0.75 1],...
        'Color', 'w');
% subplotting ERPIMs
for g = 1:numel(grp)
    for c = 1:size(cnd, 2)
        
        eeg = eval([lower(grp{g}(1)) 'EEG' cnd{1, c}]);
        dsg = find(eeg.times > startms(c), 1);
        dsg = dsg:min(eeg.pnts, dsg + eeg.srate - 1);
        
        for r = 1:size(ROI, 2)
            
            subplot(rows, cols, plt)
            if plt < 5
                ttl = sprintf('%s\n%s', cnd{2, c}, ROI{2, r});
            else
                ttl = '';
            end
            
            [~, ~, ~, ~, axh, erp{g, c, r}] =...
                erpimage( mean(eeg.data(ROI{1, r}, dsg, :), 1)...
                    , eeg_getepochevent( eeg, {sortype{c}}, [], sortfld{c})...
                    , [startms(c) numel(dsg) eeg.srate]...
                    , ttl, smth(g), 1 ...
                    , 'erp', 1 ...
                    , 'erpalpha', 0.001 ...
                    , 'yerplabel', ''...
                    , 'caxis', cax...
                    , 'cbar', cbr...
                    , 'cbar_title', '\muV'...
                    , 'img_trialax_label', ylab...
                    , 'noxlabel', 'on'); %#ok<CCAT1>

            % ERPIM plot tweaking
            if any(plt == [1 5])
                ylabel(axh{1}, sprintf('%s\n%s', upper(grp{g}), ylab)...
                    , 'FontWeight', 'bold'...
                    , 'Units', 'normalized'...
                    , 'Position', [-0.2 0.5]...
                    , 'Rotation', 90)
            else
                ylabel(axh{1}, '')
            end
            if plt == 3
                cbr = 'on';
            else 
                cbr = 'off';
            end
            % colorbar tweaking
            if plt == 4
                h = axh{2}.Position(4) / 5;
                axh{2}.Position(3:4) = axh{2}.Position(3:4) .* 0.8;
                axh{2}.Position(2) = axh{2}.Position(2) + h;
            end
            % ERP underplot tweaking
            xtickangle(axh{3}, 45)
            set(axh{3}, 'Color', 'w'...
                      , 'YLim', [-3 4]...
                      , 'YAxisLocation', 'right'...
                      , 'YTick', [-3 0 4]...
                      , 'YTickLabel', {'-3' '0' '4'})
            if any(plt == [4 8])
                ylabel(axh{3}, '\muV'...
                    , 'FontWeight', 'bold'...
                    , 'Units', 'normalized'...
                    , 'Position', [1.15 1]...
                    , 'Rotation', 0)
                xlabel(axh{3}, sprintf('Time\n(ms)')...
                    , 'FontWeight', 'bold'...
                    , 'Units', 'normalized'...
                    , 'Position', [0.9 -0.5]...
                    , 'Rotation', 45)
            end
            plt = plt + 1;
        end %ROI
    end %condition
end %group
clear eeg


%% PLOT subplots for ERPs
lgnd = '';
xlab = '';
ylab = sprintf('%s vs %s\n%s', grp{1}, grp{2}, 'amplitude (\muV)');

for c = 1:size(cnd, 2)

    for r = 1:size(ROI, 2)

        CEEG = eval(['cEEG' cnd{1, c}]);
        AEEG = eval(['aEEG' cnd{1, c}]);
        dsg = find(CEEG.times > startms(c), 1);
        dsg = dsg:min(CEEG.pnts, dsg + CEEG.srate - 1);
        
        % ERP of ADHD v CTRL aligned to RT
        % First single-chan ERP for central 'trode (Pz, Fz): at index 2 in ROIs
        CEz = ctap_get_erp(CEEG, 'pnts', dsg, 'loc', ROI{1, r}(2));
        AEz = ctap_get_erp(AEEG, 'pnts', dsg, 'loc', ROI{1, r}(2));
        % Then ERP plus 95%CI of ROIs
        AvCerp = zeros(2, numel(dsg));
        AvCitv = zeros(2, numel(dsg), 2);
        [AvCerp(1, :), AvCitv(1, :, :)] = ctap_get_erp(CEEG...
        , 'pnts', dsg, 'loc', ROI{1, r}, 'dispersion', 'bootci', 'nboot', 10); % TODO : REMEMBER TO REMOVE NBOOT=10!!
        [AvCerp(2, :), AvCitv(2, :, :)] = ctap_get_erp(AEEG...
        , 'pnts', dsg, 'loc', ROI{1, r}, 'dispersion', 'bootci', 'nboot', 10);

        % Plotting
        subplot(rows, cols, plt)
        
        if any(plt == [11 12])
            Ez = biosemi1020(ROI{1, r}(2));
            lgnd = {'Control:ROI' 'ADHD:ROI' ['Control:' Ez] ['ADHD:' Ez]};
        end
        
        sh = ctap_plot_basic_erp(...
            AvCerp, find(CEEG.times(dsg) == 0), CEEG.srate...
            , 'areas', AvCitv...
            , 'overploterp', [CEz; AEz]...
            , 'lgnd', lgnd...
            , 'vlines', tstwn(c, :)...
            , 'testps', meanps(:, r, c)...
            , 'xlbl', xlab...
            , 'ylbl', ylab...
            , 'ylimits', [-3 4]...
            , 'timeunit', 'ms'...
            , 'tkjump', 200);
        
        % tweak plot elements
        xtickangle(sh, 45)
        if plt == 9
            ylab = '';
        end
        if plt == 12
            xlabel(sh, sprintf('Time (ms)')...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [1 0]...
                , 'Rotation', 45)
        end
        plt = plt + 1;
        
    end %roi
end %condition
clear CEEG AEEG


%% PRINT
colormap(cmap)
print(fh, '-dsvg', fullfile(oud, 'erpimXrt', 'TOVA-ERPIMxRT'))
print(fh, '-dpng', fullfile(oud, 'erpimXrt', 'TOVA-ERPIMxRT'))
close


%% erp, src, pnts, zeropt, srate, tkoffset, lgnd, ttl, savename
% fh = ctap_plot_basic_erp(...
%     AvCerp, find(cEEGrt.times == 0), cEEGrt.srate...
%     , 'waves', AvCsrc...
%     , 'lgnd', {'Control' 'ADHD'}...
%     , 'ttl', ttl...
%     , 'vlines', tstwn...
%     , 'cmap', cmap);
% 
% print(fh, '-dsvg', fullfile(oud, 'erpimXrt', ttl))


Crts = eeg_getepochevent(cEEG, {'cor_rsp'}, [], 'duration');
idx = Crts < median(Crts);
cEEGloRT = pop_select(cEEG, 'trial', find(idx));
cEEGhiRT = pop_select(cEEG, 'trial', find(~idx));
Arts = eeg_getepochevent(aEEG, {'cor_rsp'}, [], 'duration');
idx = Arts < median(Arts);
aEEGloRT = pop_select(aEEG, 'trial', find(idx));
aEEGhiRT = pop_select(aEEG, 'trial', find(~idx));


%% Phase-locking tests
% INIT VARS
roi = sort([ROI{1, :}]);
tx = biosemi1020(roi);
filtSpec.range = [6 10];
filtSpec.order = 300;
nbootci = 100;

wdwinc = 100;
wdwstarts = -200:wdwinc:600;
sldngwdws = cell(3, 6, numel(wdwstarts));

for sw = 1:numel(sldngwdws)
    tms = [wdwstarts(sw) wdwstarts(sw) + wdwinc * 2];
    % calculate PLV & bootstrap 95% CIs of whole epoch...
    [Cplv, CplvCI, Ctime] = sbf_get_plv(cEEG, roi, tms, filtSpec, nbootci);
    [Aplv, AplvCI, Atime] = sbf_get_plv(aEEG, roi, tms, filtSpec, nbootci);

    % ...and for +/- median RT
    [CloRTplv, CloRTplvCI, CloTime] = ...
                        sbf_get_plv(cEEGloRT, roi, tms, filtSpec, nbootci);
    [ChiRTplv, ChiRTplvCI, ChiTime] = ...
                        sbf_get_plv(cEEGhiRT, roi, tms, filtSpec, nbootci);
    [AloRTplv, AloRTplvCI, AloTime] = ...
                        sbf_get_plv(aEEGloRT, roi, tms, filtSpec, nbootci);
    [AhiRTplv, AhiRTplvCI, AhiTime] = ...
                        sbf_get_plv(aEEGhiRT, roi, tms, filtSpec, nbootci);
    sldngwdws{sw} = {Cplv, CplvCI, Ctime
                     Aplv, AplvCI, Atime
                     CloRTplv, CloRTplvCI, CloTime
                     ChiRTplv, ChiRTplvCI, ChiTime
                	 AloRTplv, AloRTplvCI, AloTime
                     AhiRTplv, AhiRTplvCI, AhiTime};
end

tms = [-200 800];

% calculate PLV & bootstrap 95% CIs of whole epoch...
[Cplv, CplvCI, Ctime] = sbf_get_plv(cEEG, roi, tms, filtSpec, nbootci);
[Aplv, AplvCI, Atime] = sbf_get_plv(aEEG, roi, tms, filtSpec, nbootci);

% ...and for +/- median RT
[CloRTplv, CloRTplvCI, CloTime] = ...
                    sbf_get_plv(cEEGloRT, roi, tms, filtSpec, nbootci);
[ChiRTplv, ChiRTplvCI, ChiTime] = ...
                    sbf_get_plv(cEEGhiRT, roi, tms, filtSpec, nbootci);
[AloRTplv, AloRTplvCI, AloTime] = ...
                    sbf_get_plv(aEEGloRT, roi, tms, filtSpec, nbootci);
[AhiRTplv, AhiRTplvCI, AhiTime] = ...
                    sbf_get_plv(aEEGhiRT, roi, tms, filtSpec, nbootci);

sldngwdws{end + 1} = {Cplv, CplvCI, Ctime
                     Aplv, AplvCI, Atime
                     CloRTplv, CloRTplvCI, CloTime
                     ChiRTplv, ChiRTplvCI, ChiTime
                	 AloRTplv, AloRTplvCI, AloTime
                     AhiRTplv, AhiRTplvCI, AhiTime};


%% Stat testing - replace ANOVA with...? bootstrap?
% 
% mnCplv = squeeze(mean(Cplv));
% mnAplv = squeeze(mean(Aplv));
% 
% % Indexing
% [conIx(:, 1), conIx(:, 2)] = ind2sub(size(mnCplv), find(mnCplv));
% 
% pPLV = zeros(1, size(conIx, 1));
% for t = 1:size(conIx, 1)
%     tst = [Cplv(:, conIx(t, 1), conIx(t, 2)) Aplv(:, conIx(t, 1), conIx(t, 2))];
%     conNm = [tx{conIx(t, 1)} '_' tx{conIx(t, 2)}];
%     [pPLV(t), testPLV.(conNm), statsPLV.(conNm)] = anova1(tst, grp, 'off');
% end
% pPLV = mafdr(pPLV);


%% PLV time course plots
rows = 9;
cols = 9;
dsg = find(cEEG.times > tms(1), 1) + 1 : find(cEEG.times <= tms(2), 1, 'last');
lgnd = '';
pltix = zeros(9, 9);
pltix(1:81) = 1:81;
pltix = pltix';

figh = figure('Position', lbwh, 'Color', 'w');

for p = 1:size(conIx, 1)
    
    subplot(rows, cols, pltix(conIx(p, 1), conIx(p, 2) - 1))
%     erp = [Cplv(:, conIx(i, 1), conIx(i, 2)) Aplv(:, conIx(i, 1), conIx(i, 2))];
%     erp = [Cplv(:, conIx(i, 1), conIx(i, 2)) Aplv(:, conIx(i, 1), conIx(i, 2))];
%     cCIs = [CplvCI(:, conIx(p, 1), conIx(p, 2), 1) CplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
%     aCIs = [AplvCI(:, conIx(p, 1), conIx(p, 2), 1) AplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
    cCIs = [CloRTplvCI(:, conIx(p, 1), conIx(p, 2), 1) CloRTplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
    aCIs = [ChiRTplvCI(:, conIx(p, 1), conIx(p, 2), 1) ChiRTplvCI(:, conIx(p, 1), conIx(p, 2), 2)];
    sh=plot([cCIs aCIs]);
%     sh = ctap_plot_basic_erp(erp', find(cEEG.times(dsg) == 0), cEEG.srate...
%                     , 'timeunit', 'ms'...
%                     , 'overploterp', CIs'...
%                     , 'lgnd', lgnd...
%                     , 'lgndloc', 'southeast');
% %                     , 'ylimits', [0.15 0.9]...
    if p == size(conIx, 1) - 1
%         lgnd = {'Ctrl PLV' 'ADHD PLV' 'Ctrl CI1' 'Ctrl CI2'};
%         lgnd = {'Ctrl:rt<md' 'Ctrl:rt>md' 'ADHD:rt<md' 'ADHD:rt>md'};
    end
    if p < size(conIx, 1)
        xlabel(''); xticklabels({});
    else
        xlabel('Time (ms)')%; xtickangle(sh, 45)
        legend({'Ctrl CI1' 'Ctrl CI2' 'ADHD CI1' 'ADHD CI2'})
    end
    if p == 1, ylabel('PLV'); else, ylabel(''); end
    title([tx{conIx(p, 1)} '<>' tx{conIx(p, 2)}])
        % ', p=' num2str(round(pPLV(i), 4), 4)]);
end


%% PLV plotting
swidx = 2;
Cplv = sldngwdws{swidx}{7};
Aplv = sldngwdws{swidx}{11};
mnCplv = squeeze(mean(Cplv));
mnAplv = squeeze(mean(Aplv));
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
figure('Position', lbwh, 'Color', 'w')%[lbwh(1:3) lbwh(4) * 0.5])

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
title([grp{1} ' - ' grp{2}]); xticks(1:10); xticklabels(tx); yticklabels(tx);

cb = colorbar('EastOutside');
cb.Ticks = [1 64:64:256];
cb.TickLabels = round(cb.Ticks ./ 256 * (2 * mxdff) - mxdff, 3);
cb.Title.String = 'PLV diff';


chlocs = cEEG.chanlocs(roi);
[chlocs.labels] = deal(tx{:});

subplot(2, 3, 4)
topoplot_connect(Cds, chlocs, 'colormap', cmap(128:255, :), 'showlabels', 1)

subplot(2, 3, 5)
topoplot_connect(Ads, chlocs, 'colormap', cmap(128:255, :), 'showlabels', 1)

subplot(2, 3, 6)
topoplot_connect(Tds, chlocs, 'colormap', cmap, 'showlabels', 1)


%%
tightfig
svnm = 'TOVA-PLV-response-nu';
% svnm = 'TOVA-PLV-response-topo';
print(gcf, '-dpng', fullfile(oud, 'PLV', svnm))


%% Topoplot the neg+pos peaks
inc = 22;
fh = figure('Position', [1 1 lbwh(4) lbwh(3)], 'Color', 'w');
plt = 1;

ctl = 'Control';
atl = 'ADHD';

cax = [-3 1];
for w = 1:2:3
    for tx = tstwn(w):inc:tstwn(w + 1)

        % Plot the Control above
        subplot(4, 3, plt)
        data = mean(mean(cEEGrt.data(:, tx:tx + inc - 2, :), 3), 2);
        topoplot(data, cEEGrt.chanlocs...
            , 'maplimits', cax, 'electrodes', 'on'...
            , 'emarker', {'.','k',[],1}, 'emarker2', {ROI, 'o','m',3,1})
        title(sprintf('%sms\n%s'...
            , num2str(round(cEEGrt.times([tx tx + inc - 1]))), ctl))
        cb = colorbar('EastOutside'); title(cb, '\muV', 'FontSize', 9)
        
        % Plot the ADHD below
        subplot(4, 3, plt + 3)
        data = mean(mean(aEEGrt.data(:, tx:tx + inc - 2, :), 3), 2);
        topoplot(data, aEEGrt.chanlocs...
            , 'maplimits', cax...
            , 'electrodes', 'on'...
            , 'emarker', {'.', 'k', [], 1}...
            , 'emarker2', {ROI, 'o', 'm', 3, 1})
        title(atl)
        colormap(cmap)
        cb = colorbar('EastOutside'); title(cb, '\muV', 'FontSize', 9)
        
        plt = plt + 1;
    end
    plt = plt + 3;
    cax = [-1 3];
end
inc = num2str(round((1000 / cEEGrt.srate) * inc - 2));
% sgt = sgtitle(['Scalp maps at ' inc 'ms intervals around RT'], 'Color', 'r');
text(-3, -1, ['Scalp maps at ' inc 'ms intervals around RT'], ...
  'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
export_fig(fullfile(oud, 'TOVA-topoplots-ERPxRT.pdf'))
print(fh, '-dsvg', fullfile(oud, 'TOVA-topoplots-ERPxRT'))



%% 
% [STUDY ALLEEG] = std_editset( [], []...
%     , 'name', 'TOVA_adhdin_erpimXrt'...
%     , 'task', 'correct_response_locked'...
%     , 'updatedat','on', 'savedat','on', 'rmclust','on'...
%     , 'commands', {{'index', 1, 'load'...
%           , fullfile(ind, 'ADHD_MERGED_randomised_COR_RSP_TOTAL.set')}}...
%     );
% % adhdin_erpimXrt = 
% precomp_wrap(STUDY, 'alleeg', ALLEEG, 'erpim', 1 ...
%     , 'nlines', 4100, 'smoothing', 285 ...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration');
% 
% 
% [STUDY ALLEEG] = std_editset( [], []...
%     , 'name', 'TOVA_control_erpimXrt'...
%     , 'task', 'correct_response_locked'...
%     , 'updatedat','on', 'savedat','on', 'rmclust','on'...
%     , 'commands', {{'index', 1, 'load'...
%           , fullfile(ind, 'CONTROL_MERGED_randomised_COR_RSP_TOTAL.set')}}...
%     );
% precomp_wrap(STUDY, 'alleeg', ALLEEG, 'erpim', 1 ...
%     , 'nlines', 1490, 'smoothing', 100 ...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration');



%% Sub-funcs
function fh = sbf_erpim2x2(eegA, eegB...
                        , smthA, smthB...
                        , startms...
                        , cond1, cond2...
                        , ttlA, ttlB...
                        , cname1, cname2...
                        , ttl, cax...
                        , sortype, sortfld) %#ok<DEFNU>
    %
    fh = figure('Units','pixels',...
            'Position', get(0,'ScreenSize') ./ 2,...
            'Color', [1.0, 1.0, 1.0]);
    % CTRL
    timing = [startms eegA.pnts eegA.srate];
    subplot(2, 2, 1)
    pop_erpimage(eegA, 1, cond1, [], [ttlA '-' cname1], smthA...
        , 1, sortype, [], sortfld, 'times', timing...
        , 'caxis', cax, 'cbar', 'off'...
        , 'noxlabel', 'on', 'img_trialax_label', 'Reaction Time') 
    subplot(2, 2, 2)
    pop_erpimage(eegA, 1, cond2, [], [ttlA '-' cname2], smthA...
        , 1, sortype, [], sortfld, 'times', timing...
        , 'caxis', cax, 'cbar', 'on', 'cbar_title', '\muV'...
        , 'noxlabel', 'on', 'img_trialax_label', 'Reaction Time')
    % ADHD
    timing = [startms eegB.pnts eegB.srate];
    subplot(2, 2, 3)
    pop_erpimage(eegB, 1, cond1, [], [ttlB '-' cname1], smthB...
        , 1, sortype, [], sortfld, 'times', timing...
        , 'caxis', cax, 'cbar', 'off', 'img_trialax_label', 'Reaction Time')
    subplot(2, 2, 4)
    pop_erpimage(eegB, 1, cond2, [], [ttlB '-' cname2], smthB...
        , 1, sortype, [], sortfld, 'times', timing...
        , 'caxis', cax, 'cbar', 'off', 'img_trialax_label', 'Reaction Time')
    %
    sgtitle(ttl)
end

%
function [plv, plvCI, plvTime] = sbf_get_plv(eeg, roi, millis, filtSpec, nbootci)

    t1 = max(millis(1) - filtSpec.order, eeg.xmin * 1000);
    t2 = min(millis(2) + filtSpec.order, eeg.xmax * 1000);
    calcseg = find(eeg.times >= t1, 1) : find(eeg.times <= t2, 1, 'last');
    plvTime = eeg.times(calcseg);

    bootrep = {'' sprintf(' (with 95%% CIs from %d bootstraps)', nbootci)};
    fprintf('%s%s for %s;\nCalc segment:%d..%dms to extract:%d..%dms\n'...
        , 'Calculating Phase-locking value', bootrep{(nbootci > 0) + 1}...
        , eeg.setname, t1, t2, millis(1), millis(2))

    % calculate PLV
    [calcCplv, CI] = eegPLV(eeg.data(roi, calcseg, :), eeg.srate, filtSpec...
                            , 'nbootci', nbootci);
    
    use1 = find(eeg.times >= millis(1), 1) - calcseg(1) + 1;
    use2 = ceil((((millis(2) - millis(1)) / 1000) * eeg.srate) + use1);
    plv = calcCplv(use1:use2, :, :);
    if ~isempty(CI)
        plvCI = CI(use1:use2, :, :, :);
    end
    plvTime = plvTime(use1:use2);
end

% TODO - FINISH CONVERTING TO STANDALONE, AND FIX TO LOOP OVER GIVEN WINDOWS
function [meanps, latps] = sbf_erpwin_test(CEEG, AEEG, ROI, tstwn)

    diffNA = NaN(1, AEEG.trials - CEEG.trials);
    wdwlab = {'2' '3'
            'pre' 'post'};

    trodes = ROI{1, r};
    tstwdw = arrayfun(@(x) find(CEEG.times >= x, 1), tstwn(c, :));
    
    if numel(trodes) > 1
        cN = squeeze(mean(CEEG.data(trodes, tstwdw(1):tstwdw(2), :)));
        aN = squeeze(mean(AEEG.data(trodes, tstwdw(1):tstwdw(2), :)));
        cP = squeeze(mean(CEEG.data(trodes, tstwdw(3):tstwdw(4), :)));
        aP = squeeze(mean(AEEG.data(trodes, tstwdw(3):tstwdw(4), :)));
    else
        cN = squeeze(CEEG.data(trodes, tstwdw(1):tstwdw(2), :));
        aN = squeeze(AEEG.data(trodes, tstwdw(1):tstwdw(2), :));
        cP = squeeze(CEEG.data(trodes, tstwdw(3):tstwdw(4), :));
        aP = squeeze(AEEG.data(trodes, tstwdw(3):tstwdw(4), :));
    end

    [pNaov, testTN, statsN] = anova1([[mean(cN) diffNA]; mean(aN)].', grp, 'off'); %#ok<*ASGLU>
    [hN, pNks, ks2statN] = kstest2(mean(cN), mean(aN));
    [pPaov, testTP, statsP] = anova1([[mean(cP) diffNA]; mean(aP)].', grp, 'off');
    [hP, pPks, ks2statP] = kstest2(mean(cP), mean(aP));
    testTN = [testTN {'KS-test'; pNks; hN; ks2statN}]; %#ok<*AGROW>
    testTP = [testTP {'KS-test'; pPks; hP; ks2statP}];

    svnm = ['anova-wdw-mean-N' wdwlab{c, 1} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTN(2:end, :), 'VariableNames', testTN(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))
    svnm = ['anova-wdw-mean-P' wdwlab{c, 2} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTP(2:end, :), 'VariableNames', testTP(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))

    [~, lCN] = min(cN); lCN = CEEG.times(lCN + tstwdw(1) - 1);
    [~, lAN] = min(aN); lAN = AEEG.times(lAN + tstwdw(1) - 1);
    [~, lCP] = min(cP); lCP = CEEG.times(lCP + tstwdw(3) - 1);
    [~, lAP] = min(aP); lAP = AEEG.times(lAP + tstwdw(3) - 1);

    [plNaov, testTlN, statslN] = anova1([[lCN diffNA]; lAN].', grp, 'off');
    [hlN, plNks, ks2statlN] = kstest2(lCN, lAN);
    [plPaov, testTlP, statslP] = anova1([[lCP diffNA]; lAP].', grp, 'off');
    [hlP, plPks, ks2statlP] = kstest2(lCP, lAP);
    testTlN = [testTlN {'KS-test'; plNks; hlN; ks2statlN}];
    testTlP = [testTlP {'KS-test'; plPks; hlP; ks2statlP}];

    svnm = ['anova-wdw-lat-N' wdwlab{c, 1} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTlN(2:end, :), 'VariableNames', testTlN(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))
    svnm = ['anova-wdw-lat-P' wdwlab{c, 2} upper(cnd{1, c}) '-' ROI{2, r}];
    writetable(cell2table(testTlP(2:end, :), 'VariableNames', testTlP(1, :))...
        , fullfile(oud, 'erpimXrt', 'stats', [svnm '.csv']))
    
    meanps = [pNks pPks];
    latps = [plNks plPks];
end