ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
% ind = '/wrk/group/hipercog/project_TOVA/ANALYSIS/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extanal';
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

%% Split data by median response time, for..?
Crts = eeg_getepochevent(cEEG, {'cor_rsp'}, [], 'duration');
idx = Crts < median(Crts);
cEEGloRT = pop_select(cEEG, 'trial', find(idx));
cEEGhiRT = pop_select(cEEG, 'trial', find(~idx));
Arts = eeg_getepochevent(aEEG, {'cor_rsp'}, [], 'duration');
idx = Arts < median(Arts);
aEEGloRT = pop_select(aEEG, 'trial', find(idx));
aEEGhiRT = pop_select(aEEG, 'trial', find(~idx));


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
            set(axh{3}.Children(5), 'FaceColor', [0.8 0.8 0.8])
            if plt < 5
                set(axh{3}.Children(2), 'Color', 'r')
            end
            if any(plt == [4 8])
                ylabel(axh{3}, '\muV'...
                    , 'FontWeight', 'bold'...
                    , 'Units', 'normalized'...
                    , 'Position', [1.15 1]...
                    , 'Rotation', 0)
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