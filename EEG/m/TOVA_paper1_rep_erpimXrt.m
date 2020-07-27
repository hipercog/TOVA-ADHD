ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extended_anal/erpimXrt';
ROI = [7 15 19 21 23 28 36];
Pz = 19;
lbwh = get(0,'ScreenSize');
ctl = 'Control';
atl = 'ADHD';
ttl = 'ERPxRT';


%% ERPIM x RT
ctrlEEG = pop_loadset('filepath', ind, 'filename', 'CONTROL_MERGED_randomised_COR_RSP_TOTAL.set');
adhdEEG = pop_loadset('filepath', ind, 'filename', 'ADHD_MERGED_randomised_COR_RSP_TOTAL.set');


%% plot sorted by RT
fh = sbf_erpim2x2(ctrlEEG, adhdEEG, 100, 285, Pz, ROI, ctl, atl, 'Pz', 'ROI'...
                , ttl, [-8.1 8.1], {'cor_rsp'}, 'duration');
print(fh, '-dsvg', fullfile(oud, 'TOVA-ERPIMxRT'))
%% get structs for group-wise erp-images - NOT NEEDED?
% control_erpimXrt = std_erpimage(ctrlEEG, 'channels', 1:128 ...
%     , 'nlines', 1500, 'smoothing', 100, 'recompute', 'on'...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration'...
%     , 'fileout', fullfile(oud, cttl));
% adhdin_erpimXrt = std_erpimage(adhdEEG, 'channels', 1:128 ...
%     , 'nlines', 4100, 'smoothing', 285, 'recompute', 'on'...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration'...
%     , 'fileout', fullfile(oud, attl));


%% epoch aligned to RT
ix = [ctrlEEG.event(ismember({ctrlEEG.event.type}, 'cor_rsp')).duration] > 600;
cEEGrt = pop_epoch(pop_select(ctrlEEG, 'notrial', ix), {'RESP'}, [-0.6 0.4]);
ix = [adhdEEG.event(ismember({adhdEEG.event.type}, 'cor_rsp')).duration] > 600;
aEEGrt = pop_epoch(pop_select(adhdEEG, 'notrial', ix), {'RESP'}, [-0.6 0.4]);
ttl = 'ControlvsADHD-RTlocked-ERP';
tstwn = [215 280 285 350];


%% plot sorted by TARGET
fh = sbf_erpim2x2(cEEGrt, aEEGrt, 100, 285, Pz, ROI, ctl, atl, 'Pz', 'ROI'...
                , ttl, [-8.1 8.1], {'TRGT'}, 'latency');
print(fh, '-dsvg', fullfile(oud, 'TOVA-RTlocked-ERPIMxTarget'))


%% ERP of ADHD v CTRL aligned to RT
% erp, src, pnts, zeropt, srate, tkoffset, lgnd, ttl, savename
AvCerp = [ctap_get_erp(cEEGrt, 'loc', Pz); ctap_get_erp(aEEGrt, 'loc', Pz)];
AvCsrc = {ctap_get_erp(cEEGrt, 'loc', ROI, 'roi', false);...
          ctap_get_erp(aEEGrt, 'loc', ROI, 'roi', false)};
fh = ctap_plot_basic_erp(...
    AvCerp, 512, find(cEEGrt.times == 0), cEEGrt.srate...
    , 'waves', AvCsrc...
    , 'lgnd', {'Control' 'ADHD'}...
    , 'ttl', ttl...
    , 'vlines', tstwn);
print(fh, '-dsvg', fullfile(oud, ttl))


%% Topoplot the neg+pos peaks
inc = 22;
fh = figure('Position', [1 1 lbwh(4) lbwh(3)], 'Color', 'w');
plt = 1;

cax = [-3 1];
for w = 1:2:3
    for tx = tstwn(w):inc:tstwn(w + 1)

        % Plot the Control above
        subplot(4, 3, plt)
        data = mean(mean(cEEGrt.data(:, tx:tx + inc - 2, :), 3), 2);
        topoplot(data, cEEGrt.chanlocs...
            , 'maplimits', cax, 'electrodes', 'on'...
            , 'emarker', {'.','k',[],1}, 'emarker2', {ROI, 'o','m',3,1})
        title(sprintf('%sms\n%s', num2str(round(cEEGrt.times([tx tx + inc - 1]))), ctl))
        cb = colorbar('EastOutside'); title(cb, '\muV', 'FontSize', 9)
        
        % Plot the ADHD below
        subplot(4, 3, plt + 3)
        data = mean(mean(aEEGrt.data(:, tx:tx + inc - 2, :), 3), 2);
        topoplot(data, aEEGrt.chanlocs...
            , 'maplimits', cax, 'electrodes', 'on'...
            , 'emarker', {'.','k',[],1}, 'emarker2', {ROI, 'o','m',3,1})
        title(atl)
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


%% test negative and positive peaks around RT
%TODO - TEST EACH TRODE SEPARATELY, BONF-HOM ADJUSTED?
diffNA = NaN(1, aEEGrt.trials - cEEGrt.trials);

trodes = ROI;
if numel(trodes) > 1
    cN = squeeze(mean(cEEGrt.data(trodes, tstwn(1):tstwn(2), :)));
    aN = squeeze(mean(aEEGrt.data(trodes, tstwn(1):tstwn(2), :)));
    cP = squeeze(mean(cEEGrt.data(trodes, tstwn(3):tstwn(4), :)));
    aP = squeeze(mean(aEEGrt.data(trodes, tstwn(3):tstwn(4), :)));
else
    cN = squeeze(cEEGrt.data(trodes, tstwn(1):tstwn(2), :));
    aN = squeeze(aEEGrt.data(trodes, tstwn(1):tstwn(2), :));
    cP = squeeze(cEEGrt.data(trodes, tstwn(3):tstwn(4), :));
    aP = squeeze(aEEGrt.data(trodes, tstwn(3):tstwn(4), :));
end

[pN, testTN, statsN] = anova1([[mean(cN) diffNA]; mean(aN)].');
[hN, pN, ks2statN] = kstest2(mean(cN), mean(aN));
[pP, testTP, statsP] = anova1([[mean(cP) diffNA]; mean(aP)].');
[hP, pP, ks2statP] = kstest2(mean(cP), mean(aP));
testTN = [testTN {'KS-test'; pN; hN; ks2statN}];
testTP = [testTP {'KS-test'; pP; hP; ks2statP}];

writetable(cell2table(testTN(2:end, :), 'VariableNames', testTN(1, :))...
    , fullfile(oud, 'anova-wdw-mean-NpreRT-ROI.csv'))
writetable(cell2table(testTP(2:end, :), 'VariableNames', testTP(1, :))...
    , fullfile(oud, 'anova-wdw-mean-PpostRT-ROI.csv'))

[~, lCN] = min(cN); lCN = cEEGrt.times(lCN + tstwn(1) - 1);
[~, lAN] = min(aN); lAN = aEEGrt.times(lAN + tstwn(1) - 1);
[~, lCP] = min(cP); lCP = cEEGrt.times(lCP + tstwn(3) - 1);
[~, lAP] = min(aP); lAP = aEEGrt.times(lAP + tstwn(3) - 1);

[plN, testTlN, statslN] = anova1([[lCN diffNA]; lAN].');
[hlN, plN, ks2statlN] = kstest2(lCN, lAN);
[plP, testTlP, statslP] = anova1([[lCP diffNA]; lAP].');
[hlP, plP, ks2statlP] = kstest2(lCP, lAP);
testTlN = [testTlN {'KS-test'; plN; hlN; ks2statlN}];
testTlP = [testTlP {'KS-test'; plP; hlP; ks2statlP}];

writetable(cell2table(testTlN(2:end, :), 'VariableNames', testTlN(1, :))...
    , fullfile(oud, 'anova-wdw-lat-NpreRT-ROI.csv'))
writetable(cell2table(testTlP(2:end, :), 'VariableNames', testTlP(1, :))...
    , fullfile(oud, 'anova-wdw-lat-PpostRT-ROI.csv'))


%% 
% [STUDY ALLEEG] = std_editset( [], []...
%     , 'name', 'TOVA_adhdin_erpimXrt'...
%     , 'task', 'correct_response_locked'...
%     , 'updatedat','on', 'savedat','on', 'rmclust','on'...
%     , 'commands', {{'index', 1 ...
%         , 'load', fullfile(ind, 'ADHD_MERGED_randomised_COR_RSP_TOTAL.set')}}...
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
%     , 'commands', {{'index', 1 ...
%         , 'load', fullfile(ind, 'CONTROL_MERGED_randomised_COR_RSP_TOTAL.set')}}...
%     );
% precomp_wrap(STUDY, 'alleeg', ALLEEG, 'erpim', 1 ...
%     , 'nlines', 1490, 'smoothing', 100 ...
%     , 'sorttype', 'cor_rsp', 'sortfield', 'duration');



%% Sub-funcs
function fh = sbf_erpim2x2(eegA, eegB...
                        , smthA, smthB...
                        , cond1, cond2...
                        , ttlA, ttlB...
                        , cname1, cname2...
                        , ttl, cax...
                        , sortype, sortfld)
    %%
    fh = figure('Units','pixels',...
            'Position', get(0,'ScreenSize'),...
            'Color', [1.0, 1.0, 1.0]);
    %% CTRL
    subplot(2, 2, 1)
    pop_erpimage(eegA, 1, cond1, [], [ttlA '-' cname1], smthA...
        , 1, sortype, [], sortfld...
        , 'caxis', cax, 'cbar', 'off'...
        , 'noxlabel', 'on', 'img_trialax_label', 'Reaction Time')
    subplot(2, 2, 2)
    pop_erpimage(eegA, 1, cond2, [], [ttlA '-' cname2], smthA...
        , 1, sortype, [], sortfld...
        , 'caxis', cax, 'cbar', 'on', 'cbar_title', '\muV'...
        , 'noxlabel', 'on', 'img_trialax_label', 'Reaction Time')
    %% ADHD
    subplot(2, 2, 3)
    pop_erpimage(eegB, 1, cond1, [], [ttlB '-' cname1], smthB...
        , 1, sortype, [], sortfld...
        , 'caxis', cax, 'cbar', 'off', 'img_trialax_label', 'Reaction Time')
    subplot(2, 2, 4)
    pop_erpimage(eegB, 1, cond2, [], [ttlB '-' cname2], smthB...
        , 1, sortype, [], sortfld...
        , 'caxis', cax, 'cbar', 'off', 'img_trialax_label', 'Reaction Time')
    %%
    sgtitle(ttl)
end