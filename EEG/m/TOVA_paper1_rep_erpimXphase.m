ind = '/home/bcowley/Benslab/CENT/project_TOVA/TOVA-data/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extended_anal/erpimXphase';
ROI = [7 15 19 21 23 28 36];
Pz = 19;
ctl = 'Control';
atl = 'ADHD';
ttl = 'ERP x \alphaphase';
lbwh = get(0,'ScreenSize');


%% ERPIM x pre-stim alpha phase
cEEGh1 = pop_loadset('filepath', ind, 'filename', 'MERGED_H1_CORRECT_RANDOMISED_CONTROL.set');
cEEGh2 = pop_loadset('filepath', ind, 'filename', 'MERGED_H2_CORRECT_RANDOMISED_CONTROL.set');
aEEGh1 = pop_loadset('filepath', ind, 'filename', 'merged_h1_correct_randomised_adhd.set');
aEEGh2 = pop_loadset('filepath', ind, 'filename', 'merged_h2_correct_randomised_adhd.set');


%% testi
TRD = Pz;
NM = 'Pz';
cax = [-11 11];
fh = figure('Units','pixels', 'Position', [1 1 lbwh(4) lbwh(3)], 'Color', 'w');
subplot(4, 1, 1)
pop_erpimage(...
    cEEGh1, 1, TRD, [], [ctl ' - H1 - ' ttl ' - ' NM], 100, 1, '', [], []...
    , 'phasesort', [-80 25 10]...
    , 'caxis', cax, 'cbar', 'on', 'cbar_title', '\muV', 'noxlabel', 'on')
% print(fh, '-dsvg', fullfile(oud, [ctl '-H1-erpimXphase-Pz']))

subplot(4, 1, 2)
pop_erpimage(...
    cEEGh2, 1, TRD, [], [ctl ' - H2 - ' ttl ' - ' NM], 100, 1, '', [], []...
    , 'phasesort', [-80 25 10], 'cbar', 'off', 'noxlabel', 'on')
% print(fh, '-dsvg', fullfile(oud, [ctl '-H2-erpimXphase-Pz']))

subplot(4, 1, 3)
pop_erpimage(...
    aEEGh1, 1, TRD, [], [atl ' - H1 - ' ttl ' - ' NM], 285, 1, '', [], []...
    , 'phasesort', [-80 25 10], 'cbar', 'off', 'noxlabel', 'on')
% print(fh, '-dsvg', fullfile(oud, [atl '-H1-erpimXphase-Pz']))

subplot(4, 1, 4)
pop_erpimage(...
    aEEGh2, 1, TRD, [], [atl ' - H2 - ' ttl ' - ' NM], 285, 1, '', [], []...
    , 'phasesort', [-80 25 10], 'cbar', 'off')
% print(fh, '-dsvg', fullfile(oud, [atl '-H2-erpimXphase-Pz']))
print(fh, '-dsvg', fullfile(oud, ['TOVA-erpimXphase-' NM]))


%% Load EEGLAB study for pre-stim phase work
% [STUDY ALLEEG] = std_editset( [], []...
%     , 'name', 'TOVA_erpimXphase'...
%     , 'task', 'correct_rand_h1h2'...
%     , 'updatedat','on', 'savedat','on', 'rmclust','on'...
%     , 'commands'...
%         , {{'index',1,'load',fullfile(ind, 'merged_h1_correct_randomised_adhd.set'),'run',1,'group','1'}...
%         ,  {'index',2,'load',fullfile(ind, 'merged_h2_correct_randomised_adhd.set'),'run',2,'group','1'}...
%         ,  {'index',3,'load',fullfile(ind, 'MERGED_H1_CORRECT_RANDOMISED_CONTROL.set'),'run',1,'group','2'}...
%         ,  {'index',4,'load',fullfile(ind, 'MERGED_H2_CORRECT_RANDOMISED_CONTROL.set'),'run',2,'group','2'}}...
%     );

