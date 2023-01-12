% SCRIPT TO LOAD AND CLEAN TOVA 1 EEG data
pop_editoptions('option_storedisk ', 0)
%% INITIALISE
ind = '/home/bcowley/Benslab/CENT/project_TOVA/data/paper1';
oud = '/home/bcowley/Benslab/CENT/project_TOVA/ANALYSIS/paper1_extanal';
% oud = '/wrk/group/hipercog/project_TOVA/ANALYSIS/paper1/PLV';

% 16 channels spread around the scalp
ROI = {[7 5 4 32 36 17 21 30] [68 77 85 90 100 78 83 91]
        'parietal' 'frontal'};
% chs = [ROI{1, :}];
% biosemi = cellfun(@bionum2nimi, ROI(1, :), 'Un', 0)
% ten20 = cellfun(@biosemi1020, ROI(1, :), 'Un', 0);

grp = {'ctrl' 'adhd'};
Grup = {'Control' 'ADHD'};
cnd = {'rsp' 'inb'};
Cond = {'Response' 'Inhibition'};
which_g = 1:2;
which_c = 1:2;
which_h = '12';
which_r = 1:2;

lbwh = get(0,'ScreenSize');
clrmp = {'greytight'
'whatjet(''what'', [0.9 0.9 0.9], ''stops'', [0 0.15 0.25 0.1 0.1 0.25 0.15])'};
cmap = eval(clrmp{2, 1});


%% LOAD DATA IF EXISTS
% load(fullfile(oud, 'TOVA1eeg.mat'), 'aEEGinb', 'aEEGrsp', 'cEEGinb', 'cEEGrsp')

%% ...by HALVES
% load(fullfile(oud, 'TOVA1eeg.mat'), 'aEEGinb_H1', 'aEEGinb_H2', 'aEEGrsp_H1', 'aEEGrsp_H2', 'cEEGinb_H1', 'cEEGinb_H2', 'cEEGrsp_H1', 'cEEGrsp_H2')


%% FIND CORRECT DATA
% CTRL - RESPONSE
c_rsp_fs = dir(fullfile(ind, 'CTRL_COR_RSP', '*_TOVA.set'));
% CTRL - INHIBITION
c_inb_fs = dir(fullfile(ind, 'CTRL_COR_INB', '*_TOVA.set'));
% CTRL - RESPONSE H1
c_rsp_h1 = dir(fullfile(ind, 'CTRL_COR_RSP', 'H1', '*_TOVA.set'));
% CTRL - INHIBITION H1
c_inb_h1 = dir(fullfile(ind, 'CTRL_COR_INB', 'H1', '*_TOVA.set'));
% CTRL - RESPONSE H2
c_rsp_h2 = dir(fullfile(ind, 'CTRL_COR_RSP', 'H2', '*_TOVA.set'));
% CTRL - INHIBITION H2
c_inb_h2 = dir(fullfile(ind, 'CTRL_COR_INB', 'H2', '*_TOVA.set'));

% ADHD - RESPONSE
a_rsp_fs = dir(fullfile(ind, 'ADHD_COR_RSP', '*_TOVA.set'));
% ADHD - INHIBITION
a_inb_fs = dir(fullfile(ind, 'ADHD_COR_INB', '*_TOVA.set'));
% ADHD - RESPONSE H1
a_rsp_h1 = dir(fullfile(ind, 'ADHD_COR_RSP', 'H1', '*_TOVA.set'));
% ADHD - INHIBITION H1
a_inb_h1 = dir(fullfile(ind, 'ADHD_COR_INB', 'H1', '*_TOVA.set'));
% ADHD - RESPONSE H2
a_rsp_h2 = dir(fullfile(ind, 'ADHD_COR_RSP', 'H2', '*_TOVA.set'));
% ADHD - INHIBITION H2
a_inb_h2 = dir(fullfile(ind, 'ADHD_COR_INB', 'H2', '*_TOVA.set'));


%% INTERP BAD CHANNEL LIST
c_intrp_lst = {
    [80 81 89 93 100]
    NaN 
    [44 61 127] 
    [1 84 111]
    [80 81 89 93 103]
    NaN 
    65 
    [80 81 84 93]
    [93 100]
    33
    NaN
    [80 86 93]
    [69 80 81 90 92 93]
    NaN
    NaN
    [72 76 80 81 92 93 102]
    [80 81 92 93]
    44};
T_Ctrl = table(c_intrp_lst, 'RowNames', {
    '2001'
    '2002'
    '2003'
    '2004'
    '2005'
    '2006'
    '2007'
    '2008'
    '2010'
    '2011'
    '2012'
    '2016'
    '2017'
    '2018'
    '2019'
    '2020'
    '2021'
    '2022'});

a_intrp_lst = {
    [80 93]
    [1 65]
    NaN
    NaN
    [40 107 119]
    NaN
    [96 101]
    NaN
    [42 46 56 62]
    NaN
    [58 111]
    NaN
    [76 77]
    93
    21
    NaN
    [76 77 91 93]
    25
    58
    [59 64 101 113]
    74
    NaN
    NaN
    NaN
    [30 39 40 41 57 58 62 121]
    NaN
    [20 21 83 86 89]
    NaN
    [89 90 94 100]
    [30 33 53 54 55]
    [120 127 128]
    NaN
    [31 42 46 48]
    96
    NaN
    NaN
    [1 2]
    33
    [89 90]
    NaN
    33
    [76 80]
    NaN
    NaN
    52
    NaN
    NaN
    15
    [62 76 78]};
T_ADHD = table(a_intrp_lst, 'RowNames', {
    '1001'
    '1002'
    '1007'
    '1012'
    '1015'
    '1017'
    '1019'
    '1023'
    '1024'
    '1025'
    '1026'
    '1027'
    '1029'
    '1033'
    '1034'
    '1035'
    '1036'
    '1041'
    '1042'
    '1045'
    '1048'
    '1050'
    '1051'
    '1052'
    '1053'
    '1054'
    '1055'
    '1056'
    '1057'
    '1058'
    '1059'
    '1060'
    '1061'
    '1062'
    '1064'
    '1065'
    '1066'
    '1067'
    '1068'
    '1070'
    '1072'
    '1073'
    '1074'
    '1076'
    '1077'
    '1078'
    '1079'
    '1080'
    '1081'});


%% LOAD CORRECT RESPONSE DATA - CTRL
% [cSTDYrsp, cEEGrsp] = std_editset([], []...
%         , 'name', 'CTRL_COR_RSP'...
%         , 'task', 'TOVA correct response, CONTROL group'...
%         , 'filename', 'CTRL_COR_RSP.study'...
%         , 'filepath', fullfile(ind, 'CTRL_COR_RSP')...
%         , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
%                 , num2cell(1:numel(c_rsp_fs))...
%                 , fullfile({c_rsp_fs.folder}, {c_rsp_fs.name})...
%                 , cellfun(@(x) x(1:4), {c_rsp_fs.name}, 'Un', 0)...
%                 , 'Un', 0));
% % .. AND INTERPOLATE
% for i = 1:numel(cEEGrsp)
%     idx = ismember(T_Ctrl.Properties.RowNames, cEEGrsp(i).subject(1:4));
%     if ~isnan(T_Ctrl.c_intrp_lst{idx})
%         cEEGrsp(i) = pop_interp(cEEGrsp(i), T_Ctrl.c_intrp_lst{idx});
%     end
% end


%% ...BY HALVES
[cSTDYrsp_H1, cEEGrsp_H1] = std_editset([], []...
        , 'name', 'ADHD_COR_RSP_H1'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_RSP_H1.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_RSP', 'H1')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(c_rsp_h1))...
                , fullfile({c_rsp_h1.folder}, {c_rsp_h1.name})...
                , cellfun(@(x) x(1:4), {c_rsp_h1.name}, 'Un', 0)...
                , 'Un', 0));
% .. AND INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(cEEGrsp_H1)
    idx = ismember(T_Ctrl.Properties.RowNames, cEEGrsp_H1(i).subject(1:4));
    if ~isnan(T_Ctrl.c_intrp_lst{idx})
        cEEGrsp_H1(i) = pop_interp(cEEGrsp_H1(i), T_Ctrl.c_intrp_lst{idx});
    end
    cEEGrsp_H1(i) = pop_select(cEEGrsp_H1(i), 'notrial'...
        , isoutlier(squeeze(var(cEEGrsp_H1(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >50% TRIALS
[cSTDYrsp_H1, cEEGrsp_H1] = std_editset(cSTDYrsp_H1, cEEGrsp_H1, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([cEEGrsp_H1.trials] < 36)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(cEEGrsp_H1)
    shufdx = shuffle(1:cEEGrsp_H1(i).trials);
    cEEGrsp_H1(i) = pop_select(cEEGrsp_H1(i)...
                , 'trial', shufdx(1:min([cEEGrsp_H1.trials])));
end

[cSTDYrsp_H2, cEEGrsp_H2] = std_editset([], []...
        , 'name', 'ADHD_COR_RSP_H2'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_RSP_H2.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_RSP', 'H2')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(c_rsp_h2))...
                , fullfile({c_rsp_h2.folder}, {c_rsp_h2.name})...
                , cellfun(@(x) x(1:4), {c_rsp_h2.name}, 'Un', 0)...
                , 'Un', 0));
% .. AND INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(cEEGrsp_H2)
    idx = ismember(T_Ctrl.Properties.RowNames, cEEGrsp_H2(i).subject(1:4));
    if ~isnan(T_Ctrl.c_intrp_lst{idx})
        cEEGrsp_H2(i) = pop_interp(cEEGrsp_H2(i), T_Ctrl.c_intrp_lst{idx});
    end
    cEEGrsp_H2(i) = pop_select(cEEGrsp_H2(i), 'notrial'...
        , isoutlier(squeeze(var(cEEGrsp_H2(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >25% TRIALS
[cSTDYrsp_H2, cEEGrsp_H2] = std_editset(cSTDYrsp_H2, cEEGrsp_H2, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([cEEGrsp_H2.trials] < 72)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(cEEGrsp_H2)
    shufdx = shuffle(1:cEEGrsp_H2(i).trials);
    cEEGrsp_H2(i) = pop_select(cEEGrsp_H2(i)...
                , 'trial', shufdx(1:min([cEEGrsp_H2.trials])));
end


%% LOAD CORRECT INHIBITION DATA - CTRL
% [cSTDYinb, cEEGinb] = std_editset([], []...
%         , 'name', 'CTRL_COR_INB'...
%         , 'task', 'TOVA correct inhibition, CONTROL group'...
%         , 'filename', 'CTRL_COR_INB.study'...
%         , 'filepath', fullfile(ind, 'CTRL_COR_INB')...
%         , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
%                 , num2cell(1:numel(c_inb_fs))...
%                 , fullfile({c_inb_fs.folder}, {c_inb_fs.name})...
%                 , cellfun(@(x) x(1:4), {c_inb_fs.name}, 'Un', 0)...
%                 , 'Un', 0));
% % .. AND INTERPOLATE
% for i = 1:numel(cEEGinb)
%     idx = ismember(T_Ctrl.Properties.RowNames, cEEGinb(i).subject(1:4));
%     if ~isnan(T_Ctrl.c_intrp_lst{idx})
%         cEEGinb(i) = pop_interp(cEEGinb(i), T_Ctrl.c_intrp_lst{idx});
%     end
% end


%% ...BY HALVES
[cSTDYinb_H1, cEEGinb_H1] = std_editset([], []...
        , 'name', 'ADHD_COR_INB_H1'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_INB_H1.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_INB', 'H1')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(c_inb_h1))...
                , fullfile({c_inb_h1.folder}, {c_inb_h1.name})...
                , cellfun(@(x) x(1:4), {c_inb_h1.name}, 'Un', 0)...
                , 'Un', 0));
% .. AND INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(cEEGinb_H1)
    idx = ismember(T_Ctrl.Properties.RowNames, cEEGinb_H1(i).subject(1:4));
    if ~isnan(T_Ctrl.c_intrp_lst{idx})
        cEEGinb_H1(i) = pop_interp(cEEGinb_H1(i), T_Ctrl.c_intrp_lst{idx});
    end
    cEEGinb_H1(i) = pop_select(cEEGinb_H1(i), 'notrial'...
        , isoutlier(squeeze(var(cEEGinb_H1(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >25% TRIALS
[cSTDYinb_H1, cEEGinb_H1] = std_editset(cSTDYinb_H1, cEEGinb_H1, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([cEEGinb_H1.trials] < 72)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(cEEGinb_H1)
    shufdx = shuffle(1:cEEGinb_H1(i).trials);
    cEEGinb_H1(i) = pop_select(cEEGinb_H1(i)...
                , 'trial', shufdx(1:min([cEEGinb_H1.trials])));
end
            
[cSTDYinb_H2, cEEGinb_H2] = std_editset([], []...
        , 'name', 'ADHD_COR_INB_H2'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_INB_H2.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_INB', 'H2')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(c_inb_h2))...
                , fullfile({c_inb_h2.folder}, {c_inb_h2.name})...
                , cellfun(@(x) x(1:4), {c_inb_h2.name}, 'Un', 0)...
                , 'Un', 0));
% .. INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(cEEGinb_H2)
    idx = ismember(T_Ctrl.Properties.RowNames, cEEGinb_H2(i).subject(1:4));
    if ~isnan(T_Ctrl.c_intrp_lst{idx})
        cEEGinb_H2(i) = pop_interp(cEEGinb_H2(i), T_Ctrl.c_intrp_lst{idx});
    end
    cEEGinb_H2(i) = pop_select(cEEGinb_H2(i), 'notrial'...
        , isoutlier(squeeze(var(cEEGinb_H2(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >50% TRIALS
[cSTDYinb_H2, cEEGinb_H2] = std_editset(cSTDYinb_H2, cEEGinb_H2, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([cEEGinb_H2.trials] < 36)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(cEEGinb_H2)
    shufdx = shuffle(1:cEEGinb_H2(i).trials);
    cEEGinb_H2(i) = pop_select(cEEGinb_H2(i)...
                , 'trial', shufdx(1:min([cEEGinb_H2.trials])));
end


%% LOAD CORRECT RESPONSE DATA - ADHD
% [aSTDYrsp, aEEGrsp] = std_editset([], []...
%         , 'name', 'ADHD_COR_RSP'...
%         , 'task', 'TOVA correct response, ADHD group'...
%         , 'filename', 'ADHD_COR_RSP.study'...
%         , 'filepath', fullfile(ind, 'ADHD_COR_RSP')...
%         , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
%                 , num2cell(1:numel(a_rsp_fs))...
%                 , fullfile({a_rsp_fs.folder}, {a_rsp_fs.name})...
%                 , cellfun(@(x) x(1:4), {a_rsp_fs.name}, 'Un', 0)...
%                 , 'Un', 0));
% % .. AND INTERPOLATE
% for i = 1:numel(aEEGrsp)
%     idx = ismember(T_ADHD.Properties.RowNames, aEEGrsp(i).subject(1:4));
%     if ~isnan(T_ADHD.a_intrp_lst{idx})
%         aEEGrsp(i) = pop_interp(aEEGrsp(i), T_ADHD.a_intrp_lst{idx});
%     end
% end


%% ...BY HALVES
[aSTDYrsp_H1, aEEGrsp_H1] = std_editset([], []...
        , 'name', 'ADHD_COR_RSP_H1'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_RSP_H1.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_RSP', 'H1')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(a_rsp_h1))...
                , fullfile({a_rsp_h1.folder}, {a_rsp_h1.name})...
                , cellfun(@(x) x(1:4), {a_rsp_h1.name}, 'Un', 0)...
                , 'Un', 0));
% .. INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(aEEGrsp_H1)
    idx = ismember(T_ADHD.Properties.RowNames, aEEGrsp_H1(i).subject(1:4));
    if ~isnan(T_ADHD.a_intrp_lst{idx})
        aEEGrsp_H1(i) = pop_interp(aEEGrsp_H1(i), T_ADHD.a_intrp_lst{idx});
    end
    aEEGrsp_H1(i) = pop_select(aEEGrsp_H1(i), 'notrial'...
        , isoutlier(squeeze(var(aEEGrsp_H1(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >50% TRIALS
[aSTDYrsp_H1, aEEGrsp_H1] = std_editset(aSTDYrsp_H1, aEEGrsp_H1, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([aEEGrsp_H1.trials] < 36)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(aEEGrsp_H1)
    shufdx = shuffle(1:aEEGrsp_H1(i).trials);
    aEEGrsp_H1(i) = pop_select(aEEGrsp_H1(i)...
                , 'trial', shufdx(1:min([aEEGrsp_H1.trials])));
end

[aSTDYrsp_H2, aEEGrsp_H2] = std_editset([], []...
        , 'name', 'ADHD_COR_RSP_H2'...
        , 'task', 'TOVA correct response, ADHD group'...
        , 'filename', 'ADHD_COR_RSP_H2.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_RSP', 'H2')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(a_rsp_h2))...
                , fullfile({a_rsp_h2.folder}, {a_rsp_h2.name})...
                , cellfun(@(x) x(1:4), {a_rsp_h2.name}, 'Un', 0)...
                , 'Un', 0));
% .. INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(aEEGrsp_H2)
    idx = ismember(T_ADHD.Properties.RowNames, aEEGrsp_H2(i).subject(1:4));
    if ~isnan(T_ADHD.a_intrp_lst{idx})
        aEEGrsp_H2(i) = pop_interp(aEEGrsp_H2(i), T_ADHD.a_intrp_lst{idx});
    end
    aEEGrsp_H2(i) = pop_select(aEEGrsp_H2(i), 'notrial'...
        , isoutlier(squeeze(var(aEEGrsp_H2(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >25% TRIALS
[aSTDYrsp_H2, aEEGrsp_H2] = std_editset(aSTDYrsp_H2, aEEGrsp_H2, 'commands',...
cellfun(@(x) {'remove' x}, num2cell(find([aEEGrsp_H2.trials] < 72)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(aEEGrsp_H2)
    shufdx = shuffle(1:aEEGrsp_H2(i).trials);
    aEEGrsp_H2(i) = pop_select(aEEGrsp_H2(i)...
                , 'trial', shufdx(1:min([aEEGrsp_H2.trials])));
end


%% LOAD CORRECT INHIBITION DATA - ADHD
% [aSTDYinb, aEEGinb] = std_editset([], []...
%         , 'name', 'ADHD_COR_INB'...
%         , 'task', 'TOVA correct inhibition, ADHD group'...
%         , 'filename', 'ADHD_COR_INB.study'...
%         , 'filepath', fullfile(ind, 'ADHD_COR_INB')...
%         , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
%                 , num2cell(1:numel(a_inb_fs))...
%                 , fullfile({a_inb_fs.folder}, {a_inb_fs.name})...
%                 , cellfun(@(x) x(1:4), {a_inb_fs.name}, 'Un', 0)...
%                 , 'Un', 0));
% % .. AND INTERPOLATE
% for i = 1:numel(aEEGinb)
%     idx = ismember(T_ADHD.Properties.RowNames, aEEGinb(i).subject(1:4));
%     if ~isnan(T_ADHD.a_intrp_lst{idx})
%         aEEGinb(i) = pop_interp(aEEGinb(i), T_ADHD.a_intrp_lst{idx});
%     end
% end


%% ...BY HALVES
[aSTDYinb_H1, aEEGinb_H1] = std_editset([], []...
        , 'name', 'ADHD_COR_INB_H1'...
        , 'task', 'TOVA correct inhibition, ADHD group'...
        , 'filename', 'ADHD_COR_INB_H1.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_INB', 'H1')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(a_inb_h1))...
                , fullfile({a_inb_h1.folder}, {a_inb_h1.name})...
                , cellfun(@(x) x(1:4), {a_inb_h1.name}, 'Un', 0)...
                , 'Un', 0));
% .. INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(aEEGinb_H1)
    idx = ismember(T_ADHD.Properties.RowNames, aEEGinb_H1(i).subject(1:4));
    if ~isnan(T_ADHD.a_intrp_lst{idx})
        aEEGinb_H1(i) = pop_interp(aEEGinb_H1(i), T_ADHD.a_intrp_lst{idx});
    end
    aEEGinb_H1(i) = pop_select(aEEGinb_H1(i), 'notrial'...
        , isoutlier(squeeze(var(aEEGinb_H1(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >25% TRIALS
[aSTDYinb_H1, aEEGinb_H1] = std_editset(aSTDYinb_H1, aEEGinb_H1, 'commands',...
cellfun(@(x) {'remove' x}, num2cell(find([aEEGinb_H1.trials] < 72)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(aEEGinb_H1)
    shufdx = shuffle(1:aEEGinb_H1(i).trials);
    aEEGinb_H1(i) = pop_select(aEEGinb_H1(i)...
                , 'trial', shufdx(1:min([aEEGinb_H1.trials])));
end

[aSTDYinb_H2, aEEGinb_H2] = std_editset([], []...
        , 'name', 'ADHD_COR_INB_H2'...
        , 'task', 'TOVA correct inhibition, ADHD group'...
        , 'filename', 'ADHD_COR_INB_H2.study'...
        , 'filepath', fullfile(ind, 'ADHD_COR_INB', 'H2')...
        , 'commands', cellfun(@(ix, nm, no) {'index' ix 'load' nm 'subject' no}...
                , num2cell(1:numel(a_inb_h2))...
                , fullfile({a_inb_h2.folder}, {a_inb_h2.name})...
                , cellfun(@(x) x(1:4), {a_inb_h2.name}, 'Un', 0)...
                , 'Un', 0));
% .. INTERPOLATE, REMOVE OUTLIERS
for i = 1:numel(aEEGinb_H2)
    idx = ismember(T_ADHD.Properties.RowNames, aEEGinb_H2(i).subject(1:4));
    if ~isnan(T_ADHD.a_intrp_lst{idx})
        aEEGinb_H2(i) = pop_interp(aEEGinb_H2(i), T_ADHD.a_intrp_lst{idx});
    end
    aEEGinb_H2(i) = pop_select(aEEGinb_H2(i), 'notrial'...
        , isoutlier(squeeze(var(aEEGinb_H2(i).data, 0, [1 2]))'));
end
% .. KEEP ONLY THOSE WITH >50% TRIALS
[aSTDYinb_H2, aEEGinb_H2] = std_editset(aSTDYinb_H2, aEEGinb_H2, 'commands',...
 cellfun(@(x) {'remove' x}, num2cell(find([aEEGinb_H2.trials] < 36)), 'Un', 0));
% .. AND SELECT TRIAL SUBSET
for i = 1:numel(aEEGinb_H2)
    shufdx = shuffle(1:aEEGinb_H2(i).trials);
    aEEGinb_H2(i) = pop_select(aEEGinb_H2(i)...
                , 'trial', shufdx(1:min([aEEGinb_H2.trials])));
end


%% HANDLE MERGING?!
% LOAD KJ's PREMERGED CORRECT RESPONSE DATA - CTRL
% cEEGrsp0 = pop_loadset('filepath', ind...
%                  , 'filename', 'CONTROL_MERGED_randomised_COR_RSP_TOTAL.set');
% aEEG0 = pop_loadset('filepath', ind...
%                  , 'filename', 'ADHD_MERGED_randomised_COR_RSP_TOTAL.set');

% cEEGinb = pop_loadset('filepath', fullfile(ind, stdy), 'filename', fs(1).name);
% c_trials.(['trials_' cEEGinb.filename(1:5)]) = cEEGinb.trials;
% for i = 2:numel(fs)
%     eegtmp = pop_loadset('filepath', indir, 'filename', fs(i).name);
%     c_trials.(['trials_' eegtmp.filename(1:6)]) = eegtmp.trials;
%     if eegtmp.trials == 108
%         cEEGinb = pop_mergeset(cEEGinb, eegtmp);
%     end
% end

% % OLD WAY TO INTERP BAD CHANNELS FROM A MERGED DATASET
% for i = 1:108:cEEGrsp0.trials
%     cEEGtmp = pop_select(cEEGrsp0, 'trial', i:i+107);
%     if ~isnan(c_intrp_lst{1})
%         cEEGtmp = pop_interp(cEEGtmp, c_intrp_lst{1});
%         cEEGrsp.data(:, :, i:i+107) = cEEGtmp.data;
%     end
%     c_intrp_lst(1) = [];
% end

% OLD REMOVE BAD EPOCHS
% cEEGrsp = pop_select(cEEGrsp, 'notrial', [230 333 337 339 343 347 348 354 ...
%     356 359 364 369 371 373 395 398 403 408 413 424 567 719 1058 1111 1166]);
% c_rsp_
% 
% % Remove bad subject #24 + bad epochs from other subjects
% aEEG = pop_select(aEEG, 'notrial', [666 747 753 ...
%     1115 1171 1182 1513 1514 1664 1665 1693 1708 ...
%     2358 2373 2375 2593:2700 ...
%     3782 ...
%     4229 4240 4263 4298 4315]);


%% Split RESPONSE data by median response time, for..?
% Crts = eeg_getepochevent(cEEGrsp, {'cor_rsp'}, [], 'duration');
% idx = Crts < median(Crts);
% cEEGloRT = pop_select(cEEGrsp, 'trial', find(idx));
% cEEGhiRT = pop_select(cEEGrsp, 'trial', find(~idx));
% Arts = eeg_getepochevent(aEEG, {'cor_rsp'}, [], 'duration');
% idx = Arts < median(Arts);
% aEEGloRT = pop_select(aEEG, 'trial', find(idx));
% aEEGhiRT = pop_select(aEEG, 'trial', find(~idx));

clear ind idx i c_intrp_lst a_intrp_lst a_inb_* a_rsp_* c_inb_* c_rsp_* shufdx T_*