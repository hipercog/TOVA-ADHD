function ERPS = runctap_TOVA(varargin)
%% Linear CTAP script to clean CENT data
%
% OPERATION STEPS
%
% # 1 - Install / download:
% --- SOFTWARE ---
%   * Matlab R2018b or newer
%   * EEGLAB, latest version,
%       git clone https://github.com/sccn/eeglab.git
%   * Biosig toolbox (install through EEGLAB)
%   * CTAP,
%       git clone https://github.com/bwrc/ctap.git
%   * CENT analysis tools,
%       git clone https://github.com/zenBen/CENT-analysis.git
% --- DATA ---
%   * 43 files of EEG data in .bdf format
%
% # 2 - Set your working directory to CTAP root
%
% # 3 - Add EEGLAB and CTAP to your Matlab path
%
% # 4 - Set up directory to contain .bdf files, pass full path to 'proj_root'
%
% Syntax:
%   runctap_TOVA('proj_root', <some_path>, 'group', <ADHD|CTRL>)
%
% Varargin
%   user        string, an identifier string to append to output base dir
%               default, ''
%   proj_root   string, valid path to the data's project base folder, so that:
%                       <proj_root>/project_TOVA/TOVA-data/<group>
%               Default: '/wrk/group/hipercog/' (dir on ukko2 server)
%   group       string, name of the group to analyse, ADHD or CTRL
%               Default: ADHD
%   qcERPloc    string, if not null, run QC grand average ERP on this electrode
%               Default: ''
%   getsets     [1 n] | 'all', keyword 'all' selects all stepSets, or use index
%               Default: 'all'
%   sbj_filt    [1 n] | 'all', keyword 'all' for all, or some index of subjects
%               Default: 'all'
%   overwrite   logical, overwrite existing output
%               Default: true
%   debug       logical, stop on error
%               Default: false



%% Setup MAIN parameters
% use ctapID to uniquely name the base folder of the output directory tree
ctapID = 'TOVA';

% specify the file type of your data
data_type = '*.bdf';


%% Parse input parameters
p = inputParser;

p.addParameter('user', '', @ischar)
p.addParameter('proj_root', '/wrk/group/hipercog/', @ischar)
p.addParameter('group', 'Control', @(x) ismember(x, {'Control' 'Intake' 'Outtake'}))
p.addParameter('getsets', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('sbj_filt', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('qcERPloc', '', @ischar)
p.addParameter('overwrite', true, @islogical)
p.addParameter('debug', false, @islogical)

p.parse(varargin{:});
Arg = p.Results;


%% Create the CONFIGURATION struct
% First, define step sets & their parameters in sbf_cfg()
[Cfg, ctap_args] = sbf_cfg(Arg.proj_root...
                         , ctapID...
                         , Arg.user...
                         , Arg.group...
                         , Arg.getsets);

% Next, create measurement config (MC) based on folder, & select subject subset
Cfg = get_meas_cfg_MC(Cfg, Cfg.env.paths.dataRoot...
                    , 'eeg_ext', data_type...
                    , 'sbj_filt', Arg.sbj_filt...
                    , 'subject', {'[0-9]{4}[CP]', 0, 5}...
                    , 'subjectnr', {'[0-9]{4}[CP]', 0, 4}...
                    , 'session', {ctapID}...
                    , 'measurement', {Arg.group});

% Assign arguments to the selected functions, perform various checks
Cfg = ctap_auto_config(Cfg, ctap_args);


%% Run the pipe
tic
    CTAP_pipeline_looper(Cfg, 'debug', Arg.debug, 'overwrite', Arg.overwrite)
toc


%% Finally, obtain ERPs of known conditions from the processed data
if ~isempty(Arg.qcERPloc)
    try
        ERPS = ctap_manu2_oddball_erps(Cfg, 'loc_label', Arg.qcERPloc);
    catch ME
        if Arg.debug
            error('runctap_PSICAT:erps', '%s', ME.message)
        else
            warning('runctap_PSICAT:erps', '%s', ME.message)
        end
    end
else
    ERPS = [];
end

end

%% %%%%%%%%%%%%%%%%%%% CONFIGURE ANALYSIS PIPE!! %%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration Subfunction
function [Cfg, out] = sbf_cfg(proj_dir, ID, user, grp, runsets)


%% Define important directories and files
Cfg.id = ID;
Cfg.env.paths.projectRoot = fullfile(proj_dir, ['project_' ID]);
Cfg.env.paths.ctapRoot = fullfile(Cfg.env.paths.projectRoot, 'ANALYSIS');
Cfg.env.paths.analysisRoot = fullfile(Cfg.env.paths.ctapRoot, [grp  user]);
%THIS ONE IS PROBABLY OBSOLETE - Cfg.MC already gives the file locations
Cfg.env.paths.dataRoot = fullfile(Cfg.env.paths.projectRoot, 'data', grp);


%% Define stuff related to Channel locations
Cfg.eeg.chanlocs = which('chanlocs128_pist.elp');
Cfg.eeg.reference = {'L_MASTOID' 'R_MASTOID'};%{'average'};
Cfg.eeg.heogChannelNames = {'HEOG1' 'HEOG2'};
Cfg.eeg.veogChannelNames = {'VEOG1' 'VEOG2'};


%% LOAD AND RE-REF
i = 1; %stepSet 1
stepSet(i).funH = { @CTAP_load_data,...
                    @CTAP_load_chanlocs,...
                    @CTAP_load_events,...
                    @CTAP_reref_data};
stepSet(i).id = [num2str(i) '_LOAD_TOVA'];

out.load_data = struct(...
    'biosig', true);

out.load_chanlocs = struct(...
    'overwrite', true,...
    'index_match', false);
out.load_chanlocs.field = {{1:128 'type' 'EEG'}...
                           {{'EXG1' 'EXG2'} 'type' 'ECG'}...
                           {{'EXG3' 'EXG4' 'EXG5' 'EXG6'} 'type' 'EOG'}...
                           {{'EXG7' 'EXG8'} 'type' 'REF'}...
                           {'GSR1' 'type' 'EDA'}...
        {{'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp' 'Status'} 'type' 'NA'}...
        {'EXG1' 'labels' 'ECG1'} {'EXG2' 'labels' 'ECG2'}...
        {'EXG3' 'labels' 'HEOG1'} {'EXG4' 'labels' 'HEOG2'}...
        {'EXG5' 'labels' 'VEOG1'} {'EXG6' 'labels' 'VEOG2'}...
        {'EXG7' 'labels' 'L_MASTOID'} {'EXG8' 'labels' 'R_MASTOID'}};
out.load_chanlocs.tidy  = {{'type' 'NA'}};

out.load_events = struct(...
    'method', 'handle',...
    'handle', @parseTOVAevents);

out.reref_data = struct('reference', 'REF');


%% EXTRACT ANS SIGNALS
i = i+1;  %stepSet 2
stepSet(i).funH = { @CTAP_extract_signal};
stepSet(i).id = [num2str(i) '_EXTRACT'];
stepSet(i).save = false;
% extract data is only to be called once -> hence a cut here

out.extract_signal.types = {'ECG' 'EDA'};


%% FILTER
i = i+1;  %stepSet 3
stepSet(i).funH = { @CTAP_fir_filter,...
                    @CTAP_fir_filter };
stepSet(i).id = [num2str(i) '_FILTER'];
stepSet(i).srcID = '1_LOAD_TOVA'; %take source from step 1
% filtering takes ages -> hence a cut here

out.fir_filter = struct(...
    'locutoff', {2 []},...
    'hicutoff', {[] 45},...
    'filtorder', {3380 226});


%% SELECT, BAD CHANNELS, ICA
i = i+1;  %stepSet 4
stepSet(i).funH = { @CTAP_select_evdata,...
                    @CTAP_peek_data,...
                    @CTAP_detect_bad_channels,...%given bad channels
                    @CTAP_detect_bad_channels,... %faster
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_segments,... &+-80uV
                    @CTAP_reject_data,...
                    @CTAP_run_ica};
stepSet(i).id = [num2str(i) '_SUBSET_PEEK_CHANS_ICA'];
% ICA can take ages -> hence a cut here

out.select_evdata = struct(...
    'evtype', 'TOVAtest',...
    'covertype', 'own');
%     'duration', [-2000 2000]);

out.peek_data = struct(...
    'numpeeks', 20,...
    'secs', [10 30],... %start few seconds after data starts
    'peekStats', true,... %get statistics for each peek!
    'overwrite', false,...
    'savePeekData', false,...
    'savePeekICA', false);

%For all 'detect_bad_x' you should specify at least a method
out.detect_bad_channels = struct(...
     'method', 'given',...
     'badChanCsv', fullfile(Cfg.env.paths.dataRoot, ['badchans_' grp '.txt']));

out.detect_bad_channels(2).method = 'faster';
out.detect_bad_channels(2).refChannel = 'C21'; %should be 'C21' alias 'Fz'

out.detect_bad_segments = struct(...
    'amplitudeTh', [-80, 80]); %in muV [-100, 100]

out.run_ica = struct(...
    'method', 'fastica');
out.run_ica.channels = {'EEG' 'EOG'};


%% DETECT, REJECT BAD ICs
i = i+1;  %stepSet 5
stepSet(i).funH = { @CTAP_blink2event,...
                    @CTAP_detect_bad_comps,... %blinks
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_comps,... %ADJUST
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_comps,... %faster
                    @CTAP_reject_data};
stepSet(i).id = [num2str(i) '_IC_SEG_CORRECTION'];

out.detect_bad_comps = struct(...
    'method', {'blink_template' 'adjust' 'faster'},...
    'adjustarg', {'' 'horiz' ''});
% out.detect_bad_comps = struct(...
%     'method', 'recu_blink_tmpl');


%% INTERPOLATE BAD CHANNELS AND FINAL PEEK
i = i+1;  %stepSet 6
stepSet(i).funH = { @CTAP_interp_chan,...
                    @CTAP_peek_data};
stepSet(i).id = [num2str(i) '_INTERP'];

out.interp_chan = struct('missing_types', 'EEG');



%% Store to Cfg
Cfg.pipe.stepSets = stepSet; % return all step sets inside Cfg struct
% step sets to run, default: whole thing
if exist('runsets', 'var')
    Cfg.pipe.runSets = runsets;
else
    Cfg.pipe.runSets = {stepSet(:).id};
end

end
