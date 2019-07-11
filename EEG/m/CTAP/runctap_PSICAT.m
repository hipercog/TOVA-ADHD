function ERPS = runctap_PSICAT(varargin)
%% Linear CTAP script to clean CENT data
% 
% OPERATION STEPS
% # 1 --------------
% Install / download:
% --- HARDWARE ---
%   * Matlab R2018b or newer
%   * EEGLAB, latest version,
%       git clone https://github.com/sccn/eeglab.git
%   * Add Biosig 3.3.0 (the Biosemi .bdf-reader plugin) to EEGLAB,
%       https://sccn.ucsd.edu/wiki/EEGLAB_Extensions
%   * CTAP,
%       git clone https://github.com/bwrc/ctap.git
% --- DATA ---
%   * 69 files of EEG data in .bdf format, plus Presentation .log files
% # 2 --------------
% Set your working directory to CTAP root
% # 3 --------------
% Add EEGLAB and CTAP to your Matlab path
% # 4 --------------
% Set up directory to contain .bdf files, pass complete path to 'proj_root'
%
% Syntax:
%   runctap_PSICAT('proj_root', <some_path>, 'GIX', 1|2)
% 
% Varargin
% 
%   proj_root   string, valid path to the data's project base folder, so that:
%                       <proj_root>/project_PSICAT/PSICAT-data/<group>
%               Default: '/wrk/group/hipercog/' (dir on ukko2 server)
%   group       string, name of the group to analyse, ADHD or CTRL
%               Default: ADHD
%   qcERPloc    string, if not null, run QC grand average ERP on this electrode
%               Default: ''
%   set_select  [1 n] | 'all', keyword 'all' selects all stepSets, or use index
%               Default: 'all'
%   sbj_filt    [1 n] | 'all', keyword 'all' for all, or some index of subjects
%               Default: 'all'
%   overwrite   logical, overwrite existing output
%               Default: true
%   debug       logical, stop on error
%               Default: false



%% Setup MAIN parameters
% use ctapID to uniquely name the base folder of the output directory tree
ctapID = 'PSICAT';
%change this to change the protocol? 
%options: HELLO-GOODBYE, NEUROFEEDBACK, PSICAT, GESTALT, TOVA, VIGILANCE

% specify the file type of your data
data_type = '*.bdf';

%TODO - DEFINE PSICAT EVENTS FOR EPOCHING
%set of event names for epoching and export for event-related analysis
evrelevs = {''};


%% Parse input parameters
p = inputParser;

p.addParameter('proj_root', '/wrk/group/hipercog/', @ischar)
p.addParameter('group', 'ADHD', @(x) ismember(x, {'ADHD' 'CTRL'}))
p.addParameter('qcERPloc', '', @ischar)
p.addParameter('set_select', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('sbj_filt', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('overwrite', true, @islogical)
p.addParameter('debug', false, @islogical)

p.parse(varargin{:});
Arg = p.Results;


%% Create the CONFIGURATION struct
% First, define step sets & their parameters in sbf_cfg()
[Cfg, ctap_args] = sbf_cfg(Arg.proj_root, ctapID, Arg.group, evrelevs, Arg.set_select);

% Next, create measurement config (MC) based on folder, & select subject subset
Cfg = get_meas_cfg_MC(Cfg, Arg.proj_root...
                    , 'eeg_ext', data_type...
                    , 'sbj_filt', Arg.sbj_filt...
                    , 'subject', {'[0-9]{4}[CP]', 0, 5}...
                    , 'subjectnr', {'[0-9]{4}[CP]', 0, 4}...
                    , 'session', {ctapID}...
                    , 'measurement', {'intake'});

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
function [Cfg, out] = sbf_cfg(project_root_folder, ID, grp, erevts, runsets)


%% Define important directories and files
Cfg.id = ID;
Cfg.env.paths.projectRoot = fullfile(project_root_folder, ['project_' ID]);
Cfg.env.paths.ctapRoot = fullfile(Cfg.env.paths.projectRoot, 'ANALYSIS');
Cfg.env.paths.analysisRoot = fullfile(Cfg.env.paths.ctapRoot, grp);
%THIS ONE IS PROBABLY OBSOLETE - Cfg.MC already gives the file locations
Cfg.env.paths.dataRoot = fullfile(Cfg.env.paths.projectRoot, [ID '-data'], grp);


%% Define stuff related to Channel locations
Cfg.eeg.chanlocs = which('chanlocs128_pist.elp');
Cfg.eeg.reference = {'L_MASTOID' 'R_MASTOID'};%{'average'};
Cfg.eeg.heogChannelNames = {'HEOG1' 'HEOG2'};
Cfg.eeg.veogChannelNames = {'VEOG1' 'VEOG2'};


%% STEPSET 1 - Load and prepare
i = 1;
stepSet(i).funH = { @CTAP_load_data,...
                    @CTAP_load_events,...
                    @CTAP_select_evdata,...
                    @CTAP_resample_data,...
                    @CTAP_load_chanlocs,...
                    @CTAP_fir_filter,...
                    @CTAP_fir_filter,...
                    @CTAP_peek_data };
stepSet(i).id = [num2str(i) '_load'];

out.load_events = struct(...
    'method', 'handle',...
    'handle', @loadCENTevents,...
    'src', Cfg.env.paths.dataRoot);

out.select_evdata = struct(...
    'covertype', 'total',...
    'duration', [-2000 2000]);

out.load_chanlocs = struct(...
    'overwrite', true,...
    'index_match', false);
out.load_chanlocs.field = {{{'EXG1' 'EXG2'} 'type' 'ECG'}...
                           {{'EXG3' 'EXG4' 'EXG5' 'EXG6'} 'type' 'EOG'}...
                           {{'EXG7' 'EXG8'} 'type' 'REF'}...
                           {'GSR1' 'type' 'EDA'}...
                {{'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp'} 'type' 'NA'}...
                {'EXG1' 'labels' 'ECG1'} {'EXG2' 'labels' 'ECG2'}...
                {'EXG3' 'labels' 'HEOG1'} {'EXG4' 'labels' 'HEOG2'}...
                {'EXG5' 'labels' 'VEOG1'} {'EXG6' 'labels' 'VEOG2'}...
                {'EXG7' 'labels' 'L_MASTOID'} {'EXG8' 'labels' 'R_MASTOID'}};
out.load_chanlocs.tidy  = {{'type' 'NA'}};

out.fir_filter = struct(...
    'locutoff', {0.5 []},...
    'hicutoff', {[] 30},...
    'filtorder', {3380 226});

out.peek_data = struct(...
    'secs', [10 30],... %start few seconds after data starts
    'peekStats', true,... %get statistics for each peek!
    'overwrite', false,...
    'plotAllPeeks', false,...
    'savePeekData', false,...
    'savePeekICA', false);


%% STEPSET 2 - Manual bad chans, reref to mastoids, detect blinks, create ICA
i = i+1;
stepSet(i).funH = { @CTAP_detect_bad_channels,...%given bad channels
                    @CTAP_reject_data,...
                    @CTAP_interp_chan,...
                    @CTAP_reref_data,...
                    @CTAP_blink2event,...
                    @CTAP_run_ica };
stepSet(i).id = [num2str(i) '_setup_ICA'];

out.detect_bad_channels = struct(...
     'method', 'given',...
     'badChanCsv', fullfile(Cfg.env.paths.dataRoot, ['badchans_' grp '.txt']));

out.interp_chan = struct('missing_types', 'EEG');

out.blink2event = struct(...
    'classMethod', 'emgauss_asymmetric');

out.run_ica = struct(...
    'method', 'fastica',...
    'overwrite', true);
out.run_ica.channels = {'EEG' 'EOG'};


%% STEPSET 3 - Artefact correction
i = i+1;  %stepSet 3
stepSet(i).funH = { @CTAP_detect_bad_comps,... %FASTER for non-blinks
                    @CTAP_detect_bad_comps,... %detect blink related ICs
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_channels,...%bad channels by Mahalanobis
                    @CTAP_reject_data,...
                    @CTAP_interp_chan,...
                    @CTAP_epoch_data,...
                    @CTAP_detect_bad_epochs,...
                    @CTAP_reject_data,...
                    @CTAP_peek_data };
%MAYBEDO - ADD @CTAP_filter_blink_ica AFTER BLINK BAD ICs?
stepSet(i).id = [num2str(i) '_denoise_epoch'];

out.detect_bad_channels(2).method = 'maha_fast';
out.detect_bad_channels(2).factorVal = 2.8;

out.detect_bad_comps = struct(...
    'method', {'faster' 'recu_blink_tmpl'},...
    'match_measures', {{'m' 's' 'k' 'h'} ''},...
    'bounds', {[-2.75 2.75] []},...
    'match_logic', {@any 0});

out.epoch_data = struct(...
    'method', 'epoch',...
    'match',  'exact',...
    'timelim', [-250 500],...
    'evtype', erevts);%TODO - POPULATE THIS WITH NAMES OF EPOCHING EVENTS!!

out.detect_bad_epochs = struct(...
    'channels', {{'A1'}},...
    'method', 'eegthresh',...
    'uV_thresh', [-100 100]);


%% STEPSET 4 - Export
i = i + 1; %next stepSet
stepSet(i).funH = repmat({@CTAP_export_data}, 1, numel(erevts));
stepSet(i).id = [num2str(i) '_export'];
stepSet(i).save = false;

out.export_data = struct(...
    'type', 'hdf5',...
    'outdir', fullfile('exportRoot', 'HDF5_EXPORT'),...
    'lock_event', erevts);


%% Store to Cfg
Cfg.pipe.stepSets = stepSet; % return all step sets inside Cfg struct
% step sets to run, default: whole thing
if nargin == 5
    Cfg.pipe.runSets = runsets;
else
    Cfg.pipe.runSets = {stepSet(:).id};
end

end
