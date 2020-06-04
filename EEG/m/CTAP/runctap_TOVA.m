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
%   * CTAP,
%       git clone https://github.com/bwrc/ctap.git
%   * CENT analysis tools,
%       git clone https://github.com/zenBen/CENT-analysis.git
% --- DATA ---
%   * 69 files of EEG data in .bdf format, plus Presentation .log files
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
ctapID = 'TOVA';

% specify the file type of your data
data_type = '*.bdf';


%% Parse input parameters
p = inputParser;

p.addParameter('user', '', @ischar)
p.addParameter('proj_root', '/wrk/group/hipercog/', @ischar)
p.addParameter('group', 'Control', @(x) ismember(x, {'Control' 'Intake' 'Outtake'}))
p.addParameter('set_select', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('sbj_filt', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('qcERPloc', 'A1', @ischar)
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
                         , Arg.set_select);

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
Cfg.env.paths.dataRoot = fullfile(Cfg.env.paths.projectRoot, [ID '-data'], grp);


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

%could also just call parseTOVAevents() directly in the pipe?
out.load_events = struct(...
    'method', 'handle',...
    'handle', @parseTOVAevents);
% out.load_events = struct(...
%     'method', 'handle',...
%     'handle', @loadCENTevents,...
%     'src', Cfg.env.paths.dataRoot,...
%     'log_ext', 'log');

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
% TODO - CHECK KJ'S EEG FILES TO FIND OUT WHAT FILTERING USED - IIR?!

out.fir_filter = struct(...
    'locutoff', {2 []},...
    'hicutoff', {[] 45},...
    'filtorder', {3380 226});


%% SELECT, BAD CHANNELS, ICA
i = i+1;  %stepSet 4
stepSet(i).funH = { @CTAP_select_evdata,...
                    @CTAP_peek_data,...
                    @CTAP_detect_bad_channels,... %faster
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_segments,... &+-80uV
                    @CTAP_reject_data,...
                    @CTAP_run_ica};
stepSet(i).id = [num2str(i) '_SUBSET_PEEK_CHANS_ICA'];
% ICA can take ages -> hence a cut here

out.select_evdata = struct(...
    'evtype', 'TOVAtest');
%     'covertype', 'total',...
%     'duration', [-2000 2000]);

out.peek_data = struct(...
    'secs', [10 30],... %start few seconds after data starts
    'peekStats', true,... %get statistics for each peek!
    'overwrite', false,...
    'plotAllPeeks', false,...
    'savePeekData', false,...
    'savePeekICA', false);

%For all 'detect_bad_x' you should specify at least a method
out.detect_bad_channels = struct(...
    'method', 'faster',...
    'refChannel', 'C21',... %should be 'C21' alias 'Fz'
    'channelType', {'EEG'});

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


%% INTERPOLATE BAD CHANNELS AND FINAL PEEK
i = i+1;  %stepSet 6
stepSet(i).funH = { @CTAP_interp_chan,...
                    @CTAP_peek_data};
stepSet(i).id = [num2str(i) '_INTERP'];

out.interp_chan = struct('missing_types', 'EEG');



% %% STEPSET 2 - Manual bad chans, reref to mastoids, detect blinks, create ICA
% i = i+1;
% stepSet(i).funH = { @CTAP_detect_bad_channels,...%given bad channels
%                     @CTAP_detect_bad_channels,...%bad channels by variance
%                     @CTAP_reject_data,...
%                     @CTAP_interp_chan,...
%                     @CTAP_reref_data,...
%                     @CTAP_blink2event,...
%                     @CTAP_detect_bad_segments,...
%                     @CTAP_reject_data,...
%                     @CTAP_run_ica };
% stepSet(i).id = [num2str(i) '_setup_ICA'];

% out.detect_bad_channels = struct(...
%      'method', 'given',...
%      'badChanCsv', fullfile(Cfg.env.paths.dataRoot, ['badchans_' grp '.txt']));
% 
% out.detect_bad_channels(2).method = 'variance';
% 
% out.blink2event = struct(...
%     'classMethod', 'emgauss_asymmetric');
% 
% out.detect_bad_segments = struct(...
%     'coOcurrencePrc', 0.15,... %require 15% chans > AmpLimits
%     'normalEEGAmpLimits', [-150, 150]); %in muV


% %% STEPSET 3 - Artefact correction
% i = i+1;  %stepSet 3
% stepSet(i).funH = { @CTAP_detect_bad_comps,... %FASTER for non-blinks
%                     @CTAP_reject_data,...
%                     @CTAP_peek_data,... %COMPARISON PEEK POINT!
%                     @CTAP_epoch_data,...%bad epochs by 100uV threshold @vertex
%                     @CTAP_detect_bad_epochs,...
%                     @CTAP_reject_data };
% %MAYBEDO - ADD @CTAP_filter_blink_ica AFTER BLINK BAD ICs?
% stepSet(i).id = [num2str(i) '_denoise_epoch'];
% 
% out.detect_bad_comps = struct(...
%     'method', 'faster',...
%     'match_measures', {{'m' 's' 'k' 'h', 'e'}},...
%     'bounds', [-2.5 2.5],...
%     'match_logic', @any);
% 
% out.epoch_data = struct(...
%     'method', 'epoch',...
%     'match',  'exact',...
%     'timelim', [-250 500],...
%     'evtype', epoch_evt);
% 
% out.detect_bad_epochs = struct(...
%     'channels', {{'A1'}},...
%     'method', 'eegthresh',...
%     'uV_thresh', [-100 100]);


% %% STEPSET 4 - Export
% %set of event names for epoching and export for event-related analysis
% erpevts = struct('hitmiss', 'hit'...
%                 , 'type', 'target'...
%                 , 'congruency', {'Congruent' 'Congruent' 'InCon' 'InCon'}...
%                 , 'shape', {'shape' 'nonShape' 'shape' 'nonShape'});
% 
% i = i + 1; %next stepSet
% stepSet(i).funH = repmat({@CTAP_export_data}, 1, numel(erpevts));
% stepSet(i).id = [num2str(i) '_export'];
% stepSet(i).save = false;
% 
% out.export_data = struct(...
%     'type', 'hdf5',...
%     'outdir', fullfile('exportRoot', 'HDF5_EXPORT'),...
%     'lock_event', {erpevts(1) erpevts(2) erpevts(3) erpevts(4)});


%% Store to Cfg
Cfg.pipe.stepSets = stepSet; % return all step sets inside Cfg struct
% step sets to run, default: whole thing
if exist('runsets', 'var')
    Cfg.pipe.runSets = runsets;
else
    Cfg.pipe.runSets = {stepSet(:).id};
end

end
