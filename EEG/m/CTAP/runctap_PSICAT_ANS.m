function runctap_PSICAT_ANS(varargin)
%% Linear CTAP script to export ANS signals from CENT data
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
%   runctap_PSICAT_ANS('proj_root', <some_path>, 'group', 'ADHD')
%
% Varargin
%   user        string, an identifier string to append to output base dir
%               default, ''
%   proj_root   string, valid path to the data's project base folder, so that:
%                       <proj_root>/project_PSICAT/PSICAT-data/<group>
%               Default: '/wrk/group/hipercog/' (dir on ukko2 server)
%   group       string, name of the group to analyse, ADHD or CTRL
%               Default: ADHD
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


%% Parse input parameters
p = inputParser;

p.addParameter('user', '', @ischar)
p.addParameter('proj_root', '/wrk/group/hipercog/', @ischar)
p.addParameter('group', 'ADHD', @(x) ismember(x, {'ADHD' 'CTRL'}))
p.addParameter('sbj_filt', 'all', @(x) strcmp(x, 'all') || isnumeric(x))
p.addParameter('overwrite', true, @islogical)
p.addParameter('debug', false, @islogical)

p.parse(varargin{:});
Arg = p.Results;


%% Create the CONFIGURATION struct
% First, define step sets & their parameters in sbf_cfg()
[Cfg, ctap_args] = sbf_cfg(Arg.proj_root...
                         , ctapID...
                         , Arg.user...
                         , Arg.group);

% Next, create measurement config (MC) based on folder, & select subject subset
Cfg = get_meas_cfg_MC(Cfg, Cfg.env.paths.dataRoot...
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

end

%% %%%%%%%%%%%%%%%%%%% CONFIGURE ANALYSIS PIPE!! %%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration Subfunction
function [Cfg, out] = sbf_cfg(proj_dir, ID, user, grp)


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


%% STEPSET 1 - Load and prepare
i = 1;
stepSet(i).funH = { @CTAP_load_data,...
                    @CTAP_load_events,...
                    @CTAP_select_evdata,...
                    @CTAP_load_chanlocs,...
                    @CTAP_extract_signal,...%EXPORT ECG SIGNAL!
                    @CTAP_resample_data,...
                    @CTAP_extract_signal }; %EXPORT EDA SIGNAL!
stepSet(i).id = [num2str(i) '_load'];
stepSet(i).save = false;

out.load_events = struct(...
    'method', 'handle',...
    'handle', @loadCENTevents,...
    'src', Cfg.env.paths.dataRoot,...
    'log_ext', 'log');

out.select_evdata = struct(...
    'covertype', 'total',...
    'duration', [-2000 2000]);

out.load_chanlocs = struct(...
    'overwrite', true,...
    'index_match', false);
out.load_chanlocs.field = {{1:128 'type' 'EEG'}...
                           {{'EXG1' 'EXG2'} 'type' 'ECG'}...
                           {{'EXG3' 'EXG4' 'EXG5' 'EXG6'} 'type' 'EOG'}...
                           {{'EXG7' 'EXG8'} 'type' 'REF'}...
                           {'GSR1' 'type' 'EDA'}...
                {{'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp'} 'type' 'NA'}...
                {'EXG1' 'labels' 'ECG1'} {'EXG2' 'labels' 'ECG2'}...
                {'EXG3' 'labels' 'HEOG1'} {'EXG4' 'labels' 'HEOG2'}...
                {'EXG5' 'labels' 'VEOG1'} {'EXG6' 'labels' 'VEOG2'}...
                {'EXG7' 'labels' 'L_MASTOID'} {'EXG8' 'labels' 'R_MASTOID'}};
out.load_chanlocs.tidy  = {{'type' 'NA'}};

out.extract_signal = struct(...
    'types', {'ECG' 'EDA'},...
    'evflds', {{'code' 'congruency' 'shape' 'primertarget' 'vertices' 'hitmiss'}});
%     type
%     latency
%     urevent
%     code
%     congruency
%     shape
%     primertarget
%     vertices
%     stimtype
%     hitmiss
%     duration
%     'types', 'EDA');

out.resample_data = struct(...
    'newsrate', 32);



%% Store to Cfg
Cfg.pipe.stepSets = stepSet; % return all step sets inside Cfg struct
% step sets to run, default: whole thing
Cfg.pipe.runSets = {stepSet(:).id};

end
