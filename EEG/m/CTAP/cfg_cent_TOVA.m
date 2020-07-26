function [Cfg, args] = cfg_cent_TOVA(idstr)

if nargin < 1, idstr = 'CENT_TOVA'; end

% Cfg.grfx.on = false; %set false to skip heavy QA outputs & be fast :)

%% Cross-platform stuff (not part of CTAP)
if strcmp(computer,'GLNXA64')
    Cfg.env.paths.serverRoot = '/ukko';
elseif strcmp(computer,'MACI64')
    Cfg.env.paths.serverRoot = '/ukko';
else
    Cfg.env.paths.serverRoot = 'U:\';
end

Cfg.env.paths.projectRoot = fullfile(...
    Cfg.env.paths.serverRoot,'projects','ReKnow');


%% Define important directories and files
Cfg.env.paths.analysisRoot = fullfile(...
    Cfg.env.paths.projectRoot, 'analysis', 'ctap', idstr);

Cfg.env.paths.dataRoot = fullfile(...
    Cfg.env.paths.projectRoot, 'Data', 'CENT');

% Note: other canonical locations are added in ctap_auto_config.m
% You should use it in your analysis batch file.

Cfg.eeg.chanlocs = fullfile(...
    Cfg.env.paths.projectRoot,'Data','raw','neurone',...
                    'channel_locations_acticap_32.ced');

Cfg.eeg.chanlocs = fullfile(...
    Cfg.env.paths.dataRoot, 'chanlocs128_cent.elp');


%% Define other important stuff
Cfg.eeg.reference = {'L_MASTOID' 'R_MASTOID'};

% EOG channel specification for artifact detection purposes
% Allowed values: {},{'chname'},{'chname1','chname2'}
% In case of two channel names their abs(ch1-ch2) is used as the signal.
Cfg.eeg.veogChannelNames = {'VEOG1' 'VEOG2'};
Cfg.eeg.heogChannelNames = {'HEOG1' 'HEOG2'};


%% Configure analysis functions

% STEP 1 - LOADING AND PREPPING
args.load_data = struct(...
    'type', 'bdf');%data path is read from directory at Cfg.env.paths.dataRoot

%could also just call parseTOVAevents() directly in the pipe?
args.load_events = struct(...
    'method', 'handle',...
    'handle', @parseTOVAevents);

args.load_chanlocs.tidy = {{'type' ''}};
args.load_chanlocs.types = {...
    {'1:128' 'EEG'}, {'129:130' 'ECG'}, {'131:134' 'EOG'},...
    {'135:136' 'REF'}, {'137' 'EDA'}, {'138:143' ''}};

args.extract_signal.types = {'ECG' 'EDA'};

args.reref_data = struct('reference', 'REF');

% STEP 2 - FILTER
args.filter_data = struct(...
    'lowCutOff', 2,...
    'highCutOff', 45);

% STEP 3 - SELECT, BAD CHANNELS, ICA
args.select_evdata = struct(...
    'evtype', 'TOVAtest');

% args.peek_data = struct();

%For all 'detect_bad_x' you should specify at least a method
args.detect_bad_channels = struct(...
    'method', 'faster',...
    'refChannel', 'C21',... %should be 'C21' alias 'Fz'
    'channelType', {'EEG'});

args.detect_bad_segments = struct(...
    'amplitudeTh', [-80, 80]); %in muV [-100, 100]

args.run_ica = struct(...
    'method', 'fastica');
args.run_ica.channels = {'EEG' 'EOG'};

% STEP 4 - BAD COMPS, SEGS
args.detect_bad_comps = struct(...
    'method', {'blink_template' 'adjust' 'faster'},...
    'adjustarg', {'' 'horiz' ''});

