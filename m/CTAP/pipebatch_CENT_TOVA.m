%% CTAP pipebatch file
update_matlab_path_anyone

%change this to switch between measurements. options: Healthy, Intake, Outtake
% SESSION = 'Healthy';
%use this to name the folder under analysis root for output
batch_id = fullfile('CENT_TOVA', SESSION);

%% Load configurations
[Cfg, my_args] = cfg_cent_TOVA(batch_id); %define paths and arguments for pipe

Cfg.MC = read_measinfo_fromfiles(...
    fullfile(Cfg.env.paths.dataRoot, 'Data_TOVA', SESSION), ...
    ['*.' my_args.load_data.type], ...
    {'[0-9]{4}[PC][_TW]', 0, 6}, ...
    {'[0-9]{4}[PC]', 0, 4}, ...
    'session', {'[0-9]{6}', 0, 4}, ...
    'measurement', {'TOVA'});


%% Select measurements to process
% Select measurements to run
clear('Filt');
Filt.subjectnr = SUBJS;
% Filt.session = '?';
% Filt.measurement = '?';
Cfg.pipe.runMeasurements = get_measurement_id(Cfg.MC, Filt);


%% Define pipeline
clear('stepSet');
i = 1; %stepSet 1
stepSet(i).funH = { @CTAP_load_data,...
                    @CTAP_load_chanlocs,...
                    @parseTOVAevents,...
                    @CTAP_reref_data};
stepSet(i).id = [num2str(i) '_LOAD_TOVA'];

i = i+1;  %stepSet 2
stepSet(i).funH = { @CTAP_extract_signal};
stepSet(i).id = [num2str(i) '_EXTRACT'];
stepSet(i).save = false;
% extract data is only to be called once -> hence a cut here

i = i+1;  %stepSet 3
stepSet(i).funH = { @CTAP_filter_data};
stepSet(i).id = [num2str(i) '_FILTER'];
stepSet(i).srcID = '1_LOAD_TOVA';
% filtering takes ages -> hence a cut here

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

i = i+1;  %stepSet 5
stepSet(i).funH = { @CTAP_blink2event,...
                    @CTAP_detect_bad_comps,... %blinks
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_comps,... %ADJUST
                    @CTAP_reject_data,...
                    @CTAP_detect_bad_comps,... %faster
                    @CTAP_reject_data};
stepSet(i).id = [num2str(i) '_IC_SEG_CORRECTION'];

i = i+1;  %stepSet 6
stepSet(i).funH = { @CTAP_interp_chan,...
                    @CTAP_peek_data};
stepSet(i).id = [num2str(i) '_INTERP'];

Cfg.pipe.stepSets = stepSet;


%% Select sets to process
Cfg.pipe.runSets = {'all'};
% Cfg.pipe.runSets = {stepSet(4:6).id};


%% Assign arguments to the selected functions, perform various checks
Cfg = ctap_auto_config(Cfg, my_args);


%% Run the pipe
%%{
tic;
% CTAP_pipeline_looper(Cfg, 'debug', true)
CTAP_pipeline_looper(Cfg)
toc;
%}

%% Cleanup
clear i stepSet Filt my_args
