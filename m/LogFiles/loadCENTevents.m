function [EEG, varargout] = loadCENTevents(EEG, varargin)
% LOADCENTEVENTS loads events from CENT style Pres log files to eeglab data
%
%   DESCRIPTION:
%   Convert events from raw EEG triggers to log strings. 
%   Obtain an offset latency array by matching some event names to codes, 
%   then use interpolation to smoothly vary the offset.
%
% Syntax:
%   EEG = loadCENTevents(EEG, varargin)
%
% Inputs:
%   EEG             struct, EEG struct to process
%
% Varargin
%   'src'           string, directory containing, or full path to, log file
%                           If empty string is passed, 'log_src' will be set
%                           equal to directory of original EEG file, obtained
%                           from 'EEG.comments' field
%   'mkplot'        logical, plot interpolation of event times
%                   Default: false
%   'name'          string, 
%   'protocol'      string, 
%   'date'          string, 
%
% Output:
%   'EEG'       : EEG struct, with updated events structure
%
% NOTE:     
%
% CALLS:    eeglab functions
%
% Version History:
% 20.10.2014 Created (Benjamin Cowley, FIOH)
%
% Copyright 2014- Benjamin Cowley, benjamin.cowley@ttl.fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: varargout is dummy to support interface with CTAP pipeline
varargout{1} = '';

%% Parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('EEG', @isstruct)

p.addParameter('src', '', @ischar)
p.addParameter('mkplot', false, @islogical)

p.parse(EEG, varargin{:});
Arg = p.Results;

log_src = Arg.src;


%% Get the log file source and naming data for the EEG
if ~isempty(log_src)
    [origEEGdir, fname, e] = fileparts(log_src);
elseif ~isempty(EEG.filename)
    [origEEGdir, fname, e] = fileparts(EEG.filename);
elseif ~isempty(EEG.comments)
    [origEEGdir, fname, e] = fileparts(strrep(EEG.comments, 'Original file: ', ''));
end
if ~isfolder(origEEGdir)
    error('loadCENTevents:no_source',...
        [''  ' - auto-found directory does not exist.']);
end
log_src = fullfile(origEEGdir, [fname '.' strrep(e, '.', '')]);
[~, Arg.date, Arg.name, Arg.protocol, ~] = parseCENTfname( fname );

% If directory is not given, look in the home of the original EEG file
    
% Select the right values for the protocol we see
npt = pickfmatch({'salience','vigilance','gestalt','hello','goodbye','tova'}...
                , Arg.protocol);
% Find if this is test group saliency data or not
if npt ~= 1
    salt = false;
%TODO - Below breaks subject 2022C in Saliency - check is needed for gestalt?
% elseif str2double(Arg.date(1:4)) == 2013 
%     npt = 3;
%     salt = false;
else
    salt = true;
end


%% changes events into strings, after removing two superlative bytes
if all(cellfun(@ischar, {EEG.event.type}))
    evt = cellfun(@str2double, {EEG.event.type});
elseif all(cellfun(@isnumeric, {EEG.event.type}))
    evt = [EEG.event.type];
else
    error('loadCENTevents:bad_EEG_evts', 'EEG.event.type field has mixed data!')
end
if max(evt) > 32767, byte_off = 65280; else byte_off = 32512; end %#ok<SEPEX>
for k = 1:length(evt)
    EEG.event(k).type = num2str(evt(k) - byte_off);
end
idx = ismember({EEG.event.type}, '0');
EEG.event(idx) = [];


%% Grab events from Pres log for any protocol except TOVA

% Get the events from a given file or directory
if isfolder(log_src)
    logs = dir( fullfile(log_src, '*.log') );
    if isempty(logs)
        disp(['No log files found at given directory: ' log_src...
            ' Aborting Presentation event copy-write.']);
        return; 
    end
    for i=1:numel(logs)
        [~, d, n, p, ~] = parseCENTfname( logs(i).name );
        logs(i).parsedname=[d '-' n '-' p];
    end
    m = pickfmatch(....
        {logs.parsedname}, [Arg.date '-' Arg.name '-' Arg.protocol] );
    log_src = fullfile(log_src, logs(m).name);
end
% Check the names match
[~, ~, logname, ~, ~] = parseCENTfname( fname );
if strcmp( logname, Arg.name ) == 0
    disp(['Log = ' logname '; EEG = ' Arg.name '; PRES LOG WRITE FAIL']);
    return;
end
% Import the events
[~, log_data] = importPresLog( log_src, salt );
% Pres Log event time is in secs/10000, convert to latency=samples/sec
log_data.latencies = log_data.time * EEG.srate / 10000;
% Given the differences introduced by changing Pres logs, we do
if ~salt,	log_data.coords_str = log_data.code;	end
nl = numel(log_data.code);


%% Bring EEG events and Log events into alignment -
% this depends on EEG events being a subset of Pres log

offset = NaN( numel(log_data.code), 1 );
% Find log-event to trigger matches depending on protocol
lgm=[]; tgm=[];
switch npt
    case 1
        % For Saliency: match on log code 'Fixation Cross', 
        % trigger 143 or union of 143 and 142
        lgm = find( strcmp( 'Fixation Cross', log_data.coords_str ) );
        tgm = find( strcmp( '143', {EEG.event.type} ) );
        if length(lgm) ~= length(tgm)
            tgm = union( find(strcmp('142', {EEG.event.type})), tgm );
        end
    case 2
        % For Vigilance: match on log code 'Mental Arithmetic Result',
        % 'Eyes open' 'Resting end', trigger 55, 7, 100
        lgm = find(...
            strcmp('Mental Arithmetic Result', log_data.coords_str)+...
            strcmp('Eyes open', log_data.coords_str ) +...
            strcmp('Resting End', log_data.coords_str ) );
        tgm = find(strcmp( '55', {EEG.event.type} ) +...
                    strcmp( '7', {EEG.event.type} ) +...
                  strcmp( '100', {EEG.event.type} ) );
    case 3
        nty=str2double({EEG.event.type})';
        %         strport = 201:210;
        se=find(~cellfun(@isempty,strfind(log_data.coords_str,'_street')));
        sp=find(nty>200 & nty<215);
        if length(se)==length(sp),  lgm=se;  tgm=sp;  end
        %         ltrport = 217:226;
        le=find(~cellfun(@isempty,strfind(log_data.coords_str,'_letter')));
        lp=find(nty>216 & nty<230);
        if length(le)==length(lp) && isempty(intersect(le,lgm)) &&...
                isempty(intersect(lp,tgm))
            lgm=union(lgm,le);  tgm=union(tgm,lp);
        end
        %         ykiport = 233:242;
        ye=find(~cellfun(@isempty,strfind(log_data.coords_str,'_ykanji')));
        yp=find((nty>232 & nty<244));
        if length(ye)==length(yp) && isempty(intersect(ye,lgm)) &&...
                isempty(intersect(yp,tgm))
            lgm=union(lgm,ye);  tgm=union(tgm,yp);
        end
    otherwise
        % For Hello/Goodbye: match on log code 'End...', triggers 141
        lgm = find(...
            strcmp('End resting eyes open', log_data.coords_str ) +...
            strcmp('End resting eyes closed', log_data.coords_str )+...
            strcmp('End passive white noise', log_data.coords_str )+...
            strcmp('End active white noise', log_data.coords_str ) );
        tgm = find(...
            strcmp( '110', {EEG.event.type} ) +...
            strcmp( '120', {EEG.event.type} ) + ...
            strcmp( '130', {EEG.event.type} ) +...
            strcmp( '140', {EEG.event.type} ) );
end
% Give warning if #events ~= #triggers, offset is just 1 value
if isempty(lgm) || length(lgm) ~= length(tgm)
    disp( '*****WARNING: LOG EVENTS AND TRIGGERS DID NOT MATCH*****' );
    clear offset;
    switch npt
        case 1
            % For Saliency: match on log codes 'End task', triggers 255
            l1 = find( strcmp( 'End task', log_data.coords_str ) );
            t1 = find( strcmp( '255', {EEG.event.type} ) );
        case 2
            % Vigilance: match on log code 'Mental Arithmetic Result',
            % trigger 55
            l1=find(...
               strcmp('Mental Arithmetic Result',log_data.coords_str));
            t1 = find( strcmp( '55', {EEG.event.type} ) );
        case 3
            l1 = find( strcmp( 'Start task', log_data.coords_str ) );
            t1 = find( strcmp( '254', {EEG.event.type} ) );
        otherwise
            % For Hello/Goodbye: match on log code 
            % 'End resting eyes closed', trigger 120
            l1=find(...
               strcmp('End resting eyes closed',log_data.coords_str));
            t1 = find( strcmp( '120', {EEG.event.type} ) );
    end
    if length(t1) == 1 && length(l1) == 1
        offset = EEG.event(t1).latency - log_data.latencies(l1);
    else
        error('loadCENTevents:no_match', 'PRES-LOG EVENT IMPORT FAILS FOR: %s'...
            , [Arg.name '-' Arg.protocol] );
    end
else
    % Calculate the per-samples offset between equivalent events
    for lt = 1:numel(lgm)
        offset(lgm(lt)) =...
            EEG.event(tgm(lt)).latency - log_data.latencies(lgm(lt)); 
    end
    % Interpolate the sparse array of offsets & subtract from log times
    offset = inpaint_nans(offset);
end
% Add offsets to put Pres Log times in EEG time frame
log_data.latencies = log_data.latencies + offset;
if Arg.mkplot && length(lgm) == length(tgm)
    evsPosTest(nl, lgm, tgm, EEG.event, log_data.latencies);
end


%% Create event structures - with extra columns for salience test.

if salt
    temp=~cellfun(@isempty,strfind(log_data.pritar_str,'target')) |...
        ~cellfun(@isempty,strfind(log_data.pritar_str,'primer'));
    for f=1:nl
        if temp(f), log_data.code{f} = log_data.pritar_str{f};  end
    end
    temp = struct(...
    'type',log_data.code,...
    'latency',num2cell(log_data.latencies),...
    'urevent',num2cell([1:nl]'),...
    'code',num2cell(log_data.cond_num),...
    'congruency',log_data.congruency_str,...
    'shape',log_data.shape_str,...
    'primertarget',log_data.pritar_str,...
    'vertices',num2cell(log_data.vertices_num),...
    'stimtype',log_data.event_type,...
    'hitmiss',log_data.stim_type); %#ok<NBRAK>
else
    temp = struct('type', log_data.code,...
        'latency', num2cell(log_data.latencies),...
        'urevent', num2cell([1:nl]'),...
        'evntype', log_data.event_type,...
        'hitmiss', log_data.stim_type); %#ok<NBRAK>
end
clear log_data;
% Add Pres log events to EEG event list?
if offset == 0
    % Here we trust that offset cannot spontaneously = 0
    EEG.no_off_ev = temp;
else
    % If all went well, we have new log events + old triggers!
    EEG.oldev = EEG.event;
    EEG.event = temp;
    EEG.urevent = EEG.event;
end
end

function evsPosTest(nl, lgm, tgm, event, latencies)
%% EVSPOSTEST plot the latencies of triggers from the raw data and events 
%   from the log file to see how closely they match.
%
%   SYNTAX  : evsPosTest(nl, lgm, tgm, event, latencies)
%       'nl': number of events
%       'lgm': 
%   	'tgm': 
%   	'event': 
%   	'latencies': 
%
    test = NaN(nl, 1);
    i = 1;
    e = 1;
    lt = 1;
    while lt <= nl && i <= numel(tgm) && e <= numel(event)

        if e ~= tgm(i)
            test(lt) = event(e).latency;
            e = e + 1;
        else
            test(lgm(i)) = event(e).latency;
            lt = lgm(i);
            e = e + 1;
            i = i + 1;
        end
        lt = lt + 1;
    end
    test = inpaint_nans(test);
    plot(test);
    hold all;
    plot(latencies)
end