function [EEG, varargout] = parseTOVAevents(EEG, varargin)
% PARSETOVAEVENTS parses TOVA events! Currently only works for TOVA data loaded
% with pop_biosig, as ctap_readbdf doesn't find events
% Note: varargin and varargout are dummies to support interface with CTAP
% pipeline


varargout{1} = '';
full_codes = {'TOVA test START' 'NON-TARGET' 'TARGET' 'RESPONSE' 'undefined'};
codes = {'STRT' 'NONT' 'TRGT' 'RESP' 'UNDF'};
trigs = [65528 65529 65531 65532];
metas = {'COM' 'cor_rsp' 'cor_inb' 'OM'};
evtypes = [EEG.event.type];
[C, ~, ic] = unique(evtypes);
if any(~ismember(trigs(2:4), C))
    error('parseTOVAevents:badTriggers'...
        , 'Missing either NON-TARGET, TARGET, or RESPONSE')
end


%% get unambiguous start and end times of test, i.e. before edited events
%start S = max(L1 - sr, L1 / 2), where L=START marker latency, sr=sample rate
%duration D = min(P, Ln + sr) - S, where P=EEG samples, Ln=last event latency
test_start = EEG.event(1).latency;
test_start = max(test_start - EEG.srate, round(test_start / 2));
test_dur = min(EEG.pnts, EEG.event(end).latency + EEG.srate) - test_start;


%% get meta-events from eventcode sequential pairs
ev = [];
compare = @(x, y) x == y;
sbf_find_meta_event_pairs(trigs(2), metas{1})
sbf_find_meta_event_pairs(trigs(3), metas{2})

compare = @(x, y) x ~= y;
sbf_find_meta_event_pairs(trigs(2), metas{3})
sbf_find_meta_event_pairs(trigs(3), metas{4})


%% report distribution info, change types to strings
for i = 1:numel(C)
    if i < 5, idx = i; else, idx = 5; end
    fprintf('%d==%s (new code=''%s''), occurs %d times\n'...
        , C(i), full_codes{idx}, codes{idx}, sum(ic == i))
    logidx = find(evtypes == C(i));
    [EEG.event(logidx).type] = deal(codes{idx});
end


%% create one event with latency and duration to cover test period +-1 second
ev = union_struct(ev...
    , eeglab_create_event(test_start, 'TOVAtest', 'duration', {test_dur}));


%% join all event tables
[EEG.event, ~] = eeglab_merge_event_tables(ev, EEG.event, 'ignoreDiscontinuousTime');


%% subfunction
function sbf_find_meta_event_pairs(evnum, typestr)
    logidx = find(evtypes(1:end - 1) == evnum);
    metidx = compare([EEG.event(logidx + 1).type], trigs(4));
    if any(metidx)
        hit = logidx(metidx);
        event =...
            eeglab_create_event([EEG.event(hit).latency], typestr, 'duration',...
            num2cell([EEG.event(hit + 1).latency] - [EEG.event(hit).latency]));
        ev = union_struct(ev, event);
    end
end

end %parseTOVAevents()
