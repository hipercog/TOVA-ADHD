function out = precomp_wrap(indf, varargin)
% ------------------------------------------------

p = inputParser;

p.addRequired('indf', @(x) isstruct(x) || isfolder(x) || isfile(x) )
p.addParameter('ALLEEG', [], @isstruct)

p.addParameter('erp', 0, @(x) any(x == [0 1]) || islogical)
p.addParameter('rmbase', [-999 -1], @isnumeric)

p.addParameter('spec', 0, @(x) any(x == [0 1]) || islogical)

p.addParameter('erpim', 0, @(x) any(x == [0 1]) || islogical)
p.addParameter('nlines', 1000, @isnumeric)
p.addParameter('smoothing', 100, @isnumeric)
p.addParameter('sorttype', 'cor_rsp', @ischar)
p.addParameter('sortfield', 'duration', @ischar)

p.addParameter('ersp', 0, @(x) any(x == [0 1]) || islogical)
p.addParameter('itc', 0, @(x) any(x == [0 1]) || islogical)

p.parse(indf, varargin{:})
Arg = p.Results;

sw = {'off' 'on'};
out = cell(1, sum([Arg.erp Arg.spec Arg.erpim Arg.ersp Arg.itc]));


%% Get STUDY & ALLEEG ready
if isstruct(indf)
    STUDY = indf;
    if isfield(Arg, 'ALLEEG')
        ALLEEG = Arg.ALLEEG;
    else
        error('precomp_wrap:no_eeg', 'No ALLEEG data provided')
    end
elseif isfile(indf) %treat as an existing study
    [p, f, e] = fileparts(indf);
    [STUDY, ALLEEG] = pop_loadstudy('filename', [f e], 'filepath', p);
else
%TODO:  what to do with a folder?
end


%% PRECOMPUTE
out = std_precomp(STUDY, ALLEEG, {}...
    , 'savetrials', 'on'...
    , 'interp', 'on'...
    , 'recompute', 'on'...
    , 'erp', sw{Arg.erp + 1}...
    , 'erpparams', {'rmbase', Arg.rmbase}...
    , 'spec', sw{Arg.spec + 1}...
    , 'specparams', {'specmode', 'fft', 'logtrials', 'off'}...
    , 'erpim', sw{Arg.erpim + 1}...
    , 'erpimparams', {'nlines', Arg.nlines...
                    , 'smoothing', Arg.smoothing...
                    , 'sorttype', Arg.sorttype...
                    , 'sortfield', Arg.sortfield}...
    , 'ersp', sw{Arg.ersp + 1}...
    , 'itc', sw{Arg.itc + 1}...
    , 'erspparams', {'type', 'ersp&itc'...
                   , 'freqs', [3 30]...
                   , 'cycles', [3 0.8]...
                   , 'padratio', 3 ...
                   , 'alpha', 0.05});

