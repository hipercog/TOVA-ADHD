function batchproc(indir,conv,outdir,trodes)
% BATCHPROC 

if ~ischar(indir)
    alleeg=indir;
    if nargin<3,    outdir=pwd; end
else
    files=dir(indir);
    [d,~,~]=fileparts(indir);
    if nargin<3,    outdir=indir;
    elseif isempty(strfind(outdir,filesep)), outdir=fullfile(d, outdir);
    end
    if ~isempty(outdir) && ~isdir(outdir),	mkdir(outdir);    end
    alleeg=[];
end
if nargin<4,    trodes=1:128;   end
ntro=numel(trodes);

% COND={'street','letter','ykanji'};
COND={'b1','b2','b3','b4'};
conds={'P','C_'};
% conds={'street','letter'};

% choice of post-processing functions
switch conv
    case 'fixname'
        for i=1:numel(files)
            eeg=pop_loadset('filename',files(i).name,'filepath',d);
            eeg.setname=[eeg.CTAP.subject '_' eeg.CTAP.protocol];
            ctapeeg_save(eeg,'outdir',outdir);
        end
        return;
    case 'extract'
        for i=1:numel(files)
            eeg=pop_loadset('filename',files(i).name,'filepath',d);
            ok=exportEpochsByCond(eeg,outdir,COND);
            disp(['Condition extraction for ' files(i).name ' ' ok]);
        end
        return;
    case 'ersp'
        if isempty(alleeg),    alleeg=getmyeeg(indir); end
        sig=cell(ntro,1);
        for i=1:ntro
            if isempty(conds)
                sig{i}=c_gstlt_ersp(alleeg,trodes(i),outdir);
            else
                sig{i}=c_gstlt_ersp(alleeg,trodes(i),outdir,conds);
            end
        end
        save(fullfile(outdir,'sig.mat'),'sig');
    case 'erp'
        if isempty(alleeg),    alleeg=getmyeeg(indir); end
        sig=cell(ntro,1);
        for i=1:ntro
            if isempty(conds)
                sig{i}=c_gstlt_erp(alleeg,trodes(i),outdir);
            else
                sig{i}=c_gstlt_erp(alleeg,trodes(i),outdir,conds);
            end
        end
        save(fullfile(outdir,'sig.mat'),'sig');
end
end