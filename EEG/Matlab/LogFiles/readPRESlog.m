function data=readPRESlog( filename, filter, conv, hdr, savit )
% READPRESLOG takes a filename or directory and parses the given log file(s)
if nargin<5, savit=false;   end
if nargin<4,    hdr=5;      end
if nargin<3, conv='CENT';   end
if nargin<2, filter='';     end
data=struct;
[directory, fname, ext] = fileparts(filename);
s = filesep;
looper = 1;
% If filename is a directory, not a file
if isdir(filename)
    if ~isempty(ext) % if extension is not empty when filename is a directory, path contained period
        fname = strcat( fname, ext );
    end        
    % Make sure we have a well formed path with separators
    directory = fullfile(directory, fname, s);
    % Obtain all at once
    files = dir(fullfile(directory, [filter '*.log']));
    % Get complete path for files?
%         [directory,~,~]=fileparts(which(files(1).name));
    looper = numel(files);
else
    % This is a single file case - check file exists, and check file format
    if exist( filename, 'file' ) && strcmpi( ext, '.log' )==1
        files=dir(filename);
    else
        disp('Bad file or filename - check and retry');
        return;
    end
end
for i=1:looper
    filename=files(i).name;
    % Parse the log file name according to which protocol it came from
    switch conv
        case 'CENT'
            [~, ~, subj, prot, ~] = parseCENTfname( filename );
        case 'TTL'
            [~, prot, subj, ~, ~] = parseTTLfname( filename );
    end
    temp=readonelog(fullfile(directory, filename),hdr);
    data.([prot subj])=temp;
end
% cleanup and save
if savit,	save(fullfile(directory,'logs'), '-struct', 'data'); end
end

function data=readonelog( fname,hdr )
% READONELONG parses a single log file: extracts trials of interest
% and compiles various stats calling getstats(), to build 'data'.

    % Read file
    if exist(fname, 'file')
        fid = fopen(fname);
        textscan(fid,'%s',hdr,'Delimiter','\n');
        textcols=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s','Delimiter','\t');
    else
        error( [fname ' does not exist. Abort'] );
    end
    trial =cellfun(@str2double, textcols{2});
    idx=find(isnan(trial),1,'first');
    trial(idx:end)=[];
    data.type = textcols{3}(1:idx-1);
    data.code = textcols{4}(1:idx-1);
    data.time = cellfun(@str2double, textcols{5}(1:idx-1));
    data.dur = cellfun(@str2double, textcols{8}(1:idx-1));
    data.trial=trial;
end
