function batch( filetype, indir, outdir )

if nargin < 1,	filetype = 'tif';  end     % image type
if nargin < 2,	indir = pwd;  end
if nargin < 3,	outdir = fullfile(pwd, 'plots');   end

% ------------------------------------------------------------------------------------------------------------------------------------------
% Get the list of the input image names.
files = dir(fullfile(indir, ['*.' filetype]));
if isempty(files)
    disp('No files of the requested image type - abort');
    return;
end
nfiles = numel(files);
% File and folder stuff
if strcmpi( indir, outdir ) == 1,     outdir = fullfile(outdir, 'plots');    end
if ~isdir(outdir),  mkdir(outdir);  end
% Open the output txt-file.
fullOutputName = fullfile(outdir, 'some textual info?.txt');
fid = fopen(fullOutputName,'a');

% ------------------------------------------------------------------------------------------------------------------------------------------
% LOOP for the names in the list
for i=1:nfiles
    filename = files(i).name;
    % do something
end