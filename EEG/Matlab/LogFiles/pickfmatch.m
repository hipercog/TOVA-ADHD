function file = pickfmatch( filenames, prot )
%   filenames:  cell array of names
%   prot:       the element from filenames to match

%%
NF = numel(filenames);
f = zeros(1, NF);

for i = 1:NF
    [~, fname, ~] = fileparts( filenames{i} );
%    [boo boo boo boo protocol] = parsefname( files(i).name, informat );
    f(i) = strdist( fname, prot );
end

[~, file] = min(f);


% function [file] = pickfmatch( files, filename )
% 
% NF = numel(files);
% file = zeros(NF);
% 
% for i = 1:NF
%     file(i) = strdist( files(i).name, filename, 2, 1 );
% end
% 
% [~, file] = min(file);