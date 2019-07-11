function [file_name,project,subject,msrmnt,extension]=parseTTLfname(filename)
% PARSETTLFNAME parses & returns elements from a TTL experimental filename
%
% Syntax:
%   [file_name,date,subject,protocol,extension]=parseTTLfname(filename)
%

    % Return lower case filename for easier string comparisons
    file_name = lower(filename);
    % Get the filename without extension
    [~, temp, extension] = fileparts( file_name );

    temp=regexprep(temp, '[_\-]', ' ');
    temp=textscan(temp,'%s');
    project=regexpi(temp{1},'[A-z]','match');
    subject=regexpi(temp{1},'[0-9]','match');
    msrmnt='';
    for i=1:numel(temp)
        msrmnt=[msrmnt '_' temp{i}]; %#ok<AGROW>
    end
end