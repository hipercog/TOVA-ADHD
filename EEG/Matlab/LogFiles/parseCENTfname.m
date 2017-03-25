function [file_name,date,subject,protocol,extension]=parseCENTfname(filename)
% PARSECENTFNAME parses & returns elements from a CENT experiment filename
%
% Syntax: 
%   [file_name,date,subject,protocol,extension]=parseCENTfname(filename)
%

    % Return lower case filename for easier string comparisons
    file_name = lower(filename);
    % Get the filename without extension
    [~, subject, extension] = fileparts( file_name );
    %% Find a date: a date is at least 6 contiguous numerals
    [date, idx] = regexp( subject, '[0-9]{6,8}', 'match', 'end', 'once' );
    % If no date exists, just set the whole thing = name
    if isempty( idx ),        idx = 0;    end;
    %% The name is everything after the date
    subject = subject( idx+1:length(subject) );
    %% Find separators in the name - tokens are everything between seps
    seps = [regexpi(subject, '[_\-]') NaN];
    seps = unique([1 seps(diff(seps)~=1) length(subject)]);
    sepnum = numel(seps);
    if sepnum>2
        % reserve some space
        toks(sepnum-1).str = [];
        % Loop and extract tokens
        for i = 1:sepnum-1
            toks(i).str = subject( seps(i):seps(i+1) );
            % Find and remove separators in the token
            toks(i).str( regexpi(toks(i).str, '[_\-]') ) = [];
        end
        %% Name is everything up the protocol point
        subject = toks(1).str;
        % The protocol is the next token
        protocol = toks(2).str;
        protocol = regexprep(protocol, '[0-9]', '');
    else
        protocol='';
    end
end