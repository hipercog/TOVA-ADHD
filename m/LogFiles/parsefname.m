function [file_name date name prot prep] = parsefname( filename )

    % Return lower case filename for easier string comparisons
    file_name = lower(filename);
    % Get the filename without extension
    [~, name ~] = fileparts( file_name );
    
    %% Find a date: a date is at least 6 contiguous numerals
    [date idx] = regexpi( name, '[0-9]{6,8}', 'match', 'end', 'once' );
    % If no date exists, just set the whole thing = name
    if isempty( 'idx' ),        idx = 0;    end;
    
    %% The name is everything after the date
    name = name( idx+1:length(name) );
    % Find separators in the name - tokens are everything between seps
    seps = [regexpi(name, '[_\-]') NaN];
    seps = [seps(diff(seps)~=1) length(name)];
%     seps(seps==1) = [];
    sepnum = numel(seps);
    % reserve some space
    toks(sepnum-1).str = [];
    % Loop and extract tokens
    for i = 1:sepnum-1
        toks(i).str = name( seps(i):seps(i+1) );
        % Find and remove separators in the token
        toks(i).str( regexpi(toks(i).str, '[_\-]') ) = [];
    end

    % Name is everything up the protocol point
    name = toks(1).str;
    % The protocol is the next token
    prot = toks(2).str;
    % The rest of the tokens describe preprocessing history
    toks(1:2) = [];
    prep = toks;
end