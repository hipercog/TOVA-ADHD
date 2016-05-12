function anls=analysebandpwrs( fswitch, pathn )
%% FUNCTION analysebandpwrs
% PURPOSE:  Analyse something about band powers from P3-P4 C3-C4 F3-F4 ROIs for preprocessed data.
% PARAMS:   fswitch - define which analysis steps you want
% USAGE:    
% Note:     
% CALLS:    

    ops = 2;
    if nargin<1, fswitch = 0;   end;
    if nargin<2, pathn = pwd;   end;
    
    ptnt=[];
    load( fullfile( pathn, 'BandPowers.mat' ) );
    
    if any( mod( fswitch,ops )==0 )
        asym.all = zeros(3,5);
        asym.allec = zeros(3,5);
        asym.alleo = zeros(3,5);
        asym.gbyec = zeros(3,5);
        asym.gbyeo = zeros(3,5);
        asym.helec = zeros(3,5);
        asym.heleo = zeros(3,5);
        asym.vigec = zeros(3,5);
        asym.vigeo = zeros(3,5);

        % Get the aggregate band powers in each ROI
    %     aggall = ptnt(1,1);
    %     aggec = ptnt(1,1); % Start with an eyes closed
    %     aggeo = ptnt(1,2); % Start with an eyes open
    %     gbyec = ptnt(1,1); % start with first goodbye EC
    %     gbyeo = ptnt(1,2); % start with first goodbye EC
    %     helec = ptnt(1,3); % start with first hello EC
    %     heleo = ptnt(1,4); % start with first hello EO
    %     vigec = ptnt(1,5); % start with first vigi EC
    %     vigeo = ptnt(1,6); % start with first vigi EO
        aggall = zeros(6,5);
        aggec = zeros(6,5); % Start with an eyes closed
        aggeo = zeros(6,5); % Start with an eyes open
        gbyec = zeros(6,5); % start with first goodbye EC
        gbyeo = zeros(6,5); % start with first goodbye EC
        helec = zeros(6,5); % start with first hello EC
        heleo = zeros(6,5); % start with first hello EO
        vigec = zeros(6,5); % start with first vigi EC
        vigeo = zeros(6,5); % start with first vigi EO
        % Collect band powers
        for idx = 1:size( ptnt, 2 )
            for roidx = 1:6
                for bdx = 1:5
                    aggall(roidx, bdx) = aggall(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    if sum( strfind(ptnt(1,idx).prot, '_ec') ) > 0
                        aggec(roidx, bdx) = aggec(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif sum( strfind(ptnt(1,idx).prot, '_eo') ) > 0
                        aggeo(roidx, bdx) = aggeo(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    end
                    if strcmp(ptnt(1,idx).prot, 'goodbye_ec')
                        gbyec(roidx, bdx) = gbyec(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif strcmp(ptnt(1,idx).prot, 'goodbye_eo')
                        gbyeo(roidx, bdx) = gbyeo(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif strcmp(ptnt(1,idx).prot, 'hello_ec')
                        helec(roidx, bdx) = helec(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif strcmp(ptnt(1,idx).prot, 'hello_eo')
                        heleo(roidx, bdx) = heleo(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif strcmp(ptnt(1,idx).prot, 'vigilance_ec') || strcmp(ptnt(1,idx).prot, 'vigilanssi_ec')
                        vigec(roidx, bdx) = vigec(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    elseif strcmp(ptnt(1,idx).prot, 'vigilance_eo') || strcmp(ptnt(1,idx).prot, 'vigilanssi_eo')
                        vigeo(roidx, bdx) = vigeo(roidx, bdx) + ptnt(1,idx).bpws(roidx, bdx);
                    end
                end
            end
        end
        % Get asymmetries
        for bdx = 1:5
            asym.all(1,bdx) = log(aggall(1,bdx)/aggall(4,bdx));
            asym.all(2,bdx) = log(aggall(2,bdx)/aggall(5,bdx));
            asym.all(3,bdx) = log(aggall(3,bdx)/aggall(6,bdx));

            asym.allec(1,bdx) = log(aggec(1,bdx)/aggec(4,bdx));
            asym.allec(2,bdx) = log(aggec(2,bdx)/aggec(5,bdx));
            asym.allec(3,bdx) = log(aggec(3,bdx)/aggec(6,bdx));

            asym.alleo(1,bdx) = log(aggeo(1,bdx)/aggeo(4,bdx));
            asym.alleo(2,bdx) = log(aggeo(2,bdx)/aggeo(5,bdx));
            asym.alleo(3,bdx) = log(aggeo(3,bdx)/aggeo(6,bdx));

            if any(any(gbyec))
                asym.gbyec(1,bdx) = log(gbyec(1,bdx)/gbyec(4,bdx));
                asym.gbyec(2,bdx) = log(gbyec(2,bdx)/gbyec(5,bdx));
                asym.gbyec(3,bdx) = log(gbyec(3,bdx)/gbyec(6,bdx));
            end
            if any(any(gbyeo))
                asym.gbyeo(1,bdx) = log(gbyeo(1,bdx)/gbyeo(4,bdx));
                asym.gbyeo(2,bdx) = log(gbyeo(2,bdx)/gbyeo(5,bdx));
                asym.gbyeo(3,bdx) = log(gbyeo(3,bdx)/gbyeo(6,bdx));
            end
            if any(any(helec))
                asym.helec(1,bdx) = log(helec(1,bdx)/helec(4,bdx));
                asym.helec(2,bdx) = log(helec(2,bdx)/helec(5,bdx));
                asym.helec(3,bdx) = log(helec(3,bdx)/helec(6,bdx));
            end
            if any(any(heleo))
                asym.heleo(1,bdx) = log(heleo(1,bdx)/heleo(4,bdx));
                asym.heleo(2,bdx) = log(heleo(2,bdx)/heleo(5,bdx));
                asym.heleo(3,bdx) = log(heleo(3,bdx)/heleo(6,bdx));
            end
            if any(any(vigec))
                asym.vigec(1,bdx) = log(vigec(1,bdx)/vigec(4,bdx));
                asym.vigec(2,bdx) = log(vigec(2,bdx)/vigec(5,bdx));
                asym.vigec(3,bdx) = log(vigec(3,bdx)/vigec(6,bdx));
            end
            if any(any(vigeo))
                asym.vigeo(1,bdx) = log(vigeo(1,bdx)/vigeo(4,bdx));
                asym.vigeo(2,bdx) = log(vigeo(2,bdx)/vigeo(5,bdx));
                asym.vigeo(3,bdx) = log(vigeo(3,bdx)/vigeo(6,bdx));
            end
        end

        if ~any(any(asym.gbyec)),	asym.gbyec = [];    end
        if ~any(any(asym.gbyeo)),   asym.gbyeo = [];    end
        if ~any(any(asym.helec)),	asym.helec = [];    end
        if ~any(any(asym.heleo)),	asym.heleo = [];    end
        if ~any(any(asym.vigec)),	asym.vigec = [];    end
        if ~any(any(asym.vigeo)),	asym.vigeo = [];    end

        anls = asym;
    end
    
    if any( mod( fswitch,ops )==1 )
        allcell = struct2cell(ptnt);
        prots = allcell(3,1,:);
        z = size(allcell,3);
        
        % get individual protocols or parts
        allec=allcell(:,:,1:2:z);
        alleo=allcell(:,:,2:2:z);
        if any( structHasStrings(prots, 'goodbye_ec') )
            gbyec=allcell(:,:,1:6:z);
        end
        if any( structHasStrings(prots, 'goodbye_eo') )
            gbyeo=allcell(:,:,2:6:z);
        end
        if any( structHasStrings(prots, 'hello_ec') )
            helec=allcell(:,:,3:6:z);
        end
        if any( structHasStrings(prots, 'hello_eo') )
            heleo=allcell(:,:,4:6:z);
        end
        if any( structHasStrings(prots, 'vigilance_ec') ) || any( structHasStrings(prots, 'vigilanssi_ec'))
            vigec=allcell(:,:,5:6:z);
        end
        if any( structHasStrings(prots, 'vigilance_ec') ) || any( structHasStrings(prots, 'vigilanssi_ec'))
            vigeo=allcell(:,:,6:6:z);
        end
        
        allbp = cell2mat(allcell(1,1,:));
        allec = cell2mat(allec(1,1,:));
        alleo = cell2mat(alleo(1,1,:));
        gbyec = cell2mat(gbyec(1,1,:));
        gbyeo = cell2mat(gbyeo(1,1,:));
        helec = cell2mat(helec(1,1,:));
        heleo = cell2mat(heleo(1,1,:));
        vigec = cell2mat(vigec(1,1,:));
        vigeo = cell2mat(vigeo(1,1,:));
        
        % Get asymmetries
        asym.all = bands2asyms(allbp);
        asym.allec = bands2asyms(allec);
        asym.alleo = bands2asyms(alleo);
        asym.gbyec = bands2asyms(gbyec);
        asym.gbyeo = bands2asyms(gbyeo);
        asym.helec = bands2asyms(helec);
        asym.heleo = bands2asyms(heleo);
        asym.vigec = bands2asyms(vigec);
        asym.vigeo = bands2asyms(vigeo);
        
        % Return
        anls = asym;
    end
end

function asyms=bands2asyms( bands )

    r = size(bands,1)/2;
    c = size(bands,2);
    z = size(bands,3);
    asyms = zeros(r, c, z);
    for rdx = 1:r
        for cdx = 1:c
            asyms(rdx,cdx,:) = log(bands(rdx,cdx,:)./bands(rdx+r,cdx,:));            
        end
    end
end

function logarr=structHasStrings( stc, str )

    logarr=~cellfun(@isempty, strfind( stc, str ));
    logarr=logarr(:)';
end