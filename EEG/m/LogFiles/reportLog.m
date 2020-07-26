function [subj stats log]=reportLog( proto, path, log )
%% function reportLog() will spit out means and  from logs of the given protocol
%  Parameters
%   proto   -   name of the protocol to be analysed.
%   path    -   path to the Presentation log files
%   log     -   name of a mat file which has the already imported log data
if nargin < 3
    files = dir( fullfile( path, '*.log' ) );
    fn = numel(files);
    log = cell(1,fn);
else
    cd( path );
    if ischar(log)
        log = load(log);
        log = log.log;
    end
    fn = size(log, 2);
end
trials = 550;

%% Do the Saliency protocol
if strcmp( proto, 'salien' ) == 1
    % init behaviour data structs: subjects are in column order
    conds = struct('cnh', {}, 'cnm', {}, 'csh', {}, 'csm', {}, 'inh', {}, 'inm', {}, 'ish', {}, 'ism', {});
    subj = struct('name', {}, 'conds', {}, 'rtvar', {}, 'err', {}, 'errCon', {}, 'errInCon', {}, 'errs', {}, 'hitrate', {}, 'falarm', {}, 'dprime', {});
    condnames = {'cnh', 'cnm', 'csh', 'csm', 'inh', 'inm', 'ish', 'ism'};
    % RT stuff
    trg = cell(1,fn);
    rsp = cell(1,fn);
    tlt = zeros(trials,fn);
    rlt = zeros(trials,fn);
    rt = zeros(trials,fn);
    % Error/hit rate stuff
    err = ones(1,fn);
    errs = ones(4,fn);
    errCon = ones(1,fn);
    errInCon = ones(1,fn);
    hr = ones(1,fn);
    fa = ones(1,fn);
    hrc = ones(1,fn);
    fac = ones(1,fn);
    hri = ones(1,fn);
    fai = ones(1,fn);
    % Loop once for every subject
    for i = 1:fn
        if nargin < 3,	[~, log{1,i}] = importPresLog(fullfile(path, files(i).name), true); end;
        trg{i} = strcmp( 'target', log{1,i}.pritar_str );
        rsp{i} = circshift(trg{i},1);
        tlt(:,i) = log{1,i}.time(trg{i});
        rlt(:,i) = log{1,i}.time(rsp{i});
        rt(:,i) = rlt(:,i) - tlt(:,i);
        % Get trials per condition - ALL CONDS ARRANGED ALPHABETICALLY
            % Congruent non-shape hit condition
        conds(i).cnh=findConditions( log{1,i}, 'Congruent', 'nonShape', 'hit' );
            % Congruent non-shape miss condition
        conds(i).cnm=findConditions( log{1,i}, 'Congruent', 'nonShape', 'miss' );
            % Congruent shape hit condition
        conds(i).csh=findConditions( log{1,i}, 'Congruent', 'shape', 'hit' );
            % Congruent shape miss condition
        conds(i).csm=findConditions( log{1,i}, 'Congruent', 'shape', 'miss' );
            % INCongruent non-shape hit condition
        conds(i).inh=findConditions( log{1,i}, 'InCon', 'nonShape', 'hit' );
            % INCongruent non-shape miss condition
        conds(i).inm=findConditions( log{1,i}, 'InCon', 'nonShape', 'miss' );
            % INCongruent shape hit condition
        conds(i).ish=findConditions( log{1,i}, 'InCon', 'shape', 'hit' );
            % INCongruent shape miss condition
        conds(i).ism=findConditions( log{1,i}, 'InCon', 'shape', 'miss' );
        
        % Subject-level deets & stats
        subj(i).name = log{1,i}.subject{1}(10:end);
        subj(i).conds = conds(i);
        subj(i).rtvar = mean([conds(i).cnh.rtvar conds(i).csh.rtvar conds(i).inh.rtvar conds(i).ish.rtvar]);
        % Errors are in proportionate form: as a percentage of ALL trials, not as a ratio of misses to hits.
        subj(i).err = (conds(i).cnm.size + conds(i).csm.size + conds(i).inm.size + conds(i).ism.size)*100/trials;
        err(i) = subj(i).err;
        subj(i).errs = [(conds(i).cnm.size*100)/(conds(i).cnh.size + conds(i).cnm.size)...
                        (conds(i).csm.size*100)/(conds(i).csh.size + conds(i).csm.size)...
                        (conds(i).inm.size*100)/(conds(i).inh.size + conds(i).inm.size)...
                        (conds(i).ism.size*100)/(conds(i).ish.size + conds(i).ism.size)];
        errs(:,i) = subj(i).errs;
        subj(i).errCon = (subj(i).errs(1)+subj(i).errs(2))/2;
        errCon(i) = subj(i).errCon;
        subj(i).errInCon = (subj(i).errs(3)+subj(i).errs(4))/2;
        errInCon(i) = subj(i).errInCon;
        % d prime stuff: hit rate from shapes, false alarm from nonshapes (gives equivalent dprime)
        hr(i) = ((conds(i).csh.size + conds(i).ish.size))/(conds(i).csh.size + conds(i).csm.size + conds(i).ish.size + conds(i).ism.size);
        fa(i) = ((conds(i).cnm.size + conds(i).inm.size))/(conds(i).cnh.size + conds(i).cnm.size + conds(i).inh.size + conds(i).inm.size);
        hrc(i) = conds(i).csh.size/(conds(i).csh.size + conds(i).csm.size);
        fac(i) = conds(i).cnm.size/(conds(i).cnh.size + conds(i).cnm.size);
        hri(i) = conds(i).ish.size/(conds(i).ish.size + conds(i).ism.size);
        fai(i) = conds(i).inm.size/(conds(i).inh.size + conds(i).inm.size);
        subj(i).hitrate = hr(i);
        subj(i).falarm = fa(i);
    end
    % Extract some subject stats
    hr = hr-0.00001;    fa = fa+0.00001;
    hrc = hrc-0.00001;   fac = fac+0.00001;
    hri = hri-0.00001;   fai = fai+0.00001;
    mdprime = dprime(hr, fa);
    mdprimec = dprime(hrc, fac);
    mdprimei = dprime(hri, fai);
    for c=1:fn
        subj(c).dprime = mdprime(c);
    end
    % Extract the overall RT stats
    stats.rt=rt;
    stats.meanrt = mean(rt);
    stats.medrt = median(rt);
    stats.rtvar = std(rt);
    stats.hitrate = hr;
    stats.falarm = fa;
    stats.dprime = mdprime;
    stats.dprimec = mdprimec;
    stats.dprimei = mdprimei;
    stats.err = err;
    stats.errs = errs;
    stats.errCon = errCon;
    stats.errInCon = errInCon;
    % Extract the condition stats
    for k=1:8
        stats.(condnames{k}).meanrts = zeros(1,fn);
        stats.(condnames{k}).medrts = zeros(1,fn);
        stats.(condnames{k}).rtvars = zeros(1,fn);
        for c=1:fn
            if subj(c).conds.(condnames{k}).size > 0
                stats.(condnames{k}).meanrts(1,c) = subj(c).conds.(condnames{k}).meanrts;
                stats.(condnames{k}).medrts(1,c) =  subj(c).conds.(condnames{k}).medrts;
                stats.(condnames{k}).rtvars(1,c) =  subj(c).conds.(condnames{k}).rtvar;
            end
        end
        stats.(condnames{k}).meanrt = mean(stats.(condnames{k}).meanrts);
        stats.(condnames{k}).medrt = median(stats.(condnames{k}).medrts);
        stats.(condnames{k}).rtvar = std(stats.(condnames{k}).rtvars);
    end
    
end

end

function conds=findConditions( log, congruency, shape, hitmiss )
%% function success=findConditions( name, congruency, shape, hitmiss, prelat, lat )
%
%   Get condition data from saliency protocol log files based on given
%   params. Returns conds.trg for target trials and conds.rsp for responses
%   If no trials are found, conds returns with only its name, so
%   check for this with isfield()
%
    conds.name = [congruency '_' shape '_' hitmiss];
    idxtrg = strcmp('target', log.pritar_str) & ...
             strcmp(congruency, log.congruency_str) & ...
             strcmp(shape, log.shape_str) & ...
             strcmp(hitmiss, log.stim_type);
    idxrsp = circshift(idxtrg,1);
    if sum(idxtrg) > 0
        % skip subject (redundant)
        conds.trg.trial=log.trial(idxtrg);                      conds.rsp.trial=log.trial(idxrsp);
        % skip event_type (responses are discovered by right shift of trgs)
        conds.trg.code=log.code(idxtrg);                        conds.rsp.code =log.code(idxrsp);
        % skip coords_str (it is given in code)
        conds.trg.cond_num=log.cond_num(idxtrg);                conds.rsp.cond_num=log.cond_num(idxrsp);
        conds.trg.congruency_str=log.congruency_str(idxtrg);    conds.rsp.congruency_str=log.congruency_str(idxrsp);
        conds.trg.shape_str=log.shape_str(idxtrg);              conds.rsp.shape_str=log.shape_str(idxrsp);
        conds.trg.pritar_str=log.pritar_str(idxtrg);            conds.rsp.pritar_str=log.pritar_str(idxrsp);
        conds.trg.vertices_num=log.vertices_num(idxtrg);        conds.rsp.vertices_num=log.vertices_num(idxrsp);
        conds.trg.time=log.time(idxtrg);                        conds.rsp.time=log.time(idxrsp);
        conds.trg.uncertainty=log.uncertainty(idxtrg);          conds.rsp.uncertainty=log.uncertainty(idxrsp);
        conds.trg.duration=log.duration(idxtrg);                conds.rsp.duration=log.duration(idxrsp);
        conds.trg.stim_type=log.stim_type(idxtrg);              conds.rsp.stim_type=log.stim_type(idxrsp);
        % Get the RTs!
        conds.rt = conds.rsp.time - conds.trg.time;
        % Get stats
        conds.meanrts = mean(conds.rt);
        conds.medrts = median(conds.rt);
        conds.rtvar = std(conds.rt);
        conds.size = numel(conds.rt);
%         if strcmp('hit', hitmiss) == 1
%             conds.err = ;
%         end
    else
%         disp( 'No trials found under given conditions' );
        conds.size = 0;
    end
end


%     conds = cell(8,fn);
%     subj = struct('name', {}, 'conds', {}, 'rtvar', {}, 'err', {}, 'errCon', {}, 'errInCon', {}, 'errs', {}, 'hitrate', {}, 'falarm', {}, 'dprime', {});
%     hr = ones(1,fn);
%     fa = ones(1,fn);
%     % Loop once for every subject
%     for i = 1:fn
%         if nargin < 3,	[~, log{1,i}] = importPresLog(fullfile(path, files(i).name), true); end;
%         trg{i} = strcmp( 'target', log{1,i}.pritar_str );
%         rsp{i} = circshift(trg{i},1);
%         tlt(:,i) = log{1,i}.time(trg{i});
%         rlt(:,i) = log{1,i}.time(rsp{i});
%         rt(:,i) = rlt(:,i) - tlt(:,i);
%         % Get trials per condition - ALL CONDS ARRANGED ALPHABETICALLY
%             % Congruent non-shape hit condition
%         conds{1,i}=findConditions( log{1,i}, 'Congruent', 'nonShape', 'hit' );
%             % Congruent non-shape miss condition
%         conds{2,i}=findConditions( log{1,i}, 'Congruent', 'nonShape', 'miss' );
%             % Congruent shape hit condition
%         conds{3,i}=findConditions( log{1,i}, 'Congruent', 'shape', 'hit' );
%             % Congruent shape miss condition
%         conds{4,i}=findConditions( log{1,i}, 'Congruent', 'shape', 'miss' );
%             % INCongruent non-shape hit condition
%         conds{5,i}=findConditions( log{1,i}, 'InCon', 'nonShape', 'hit' );
%             % INCongruent non-shape miss condition
%         conds{6,i}=findConditions( log{1,i}, 'InCon', 'nonShape', 'miss' );
%             % INCongruent shape hit condition
%         conds{7,i}=findConditions( log{1,i}, 'InCon', 'shape', 'hit' );
%             % INCongruent shape miss condition
%         conds{8,i}=findConditions( log{1,i}, 'InCon', 'shape', 'miss' );
%         % Subject-level deets & stats
%         subj(i).name = log{1,i}.subject{1}(10:end);
%         subj(i).conds = conds{:,i};
%         subj(i).rtvar = mean([conds{1,i}.rtvar conds{3,i}.rtvar conds{5,i}.rtvar conds{7,i}.rtvar]);
%         % Errors are in proportionate form: as a percentage of ALL trials, not as a ratio of misses to hits.
%         subj(i).err = (conds{2,i}.size+conds{4,i}.size+conds{6,i}.size+conds{8,i}.size)*100/trials;
%         subj(i).errs = [(conds{2,i}.size*100)/(conds{1,i}.size+conds{2,i}.size)...
%                           (conds{4,i}.size*100)/(conds{3,i}.size+conds{4,i}.size)...
%                           (conds{6,i}.size*100)/(conds{5,i}.size+conds{6,i}.size)...
%                           (conds{8,i}.size*100)/(conds{7,i}.size+conds{8,i}.size)];
%         subj(i).errCon = (subj(i).errs(1)+subj(i).errs(2))/2;
%         subj(i).errInCon = (subj(i).errs(3)+subj(i).errs(4))/2;
%         hr(i) = ((conds{3,i}.size+conds{7,i}.size))/(conds{3,i}.size+conds{4,i}.size+conds{7,i}.size+conds{8,i}.size);
%         fa(i) = ((conds{2,i}.size+conds{6,i}.size))/(conds{1,i}.size+conds{2,i}.size+conds{5,i}.size+conds{6,i}.size);
%         subj(i).hitrate = hr(i);
%         subj(i).falarm = fa(i);