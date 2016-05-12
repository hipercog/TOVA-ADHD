function [answ sliceAns]=statCheckEpochs( redo, dirc, cndn, trodes )

if nargin < 1,  redo=false; end;
if nargin < 2,  dirc=pwd;   end;
if nargin < 3
    cndn={'Congr shape' 'Congr nonShape' 'InCon shape' 'InCon nonShape'};
end
if nargin < 4
    % For the classic range of 10-20 trodes: 'O1' 'O2' 'P3' 'P4' 'C3' 'C4' 'F3' 'F4'
    trodes={'A15' 'A28' 'A7' 'B4' 'D19' 'B22' 'D4' 'C4'};
end
% Parameter values
cn = numel(cndn);
tn = numel(trodes);
% these are hardcoded right now - could fix by picking epoch length, srate,
% from first .set file in given dirc; maybe pass BL remove length?
datl = 440;     % number of data points in my epochs right now
os = 51;    % number of data points in 100ms @ 512Hz
sliceL = 26;
% Check what to do
if redo
    % make space for avg curves
    pavgs = cell(cn,tn);
    cavgs = cell(cn,tn);
    % go through conditions at the top level
    for c=1:cn
        trgstr = [cndn{c} ' hit'];
        aggr=zeros(tn,datl-2*os);
        % Load all patient eegs (in this cond) and get their average data.
        Ptnts = dir( fullfile(dirc, ['1*' trgstr '.set']) );
        nump = numel(Ptnts);
        for p=1:nump
            EEG = pop_loadset('filename', Ptnts(p).name);
            for e=1:tn
                idx=strcmp(trodes{e}, {EEG.chanlocs.labels});
                % get desired data points, all epochs, baseline remove 1:100ms
                t1 = mean(EEG.data(idx,1:os,:),2);
                temp = EEG.data(idx,os+1:datl-os,:)-t1(3);
                aggr(e,:)=aggr(e,:)+mean(temp, 3);
            end
        end
        for e=1:tn, pavgs{c,e} = aggr(e,:)./nump;	end
        % Then do the same for Ctrls
        Ctrls = dir( fullfile(dirc, ['2*' trgstr '.set']) );
        numc = numel(Ctrls);
        for t=1:numc
            EEG = pop_loadset('filename', Ctrls(t).name);
            for e=1:tn
                idx=strcmp(trodes{e}, {EEG.chanlocs.labels});
                % get desired data points, all epochs, baseline remove 1:100ms
                t1 = mean(EEG.data(idx,1:os,:),2);
                temp = EEG.data(idx,os+1:datl-os,:)-t1(3);
                aggr(e,:)=aggr(e,:)+mean(temp, 3);            
            end
        end
        for e=1:tn, cavgs{c,e} = aggr(e,:)./numc;	end
    end
    % Save out the just gathered data...
    save('statCheckable.mat', 'pavgs', 'cavgs');
else
    if exist('statCheckable.mat', 'file') == 2
        load('statCheckable.mat');
    else
        disp('No previous values loaded: try statCheckEpochs( redo=true )');
    end
end
% Then test each condition x trode with statcond()
% answ=zeros(cn,tn,4);
answ=struct('Tval', {}, 'df', {}, 'pVal', {});
sliceAns=struct('Tvals', {}, 'dfs', {}, 'pVals', {});
for c=1:cn
    for e=1:tn
        [answ(c,e).Tval answ(c,e).df answ(c,e).pVal] = ...
            statcond( {pavgs{c,e} cavgs{c,e}}, 'mode', 'bootstrap' );
        % And do the 50ms slices of this trode
        pslice = mat2tiles(pavgs{c,e},1,sliceL);
        cslice = mat2tiles(cavgs{c,e},1,sliceL);
        for s=1:numel(pslice)
            [sliceAns(c,e,s).Tvals sliceAns(c,e,s).dfs sliceAns(c,e,s).pVals] = ...
                statcond( {pslice{s} cslice{s}}, 'mode', 'bootstrap' );            
        end
    end
end
pVals=zeros(4,8,13);
Tvals=zeros(4,8,13);
dfs=zeros(4,8,13);
for s=1:13
    pVals(:,:,s)=[[sliceAns(1,:,s).pVals]; [sliceAns(2,:,s).pVals]; [sliceAns(3,:,s).pVals]; [sliceAns(4,:,s).pVals]];
    Tvals(:,:,s)=[[sliceAns(1,:,s).Tvals]; [sliceAns(2,:,s).Tvals]; [sliceAns(3,:,s).Tvals]; [sliceAns(4,:,s).Tvals]];
    dfs(:,:,s)=[[sliceAns(1,:,s).dfs]; [sliceAns(2,:,s).dfs]; [sliceAns(3,:,s).dfs]; [sliceAns(4,:,s).dfs]];
end
save('b00tTslices', 'sliceAns', 'pVals', 'Tvals', 'dfs');

% Save the global stat answers
pVal=[answ(1,:).pVal; answ(2,:).pVal; answ(3,:).pVal; answ(4,:).pVal]; %#ok<NASGU>
Tval=[answ(1,:).Tval; answ(2,:).Tval; answ(3,:).Tval; answ(4,:).Tval]; %#ok<NASGU>
df=[answ(1,:).df; answ(2,:).df; answ(3,:).df; answ(4,:).df]; %#ok<NASGU>
save('b00tTepochs', 'answ', 'pVal', 'Tval', 'df');
