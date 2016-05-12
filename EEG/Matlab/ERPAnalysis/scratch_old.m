direc = pwd;
outd = fullfile( direc, 'Problem with chanlocs' );
if ~isdir(outd),    mkdir(outd);    end;
files = dir( fullfile(direc, '20*.set') );
fn = numel(files);

goldifile = '2021C_InCon shape hit.set';
EEG = pop_loadset( 'filename', goldifile, 'filepath', direc);
goldilocs = EEG.chanlocs;

for i = 1:fn
    [~, savename, ~] = fileparts(files(i).name);
    EEG = pop_loadset( 'filename', [savename '.set'], 'filepath', direc);
    if strcmp( EEG.chanlocs(1).labels, 'A1' ) == 0
        EEG.chanlocs = goldilocs;
        % Save the re-chanloc'd .set
        pop_saveset( EEG, 'filename', [savename '.set'], 'filepath', outd );
    end
end

% direc = pwd;
% outd = fullfile( direc, '..', 'EVEN' );
% if ~isdir(outd),    mkdir(outd);    end;
% files = dir( fullfile(direc, '*.set') );
% fn = numel(files);
% lowest = 487;
% times = [-0.5996 0.3496];
% for i = 1:fn
%     [~, savename, ~] = fileparts(files(i).name);
%     EEG = pop_loadset( 'filename', [savename '.set'], 'filepath', direc);
% %     for k=1:EEG.trials
% %         if ~cellfun(@isempty, strfind(EEG.epoch(1,k).eventtype(:), 'target')) ==...
% %            cell2mat(EEG.epoch(1,k).eventlatency(:))==0
% %             idx = find(~cellfun(@isempty, strfind(EEG.epoch(1,k).eventtype(:), 'target')));
% %             EEG.epoch(1,k).eventtype{idx} = 'target';
% %         end
% %     end
%     % Trim all epochs down to 487 points
%     EEG = pop_select(EEG, 'time', times);
%     % Save this slimmed down .set
%     pop_saveset( EEG, 'filename', [savename '.set'], 'filepath', outd );
% end

% direc = pwd;
% outd = fullfile( direc, '..', 'REpoc' );
% if ~isdir(outd),    mkdir(outd);    end;
% files = dir( fullfile(direc, '*.set') );
% fn = numel(files);
% for i = 1:fn
%     [~, savename, ~] = fileparts(files(i).name);
%     EEG = pop_loadset( 'filename', [savename '.set'], 'filepath', direc);
%     % Find index of correct-trial epochs
%     idx = zeros(1,EEG.trials);
%     if EEG.trials > 1
%         for k=1:EEG.trials
%             if ~cellfun(@isempty, strfind(EEG.epoch(1,k).eventtype(:), 'target')) ==...
%                cell2mat(EEG.epoch(1,k).eventlatency(:))==0
%                 idx(k) = 1;
%             end
%         end
%     end
%     % Select only those epochs
%     if sum(idx) > 0
%         EEG = pop_select(EEG, 'trial', find(idx));
%         % Save this slimmed down .set
%         pop_saveset( EEG, 'filename', [savename '.set'], 'filepath', outd );
%     else
%         disp(['No valid trials found in: ' savename]);
%     end
% end