%% ERPIM x pre-stim alpha phase

% Old data prep
% cEEGh1{1} = pop_loadset('filepath', ind, 'filename', 'MERGED_H1_CORRECT_RANDOMISED_CONTROL.set');
% cEEGh2{1} = pop_loadset('filepath', ind, 'filename', 'MERGED_H2_CORRECT_RANDOMISED_CONTROL.set');
% aEEGh1{1} = pop_loadset('filepath', ind, 'filename', 'merged_h1_correct_randomised_adhd.set');
% aEEGh2{1} = pop_loadset('filepath', ind, 'filename', 'merged_h2_correct_randomised_adhd.set');
% 
% % sub-set epoch
% times = [-200 400];
% lockev = {'TRGT' 'NONT'};
% lockevent = {'Target' 'NonTarget'};
% for l = 1:numel(lockev)
%     eeg = pop_epoch(cEEGh1{1}, lockev(l), times / 1000);
%     sub = shuffle(1:eeg.trials);
%     cEEGh1{l + 1} = pop_select(eeg, 'trial', sub(1:457)); %#ok<*SAGROW>
%     eeg = pop_epoch(cEEGh2{1}, lockev(l), times / 1000);
%     sub = shuffle(1:eeg.trials);
%     cEEGh2{l + 1} = pop_select(eeg, 'trial', sub(1:457));
%     
%     eeg = pop_epoch(aEEGh1{1}, lockev(l), times / 1000);
%     sub = shuffle(1:eeg.trials);
%     aEEGh1{l + 1} = pop_select(eeg, 'trial', sub(1:862));
%     eeg = pop_epoch(aEEGh2{1}, lockev(l), times / 1000);
%     sub = shuffle(1:eeg.trials);
%     aEEGh2{l + 1} = pop_select(eeg, 'trial', sub(1:862));
% end


%% INIT
times = [-200 400];
lockev = {'TRGT' 'NONT'};
lockevent = {'Target' 'NonTarget'};
cax = [-11 11];
ylm = [-7 7];
ylab = 'Phase-sorted Trials';
fh = struct(ROI{2,1}...
    , figure('Units', 'pixels', 'Posi', lbwh .* [1 1 0.75 0.66], 'Color', 'w')...
    , ROI{2,2}...
    , figure('Units', 'pixels', 'Posi', lbwh .* [1 1 0.75 0.66], 'Color', 'w'));
plt = 1;
erpN = 2;
clrs = {[1 0 0 0.7]; [0.8 0.2 0.1 1]};
xtl = '';


%% LOOP to create phase-locked ERPIMs
for g = which_g
for c = which_c
for h = which_h

    eeg = eval([grp{g}(1) 'EEG' cnd{c} '_H' h]);
    eeg = pop_mergeset(eeg, 1:numel(eeg));
    sub = shuffle(1:eeg.trials);
    eeg = pop_select(eeg, 'trial', sub(1:eeg.trials / g));
    eeg = pop_epoch(eeg, lockev(c), times / 1000);
    sr = eeg.srate;
    pts = eeg.pnts;
    smth = eeg.trials / 30;

    for r = which_r
        if plt == 4
            cbr = 'on';
        else 
            cbr = 'off';
        end
        if plt < 5
            ttl = ['H' h ' - ' lockevent{c}];
        else
            ttl = '';
        end
        figure(fh.(ROI{2,r}))
        subplot(2, 4, plt)
        figure(fh.(ROI{2,r}))
        
        [~, ~, ~, ~, axh, ~] =...
            erpimage( mean(eeg.data(ROI{1,r}, :, :), 1) ...
                , [] ...
                , [times(1) pts sr] ...
                , ''...
                , smth...
                , 1 ...
                , 'avg_type', 'Gaussian'...
                , 'erp', erpN ...
                , 'yerplabel', ''...
                , 'caxis', cax...
                , 'cbar', cbr...
                , 'cbar_title', '\muV'...
                , 'phasesort', [-80 5 8 12 0]...
                , 'noxlabel', 'on');
%                 , 'erpalpha', 0.001 ...

        % ERP underplot tweaking
        title(ttl)
        if plt > 4
            xtickangle(axh{3}, 45)
            set(axh{1}, 'Position', axh{1}.Position + [0 0.06 0 0])
            set(axh{3}, 'Position', axh{3}.Position + [0 0.06 0 0])
            xtl = times(1):100:times(2);
            clrs = {[0 0 1 0.7]; [0.1 0.2 0.8 1]};
        end
        if erpN == 2
            set(axh{3}.Legend, 'visible','off')
            set(axh{3}.Children(2:3), 'LineWidth', 1.4, {'LineStyle'}, {'-.'; '-'}, {'Color'}, clrs)
        else
            set(axh{3}.Children(2), 'Color', clrs{erpN})
            set(axh{3}.Children(5), 'FaceColor', [0.8 0.8 0.8])
        end
        set(axh{3}, 'Color', 'w'...
                  , 'XTick', times(1):100:times(2)...
                  , 'XTickLabel', xtl...
                  , 'FontWeight', 'bold'...
                  , 'YLim', ylm...
                  , 'YAxisLocation', 'right'...
                  , 'YTick', [ylm(1) 0 ylm(2)]...
                  , 'YTickLabel', {num2str(ylm(1)) '0' num2str(ylm(2))})
        if any(plt == [4 8])
            ylabel(axh{3}, '\muV'...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [1.15 1]...
                , 'Rotation', 0)
        end
        if plt == 8
            xlabel(axh{3}, sprintf('Time\n(ms)')...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [0.9 -0.5]...
                , 'Rotation', 45)
        end
        % color tweaking
        if plt == 4
            cbrh = axh{2}.Position(4) / 5;
            axh{2}.Position(3:4) = axh{2}.Position(3:4) .* 0.8;
            axh{2}.Position(2) = axh{2}.Position(2) + cbrh;
        end
        % ERPIM plot tweaking
        if any(plt == [1 5])
            lg = legend(axh{3}.Children(2:3)...
               , {sprintf('phase >\nmedian') sprintf('phase <\nmedian')}...
               , 'Location', 'SouthWest'...
               , 'FontWeight', 'normal'...
               , 'Box', 'off');
            lg.Position = lg.Position - [0.07 0 0 0];
            lg.ItemTokenSize = [18 3];
            
            ylabel(axh{1}, sprintf('%s\n%s', upper(Grup{g}), ylab)...
                , 'FontWeight', 'bold'...
                , 'Units', 'normalized'...
                , 'Position', [-0.2 0.5]...
                , 'Rotation', 90)
        else
            ylabel(axh{1}, '')
        end

    end

    plt = plt + 1;
end
end
end

%%
figure(fh.(ROI{2,1}))
colormap(cmap)
print(fh.(ROI{2,1}), '-dsvg', fullfile(oud, 'erpimXphase', ['TOVA-erpimXphase-' ROI{2,1}]))
figure(fh.(ROI{2,2}))
colormap(cmap)
print(fh.(ROI{2,2}), '-dsvg', fullfile(oud, 'erpimXphase', ['TOVA-erpimXphase-' ROI{2,2}]))

%%
close all
clear times locke* cax cbr* plt ylab g c h r eeg sr pts fh axh ttl ylm sub