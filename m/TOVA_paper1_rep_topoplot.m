condition = {'response' 'inhibition'};

cEEGinb_H1 = pop_mergeset(cEEGinb_H1, 1:numel(cEEGinb_H1));
cEEGinb_H2 = pop_mergeset(cEEGinb_H2, 1:numel(cEEGinb_H2));
cEEGrsp_H1 = pop_mergeset(cEEGrsp_H1, 1:numel(cEEGrsp_H1));
cEEGrsp_H2 = pop_mergeset(cEEGrsp_H2, 1:numel(cEEGrsp_H2));

aEEGinb_H1 = pop_mergeset(aEEGinb_H1, 1:numel(aEEGinb_H1));
aEEGinb_H2 = pop_mergeset(aEEGinb_H2, 1:numel(aEEGinb_H2));
aEEGrsp_H1 = pop_mergeset(aEEGrsp_H1, 1:numel(aEEGrsp_H1));
aEEGrsp_H2 = pop_mergeset(aEEGrsp_H2, 1:numel(aEEGrsp_H2));


%% Topoplot the neg+pos peaks
% options for colour bar calibration
cax = {[-2 0] 'maxmin' 'absmax'};
% tstwn = [-1000 0; 0 1000];% in ms
tstwn = [150 250; 330 430];% in ms
% tstwn = [flip(-tstwn) tstwn];
rws = numel(tstwn) / 2;
cls = numel(which_g) * numel(which_c) * numel(which_h);
% itv = 565:616;

% set figure
fh = figure('Position', lbwh .* [1 1 1 0.4], 'Color', 'w');
plt = 0;

for s = 1%:numel(cax) %using only matched scales since they seem best
for c = which_c %conditions on Fig columns
for h = which_h
    
    % merge datasets
    cEEG = eval(['cEEG' cnd{c} '_H' h]);
    aEEG = eval(['aEEG' cnd{c} '_H' h]);
    cdnm = sprintf('%s,%s', condition{c}, ['H' h]);

    for t = 1:size(tstwn, 2) % time slices

        % get sample latencies from ms times
        t0 = find(cEEG.times > tstwn(t, 1), 1);
        t1 = find(cEEG.times < tstwn(t, 2), 1, 'last');
        sgt = sprintf('%d-%dms', tstwn(t, 1), tstwn(t, 2));

        % get shared colorbar values from joint data
        if s == 1
            dat = zeros(2, 128);
            for g = 1:2
                dat(g, :) =...
                    mean(mean(eval([grp{g}(1) 'EEG']).data(:, t0:t1, :), 3), 2);
            end
            scl = [min(dat, [], 'all') max(dat, [], 'all')];
        else
            scl = cax{s};
        end

        % plot groups across subplot columns
        for g = which_g

            subidx = plt + g + (cls ^ (t - 1) - (2 - t));
            subplot(rws, cls, subidx)
            
            dat = mean(mean(eval([grp{g}(1) 'EEG']).data(:, t0:t1, :), 3), 2);
            chlocs = eval([grp{g}(1) 'EEG']).chanlocs;
            
            topoplot(dat, chlocs...
                , 'electrodes', 'on'...
                , 'headrad', 0 ...
                , 'plotrad', 0.6 ...
                , 'emarker', {'.', 'k', 2, 1}...
                , 'emarker2', {[ROI{1, 1} ROI{1, 2}], 'o', 'k', 4, 1}...
                , 'maplimits', scl...
                , 'colormap', eval(clrmp{2}))
            hold on
            topoplot([], chlocs...
                , 'style', 'blank'...
                , 'headrad', 0.5 ...
                , 'plotrad', 0.6 ...
                , 'plotchans', [ROI{1, 1} ROI{1, 2}]...
                , 'electrodes', 'on'...
                , 'emarker', {'x', 'g', 3, 1}...
                , 'colormap', eval(clrmp{2}))
%             topoplot([], chlocs...
%                 , 'style', 'blank'...
%                 , 'headrad', 0.5 ...
%                 , 'plotrad', 0.6 ...
%                 , 'plotchans', ROI{1, 2}...
%                 , 'electrodes', 'on'...
%                 , 'emarker', {'x', 'm', 6, 2}...
%                 , 'colormap', eval(clrmp{2}))

            if mod(subidx, cls) == 1
                text(-0.6, -0.3, sprintf('%s', sgt)...
                    , 'FontWeight', 'bold', 'Rotation', 90, 'FontSize', 10)
            end
            if mod(subidx, 2) == 0
                pltpos = get(gca, 'Position');
                cb = colorbar('EastOutside', 'FontWeight', 'bold');
                title(cb, '\muV', 'FontSize', 9)
                cbpos = get(gca, 'Position');
                set(gca, 'Position', [cbpos(1:2) - [0.02 0] pltpos(3) cbpos(4)])
%                 set(cb, 'Position', cbpos - [0.02 0 0 0])
            end
            if subidx <= cls
                title(sprintf('%s', Grup{g}))
            else
                if mod(subidx, 2) == 1
                    text(0.1, 0.75, cdnm)
                end
            end
        end
%         plt = plt - 1;
    end
    plt = plt + 2;
end
end
end

svnm = ['TOVA_topoplot_' sprintf('%d:', tstwn) 'ms'];
% sgtitle(strrep(svnm, '_', ' '), 'Color', 'r')
print(fh, '-dsvg', fullfile(oud, 'topos', svnm))
%         export_fig(fullfile(oud, 'topos', svnm), '-png')

%%
close
% clear fh grp cnd cax plt itv s m c g dat cb
