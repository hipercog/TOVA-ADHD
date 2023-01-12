%% TEST RANDOM DRAWS

fs = dir(fullfile(oud, 'spectra', 'tova_spectra_*.mat'));

testdat = cell(1, numel(fs));

for i = 1:numel(fs)
%     testdat.(sprintf('draw%d', i)) = load(fs(i).name);
    testdat{i} = load(fs(i).name);
end

fnms = fieldnames(testdat{1}.specdat);

for i = 1:numel(fnms)
    tmp = cellfun(@(x) x.specdat.(fnms{i}), testdat, 'Un', 0);
    tmp = cellfun(@(x) mean(x(:, 4:16))', tmp, 'Un', 0);
    vars.(fnms{i}) = mean(var([tmp{:}], 0, 2));
end

varmat = struct2mat(vars);
fprintf('BL variance: %1.4f\n', mean(varmat(1:2:32)))
fprintf('PS variance: %1.4f\n', mean(varmat(2:2:32)))