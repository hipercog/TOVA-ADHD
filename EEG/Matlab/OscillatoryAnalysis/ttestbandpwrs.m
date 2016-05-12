function test=ttestbandpwrs( eeg, roi1, roi2, iep )
%% FUNCTION ttestbandpwrs
% PURPOSE:  Test the difference between given ROIs in given EEG for power in given bands.
% PARAMS:   'eeg' - eeglab data structure
%           'roi1' - channels from first region of interest to compare, row order.
%           'roi2' - channels from second region of interest to compare, row order.
%           'iep' - band division points (in Hz) for individualisation
% USAGE:    
% Note:     
% CALLS:    
    nbands = numel(iep)-1;
    bdpw1 = zeros(numel(roi1), nbands);
    bdpw2 = zeros(numel(roi2), nbands);
%     test = zeros(1,nbands);
    
    % Limit the data to specific channels and trim excess data 
    data1 = eeg.data(roi1,:);
    seg = floor(size(data1,2)/(5*eeg.srate));
    data1 = data1(:,1:seg*5*eeg.srate);
    
    data2 = eeg.data(roi2,:);
    seg = floor(size(data2,2)/(5*eeg.srate));
    data2 = data2(:,1:seg*5*eeg.srate);

    % Design 4th order IIR (butterworth) bandpass filters according to the IAF-corrected frequency bands
    Fs = eeg.srate;
    N = 4;
    for p = 1:nbands
        Hd = design( fdesign.bandpass('N,F3dB1,F3dB2', N, iep(p), iep(p+1), Fs), 'butter' );
        for roidx = 1:numel(roi1)
            Fd1 = filter(Hd, data1(roidx,:)')';
            Fd2 = filter(Hd, data2(roidx,:)')';
            bdpw1(roidx,p) = median( log( median( reshape( mean(Fd1.^2,1), seg, 5*eeg.srate ), 2 )+1 ) );
            bdpw2(roidx,p) = median( log( median( reshape( mean(Fd2.^2,1), seg, 5*eeg.srate ), 2 )+1 ) );
        end
%         [test(1,p) test(2,p)] = ttest( bdpw1(:,p)', bdpw2(:,p)' );
    end
%     [test(1,:) test(2,:)]
    test = ttest( bdpw1, bdpw2 );
end